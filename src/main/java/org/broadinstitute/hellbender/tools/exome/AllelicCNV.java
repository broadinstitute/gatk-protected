package org.broadinstitute.hellbender.tools.exome;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionData;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionInitializer;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionState;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Detects copy-number events using allelic-count data and GATK CNV output.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Detect copy-number events in a tumor sample using allelic-count data and GATK CNV output. " +
                "Allelic-count data (reference/alternate counts from the GetHetCoverage tool) is segmented using " +
                "circular binary segmentation; the result is combined with the target coverages " +
                "and called segments found by the GATK CNV tool. Bayesian parameter estimation of models for the " +
                "copy ratios and minor allele fractions in each segment is performed using Markov chain Monte Carlo.",
        oneLineSummary = "Detect copy-number events using allelic-count data and GATK CNV output.",
        programGroup = CopyNumberProgramGroup.class
)
public class AllelicCNV extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1l;

    private static final double INITIAL_ALLELIC_BIAS_GUESS = 1.;
    private static final double SNP_EVENT_MAF_THRESHOLD = 0.475;

    //filename tags for output
    protected static final String SNP_MAF_SEG_FILE_TAG = "MAF";
    protected static final String FINAL_SEG_FILE_TAG = "final";
    protected static final String GATK_SEG_FILE_TAG = "cnv";
    protected static final String CGA_ACS_SEG_FILE_TAG = "acs";

    //CLI arguments
    protected static final String OUTPUT_PREFIX_LONG_NAME = "outputPrefix";
    protected static final String OUTPUT_PREFIX_SHORT_NAME = "pre";

    protected static final String NUM_SAMPLES_COPY_RATIO_LONG_NAME = "numSamplesCopyRatio";
    protected static final String NUM_SAMPLES_COPY_RATIO_SHORT_NAME = "numSampCR";

    protected static final String NUM_BURN_IN_COPY_RATIO_LONG_NAME = "numBurnInCopyRatio";
    protected static final String NUM_BURN_IN_COPY_RATIO_SHORT_NAME = "numBurnCR";

    protected static final String NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME = "numSamplesAlleleFraction";
    protected static final String NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME = "numSampAF";

    protected static final String NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME = "numBurnInAlleleFraction";
    protected static final String NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME = "numBurnAF";

    protected static final String MAX_NUM_SNP_SEGMENTATION_ITERATIONS_LONG_NAME = "maxNumIterationsSNPSeg";
    protected static final String MAX_NUM_SNP_SEGMENTATION_ITERATIONS_SHORT_NAME = "maxIterSNP";

    protected static final String USE_ALL_COPY_RATIO_SEGMENTS_LONG_NAME = "useAllCopyRatioSegments";
    protected static final String USE_ALL_COPY_RATIO_SEGMENTS_SHORT_NAME = "useAllCRSeg";

    @Argument(
            doc = "Input file for tumor-sample ref/alt read counts at normal-sample heterozygous-SNP sites (output of GetHetCoverage tool).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpCountsFile;

    @Argument(
            doc = "Input file for tumor-sample tangent-normalized target coverages (.tn.tsv output of GATK CNV tool).",
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File targetCoveragesFile;

    @Argument(
            doc = "Input file for tumor-sample target-coverage segments with calls (output of GATK CNV tool).",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File targetSegmentsFile;

    @Argument(
            doc = "Input file for allelic-bias panel of normals.",
            fullName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_SHORT_NAME,
            optional = true
    )
    protected File allelicPONFile;

    @Argument(
            doc = "Prefix for output files. Will also be used as the sample name if that is not provided." +
                    "(Note: if this is a file path or contains slashes (/), " +
                    "the string after the final slash will be used as the sample name if that is not provided.)",
            fullName = OUTPUT_PREFIX_LONG_NAME,
            shortName = OUTPUT_PREFIX_SHORT_NAME,
            optional = false
    )
    protected String outputPrefix;

    @Argument(
            doc = "Total number of MCMC samples for copy-ratio model.",
            fullName = NUM_SAMPLES_COPY_RATIO_LONG_NAME,
            shortName = NUM_SAMPLES_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    protected int numSamplesCopyRatio = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for copy-ratio model.",
            fullName = NUM_BURN_IN_COPY_RATIO_LONG_NAME,
            shortName = NUM_BURN_IN_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    protected int numBurnInCopyRatio = 50;

    @Argument(
            doc = "Total number of MCMC samples for allele-fraction model.",
            fullName = NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected int numSamplesAlleleFraction = 200;

    @Argument(
            doc = "Number of burn-in samples to discard for allele-fraction model.",
            fullName = NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected int numBurnInAlleleFraction = 100;

    @Argument(
            doc = "Maximum number of iterations allowed for SNP segmentation.",
            fullName = MAX_NUM_SNP_SEGMENTATION_ITERATIONS_LONG_NAME,
            shortName = MAX_NUM_SNP_SEGMENTATION_ITERATIONS_SHORT_NAME,
            optional = true
    )
    protected int maxNumSNPSegmentationIterations = 25;

    @Argument(
            doc = "Enable use of all copy-ratio--segment breakpoints. " +
                    "(Default behavior uses only breakpoints from segments not called copy neutral, " +
                    "if calls are available in output of GATK CNV provided, and none otherwise.)",
            fullName = USE_ALL_COPY_RATIO_SEGMENTS_LONG_NAME,
            shortName = USE_ALL_COPY_RATIO_SEGMENTS_SHORT_NAME,
            optional = true
    )
    protected boolean useAllCopyRatioSegments = false;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        final String originalLogLevel =
                (ctx.getLocalProperty("logLevel") != null) ? ctx.getLocalProperty("logLevel") : "INFO";
        ctx.setLogLevel("WARN");

        //the string after the final slash in the output prefix (which may be an absolute file path) will be used as the sample name
        final String sampleName = outputPrefix.substring(outputPrefix.lastIndexOf("/") + 1);

        logger.info("Starting workflow for " + sampleName + "...");

        //make Genome from input target coverages and SNP counts
        logger.info("Loading input files...");
        final Genome genome = new Genome(targetCoveragesFile, snpCountsFile, sampleName);

        //load allelic-bias panel of normals if provided
        final AllelicPanelOfNormals allelicPON =
                allelicPONFile != null ? new AllelicPanelOfNormals(allelicPONFile) : AllelicPanelOfNormals.EMPTY_PON;

        //load called target-coverage segments from input file
        final List<ModeledSegment> targetSegmentsWithCalls =
                SegmentUtils.readModeledSegmentsFromSegmentFile(targetSegmentsFile);
        logger.info("Number of target-coverage segments from CNV output: " + targetSegmentsWithCalls.size());

        //merge copy-neutral and uncalled segments (unless disabled) and fix up target-segment start breakpoints
        logger.info("Preparing target-coverage segments...");
        final List<SimpleInterval> targetSegments = prepareTargetSegments(genome, targetSegmentsWithCalls);

        //segment SNPs using CBS on per-SNP MLE minor allele fraction iteratively until convergence
        //(final segment file is output as a side effect)
        logger.info("Performing SNP segmentation...");
        final String snpSegmentsOutputPrefix = outputPrefix + "-" + SNP_MAF_SEG_FILE_TAG;
        final List<SimpleInterval> snpSegments = performSNPSegmentationStep(sampleName, genome, snpSegmentsOutputPrefix);

        //call SNP segments generated by previous step
        final List<ModeledSegment> snpSegmentsWithCalls = callSNPSegments(snpSegmentsOutputPrefix);

        logger.info("Unioning called events...");
        //add breakpoints from called AF events to target-coverage segments
        final List<SimpleInterval> targetSegmentsWithAFEvents = unionTargetSegmentsWithAFEvents(targetSegments, genome, snpSegmentsWithCalls);
        logger.info("Number of target-coverage segments after unioning SNP events: " + targetSegmentsWithAFEvents.size());
        //add breakpoints from called CR events to SNP segments
        final List<SimpleInterval> snpSegmentsWithCREvents = unionSNPSegmentsWithCREvents(snpSegments, genome, targetSegmentsWithCalls);
        logger.info("Number of SNP segments after unioning target-coverage events: " + snpSegmentsWithCREvents.size());

        final SegmentedGenome targetSegmentedGenome = new SegmentedGenome(targetSegmentsWithAFEvents, genome);
        final SegmentedGenome snpSegmentedGenome = new SegmentedGenome(snpSegmentsWithCREvents, genome);

        //initial MCMC model fitting performed by ACNVModeller constructor
        final ACNVModeller modeller = new ACNVModeller(targetSegmentedGenome, snpSegmentedGenome, allelicPON,
                numSamplesCopyRatio, numBurnInCopyRatio, numSamplesAlleleFraction, numBurnInAlleleFraction, ctx);

        //write final model fit to file
        final File finalModeledSegmentsFile = new File(outputPrefix + "-" + FINAL_SEG_FILE_TAG + ".seg");
        modeller.writeACNVModeledSegmentFile(finalModeledSegmentsFile);

        //write file for GATK CNV formatted seg file
        final File finalModeledSegmentsFileAsGatkCNV = new File(outputPrefix + "-" + FINAL_SEG_FILE_TAG + "." + GATK_SEG_FILE_TAG + ".seg");
        modeller.writeModeledSegmentFile(finalModeledSegmentsFileAsGatkCNV);

        //write file for ACS-compatible output to help Broad CGA
        final File finalACSModeledSegmentsFile = new File(outputPrefix + "-" + FINAL_SEG_FILE_TAG + "." + CGA_ACS_SEG_FILE_TAG + ".seg");
        modeller.writeAllelicCapSegFile(finalACSModeledSegmentsFile);

        ctx.setLogLevel(originalLogLevel);
        logger.info("SUCCESS: Allelic CNV run complete for sample " + sampleName + ".");
    }

    //merge copy-neutral and uncalled segments (unless disabled) and fix up target-segment start breakpoints
    private List<SimpleInterval> prepareTargetSegments(final Genome genome,
                                                       final List<ModeledSegment> targetSegmentsWithCalls) {
        final List<SimpleInterval> targetSegmentsWithUnfixedStarts;
        if (!useAllCopyRatioSegments) {
            logger.info("Merging copy-neutral and uncalled target-coverage segments...");
            targetSegmentsWithUnfixedStarts = SegmentMergeUtils.mergeNeutralSegments(targetSegmentsWithCalls);
            logger.info("Number of segments after copy-neutral merging: " + targetSegmentsWithUnfixedStarts.size());
        } else {
            targetSegmentsWithUnfixedStarts = targetSegmentsWithCalls.stream().map(ModeledSegment::getSimpleInterval)
                    .collect(Collectors.toList());
        }
        //fix up target-segment start breakpoints (convert from target-end--target-end to target-start--target-end)
        return SegmentUtils.fixTargetSegmentStarts(targetSegmentsWithUnfixedStarts, genome.getTargets());
    }

    //segment SNPs using CBS on per-SNP MLE minor allele fraction iteratively until convergence
    //(final segment file is output as a side effect)
    private List<SimpleInterval> performSNPSegmentationStep(final String sampleName, final Genome genome, final String snpSegmentsOutputPrefix) {
        if (maxNumSNPSegmentationIterations <= 0) {
            throw new UserException.BadInput("Maximum number of SNP-segmentation iterations must be positive.");
        }
        //initial segmentation of SNPs on per-SNP MLE minor allele fraction, assuming no allelic bias
        //(segments are written to temporary file)
        logger.info("Performing initial SNP segmentation (assuming no allelic bias)...");
        final File snpSegmentFile;
        try {
            snpSegmentFile = File.createTempFile(snpSegmentsOutputPrefix + "-tmp", ".seg");
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not create temporary SNP segmentation file.", e);
        }
        List<SimpleInterval> snpSegments =
                calculateAndWriteSNPSegmentation(sampleName, genome, snpSegmentFile, INITIAL_ALLELIC_BIAS_GUESS);
        logger.info("Initial number of SNP segments: " + snpSegments.size());
        final Set<List<SimpleInterval>> snpSegmentationsFound = new HashSet<>();
        //perform iterations of SNP segmentation until convergence
        for (int numIterations = 1; numIterations <= maxNumSNPSegmentationIterations; numIterations++) {
            logger.info("SNP-segmentation iteration: " + numIterations);
            snpSegmentationsFound.add(snpSegments);
            //use AlleleFractionInitializer to determine mean of global allelic-bias distribution given current SNP segmentation
            final double allelicBias = calculateAllelicBias(genome, snpSegmentFile);
            logger.info(String.format("Mean allelic bias: %.3f", allelicBias));
            //resegment SNPs on per-SNP MLE minor allele fraction assuming mean allelic bias found by AlleleFractionInitializer
            //(segments are written to temporary file)
            snpSegments = calculateAndWriteSNPSegmentation(sampleName, genome, snpSegmentFile, allelicBias);
            logger.info("Number of SNP segments: " + snpSegments.size());
            //stop if a previously found SNP segmentation is found again
            if (snpSegmentationsFound.contains(snpSegments)) {
                //take final SNP segmentation to be the same as the last one found
                logger.info("Final number of SNP segments: " + snpSegments.size());
                final File finalSNPSegmentFile = new File(snpSegmentsOutputPrefix + ".seg");
                snpSegments = calculateAndWriteSNPSegmentation(sampleName, genome, finalSNPSegmentFile, allelicBias);
                break;
            }
        }
        return snpSegments;
    }

    //return SNP segments from CBS (writes file as a side effect)
    private List<SimpleInterval> calculateAndWriteSNPSegmentation(final String sampleName, final Genome genome,
                                                                  final File snpSegmentFile, final double allelicBias) {
        SNPSegmenter.writeSegmentFile(genome.getSNPs(), sampleName, snpSegmentFile, allelicBias);
        return SegmentUtils.readIntervalsFromSegmentFile(snpSegmentFile);
    }

    //use AlleleFractionInitializer to determine mean of global allelic-bias distribution given SNP segmentation file
    private double calculateAllelicBias(final Genome genome, final File snpSegmentFile) {
        final SegmentedGenome segmentedGenomeForSNPSegmentation = new SegmentedGenome(snpSegmentFile, genome);
        final AlleleFractionData data = new AlleleFractionData(segmentedGenomeForSNPSegmentation);
        final AlleleFractionState initialState = new AlleleFractionInitializer(data).getInitializedState();
        return initialState.meanBias();
    }

    //call SNP segments generated by previous step
    private List<ModeledSegment> callSNPSegments(final String snpSegmentsOutputPrefix) {
        final List<ModeledSegment> snpSegmentsWithMeans = SegmentUtils.readModeledSegmentsFromSegmentFile(new File(snpSegmentsOutputPrefix + ".seg"));
        return snpSegmentsWithMeans.stream().map(this::callSNPSegment).collect(Collectors.toList());
    }

    private ModeledSegment callSNPSegment(final ModeledSegment snpSegmentWithMean) {
        if (snpSegmentWithMean.getSegmentMean() <= SNP_EVENT_MAF_THRESHOLD) {
            return new ModeledSegment(snpSegmentWithMean.getSimpleInterval(), ReCapSegCaller.AMPLIFICATION_CALL,
                    snpSegmentWithMean.getTargetCount(), ParamUtils.log2(snpSegmentWithMean.getSegmentMean()));
        }
        return new ModeledSegment(snpSegmentWithMean.getSimpleInterval(), ReCapSegCaller.NEUTRAL_CALL,
                snpSegmentWithMean.getTargetCount(), ParamUtils.log2(snpSegmentWithMean.getSegmentMean()));
    }

    private List<SimpleInterval> unionTargetSegmentsWithAFEvents(final List<SimpleInterval> targetSegments,
                                                                 final Genome genome,
                                                                 final List<ModeledSegment> snpSegmentsWithCalls) {
        final List<SimpleInterval> afEvents =
                snpSegmentsWithCalls.stream().filter(s -> s.getCall().equals(CnvCaller.AMPLIFICATION_CALL)).map(ModeledSegment::getSimpleInterval).collect(Collectors.toList());
        return SegmentUtils.unionSegmentsNaively(targetSegments, afEvents, genome);
    }

    private List<SimpleInterval> unionSNPSegmentsWithCREvents(final List<SimpleInterval> snpSegments,
                                                                 final Genome genome,
                                                                 final List<ModeledSegment> targetSegmentsWithCalls) {
        final List<SimpleInterval> crEventsWithUnfixedStarts =
                targetSegmentsWithCalls.stream().filter(s -> s.getCall().equals(CnvCaller.AMPLIFICATION_CALL)).map(ModeledSegment::getSimpleInterval).collect(Collectors.toList());
        final List<SimpleInterval> crEvents = SegmentUtils.fixTargetSegmentStarts(crEventsWithUnfixedStarts, genome.getTargets());
        return SegmentUtils.unionSegmentsNaively(snpSegments, crEvents, genome);
    }
}