package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Iterables;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Decile;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Detects copy-number events using allelic-count data and GATK CNV output.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Model tumor heterogeneity using a Dirichlet mixture of subclones with copy number variation.",
        oneLineSummary = "Model tumor heterogeneity using a Dirichlet mixture of subclones with copy number variation",
        programGroup = CopyNumberProgramGroup.class
)
public class TumorHeterogeneity extends SparkCommandLineProgram {
    private static final long serialVersionUID = 19738246L;

    private static final long RANDOM_SEED = 13;
    private static final double CREDIBLE_INTERVAL_ALPHA = 0.95;
    private static final int NUM_POPULATIONS_CLONAL = 2;
    private static final PloidyState NORMAL_PLOIDY_STATE = new PloidyState(1, 1);

    //filename tags for output
    protected static final String CLONAL_RESULT_FILE_SUFFIX = ".th.clonal.tsv";
    protected static final String RESULT_FILE_SUFFIX = ".th.tsv";
    protected static final String FILTERED_SEGMENTS_FILE_SUFFIX = ".filtered.seg";
    protected static final String MAF_CR_FILE_SUFFIX = ".maf-cr.tsv";

    //CLI arguments
    protected static final String OUTPUT_PREFIX_LONG_NAME = "outputPrefix";
    protected static final String OUTPUT_PREFIX_SHORT_NAME = "pre";

    protected static final String DO_FILTER_SEGMENTS_LONG_NAME = "filterSegments";
    protected static final String DO_FILTER_SEGMENTS_SHORT_NAME = "filterSeg";

    protected static final String LENGTH_PERCENTILE_LONG_NAME = "lengthPercentile";
    protected static final String LENGTH_PERCENTILE_SHORT_NAME = "lenPct";

    protected static final String LOG2_COPY_RATIO_CREDIBLE_INTERVAL_PERCENTILE_LONG_NAME = "log2CopyRatioCredibleIntervalPercentile";
    protected static final String LOG2_COPY_RATIO_CREDIBLE_INTERVAL_PERCENTILE_SHORT_NAME = "log2CRCredIntPct";

    protected static final String MINOR_ALLELE_FRACTION_CREDIBLE_INTERVAL_PERCENTILE_LONG_NAME = "minorAlleleFractionCredibleIntervalPercentile";
    protected static final String MINOR_ALLELE_FRACTION_CREDIBLE_INTERVAL_PERCENTILE_SHORT_NAME = "MAFCredIntPct";

    protected static final String MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME = "maxAllelicCopyNumberClonal";
    protected static final String MAX_ALLELIC_COPY_NUMBER_CLONAL_SHORT_NAME = "maxACNClonal";

    protected static final String MAX_ALLELIC_COPY_NUMBER_LONG_NAME = "maxAllelicCopyNumber";
    protected static final String MAX_ALLELIC_COPY_NUMBER_SHORT_NAME = "maxACN";

    protected static final String MAX_NUM_POPULATIONS_LONG_NAME = "maxNumPopulations";
    protected static final String MAX_NUM_POPULATIONS_SHORT_NAME = "maxNumPop";

    protected static final String NUM_CELLS_LONG_NAME = "numCells";
    protected static final String NUM_CELLS_SHORT_NAME = "numCells";

    protected static final String SWAP_ITERATION_DIVISOR_CLONAL_LONG_NAME = "swapIterationDivisorClonal";
    protected static final String SWAP_ITERATION_DIVISOR_CLONAL_SHORT_NAME = "swapIterDivisorClonal";

    protected static final String SWAP_ITERATION_DIVISOR_LONG_NAME = "swapIterationDivisor";
    protected static final String SWAP_ITERATION_DIVISOR_SHORT_NAME = "swapIterDivisor";

    protected static final String NUM_SAMPLES_CLONAL_LONG_NAME = "numSamplesClonal";
    protected static final String NUM_SAMPLES_CLONAL_SHORT_NAME = "numSampClonal";

    protected static final String NUM_BURN_IN_CLONAL_LONG_NAME = "numBurnInClonal";
    protected static final String NUM_BURN_IN_CLONAL_SHORT_NAME = "numBurnClonal";

    protected static final String NUM_SAMPLES_LONG_NAME = "numSamples";
    protected static final String NUM_SAMPLES_SHORT_NAME = "numSamp";

    protected static final String NUM_BURN_IN_LONG_NAME = "numBurnIn";
    protected static final String NUM_BURN_IN_SHORT_NAME = "numBurn";
    
    protected static final String CONCENTRATION_PRIOR_ALPHA_LONG_NAME = "concentrationPriorAlpha";
    protected static final String CONCENTRATION_PRIOR_ALPHA_SHORT_NAME = "concAlpha";

    protected static final String CONCENTRATION_PRIOR_BETA_LONG_NAME = "concentrationPriorBeta";
    protected static final String CONCENTRATION_PRIOR_BETA_SHORT_NAME = "concBeta";

    protected static final String VARIANT_SEGMENT_FRACTION_PRIOR_ALPHA_LONG_NAME = "variantSegmentFractionPriorAlpha";
    protected static final String VARIANT_SEGMENT_FRACTION_PRIOR_ALPHA_SHORT_NAME = "varSegFracAlpha";

    protected static final String VARIANT_SEGMENT_FRACTION_PRIOR_BETA_LONG_NAME = "variantSegmentFractionPriorBeta";
    protected static final String VARIANT_SEGMENT_FRACTION_PRIOR_BETA_SHORT_NAME = "varSegFracBeta";

    protected static final String PLOIDY_STATE_PRIOR_COMPLETE_DELETION_PENALTY_LONG_NAME = "completeDeletionPenalty";
    protected static final String PLOIDY_STATE_PRIOR_COMPLETE_DELETION_PENALTY_SHORT_NAME = "compDelPen";

    protected static final String PLOIDY_STATE_PRIOR_CHANGE_PENALTY_LONG_NAME = "changePenalty";
    protected static final String PLOIDY_STATE_PRIOR_CHANGE_PENALTY_SHORT_NAME = "changePen";

    @Argument(
            doc = "Input file for AllelicCNV result.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File allelicCNVFile;

    @Argument(
            doc = "Prefix for output files.",
            fullName = OUTPUT_PREFIX_LONG_NAME,
            shortName = OUTPUT_PREFIX_SHORT_NAME,
            optional = false
    )
    protected String outputPrefix;

    @Argument(
            doc = "Filter AllelicCNV segments based on length and credible interval.",
            fullName = DO_FILTER_SEGMENTS_LONG_NAME,
            shortName = DO_FILTER_SEGMENTS_SHORT_NAME,
            optional = true
    )
    protected boolean doFilterSegments = false;

    @Argument(
            doc = "Percentile of segment length below which to filter out segments.",
            fullName = LENGTH_PERCENTILE_LONG_NAME,
            shortName = LENGTH_PERCENTILE_SHORT_NAME,
            optional = true
    )
    protected double lengthPercentile = 0.;

    @Argument(
            doc = "Percentile of log2 copy-ratio credible-interval size above which to filter out segments.",
            fullName = LOG2_COPY_RATIO_CREDIBLE_INTERVAL_PERCENTILE_LONG_NAME,
            shortName = LOG2_COPY_RATIO_CREDIBLE_INTERVAL_PERCENTILE_SHORT_NAME,
            optional = true
    )
    protected double log2CrCredibleIntervalPercentile = 100.;

    @Argument(
            doc = "Percentile of minor-allele-fraction credible-interval size above which to filter out segments.",
            fullName = MINOR_ALLELE_FRACTION_CREDIBLE_INTERVAL_PERCENTILE_LONG_NAME,
            shortName = MINOR_ALLELE_FRACTION_CREDIBLE_INTERVAL_PERCENTILE_SHORT_NAME,
            optional = true
    )
    protected double mafCredibleIntervalPercentile = 100.;

    @Argument(
            doc = "Maximum allelic copy number for clonal model.",
            fullName = MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME,
            shortName = MAX_ALLELIC_COPY_NUMBER_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int maxAllelicCopyNumberClonal = 5;

    @Argument(
            doc = "Maximum allelic copy number for full model.",
            fullName = MAX_ALLELIC_COPY_NUMBER_LONG_NAME,
            shortName = MAX_ALLELIC_COPY_NUMBER_SHORT_NAME,
            optional = true
    )
    protected int maxAllelicCopyNumber = 5;

    @Argument(
            doc = "Number of auxilliary cell indicators.",
            fullName = NUM_CELLS_LONG_NAME,
            shortName = NUM_CELLS_SHORT_NAME,
            optional = true
    )
    protected int numCells = 50;

    @Argument(
            doc = "Swap step for clonal model will be proposed if randomly-generated integer is divisible by this number.",
            fullName = SWAP_ITERATION_DIVISOR_CLONAL_LONG_NAME,
            shortName = SWAP_ITERATION_DIVISOR_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int swapIterationDivisorClonal = 2;

    @Argument(
            doc = "Swap step for full model will be proposed if randomly-generated integer is divisible by this number.",
            fullName = SWAP_ITERATION_DIVISOR_LONG_NAME,
            shortName = SWAP_ITERATION_DIVISOR_SHORT_NAME,
            optional = true
    )
    protected int swapIterationDivisor = 2;

    @Argument(
            doc = "Maximum number of populations for full model.",
            fullName = MAX_NUM_POPULATIONS_LONG_NAME,
            shortName = MAX_NUM_POPULATIONS_SHORT_NAME,
            optional = true
    )
    protected int maxNumPopulations = 4;

    @Argument(
            doc = "Total number of MCMC samples for clonal model.",
            fullName = NUM_SAMPLES_CLONAL_LONG_NAME,
            shortName = NUM_SAMPLES_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int numSamplesClonal = 50;

    @Argument(
            doc = "Number of burn-in samples to discard for clonal model.",
            fullName = NUM_BURN_IN_CLONAL_LONG_NAME,
            shortName = NUM_BURN_IN_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int numBurnInClonal = 25;

    @Argument(
            doc = "Total number of MCMC samples for full model.",
            fullName = NUM_SAMPLES_LONG_NAME,
            shortName = NUM_SAMPLES_SHORT_NAME,
            optional = true
    )
    protected int numSamples = 500;

    @Argument(
            doc = "Number of burn-in samples to discard for full model.",
            fullName = NUM_BURN_IN_LONG_NAME,
            shortName = NUM_BURN_IN_SHORT_NAME,
            optional = true
    )
    protected int numBurnIn = 400;

    @Argument(
            doc = "Alpha hyperparameter for Gamma-distribution prior on concentration parameter.",
            fullName = CONCENTRATION_PRIOR_ALPHA_LONG_NAME,
            shortName = CONCENTRATION_PRIOR_ALPHA_SHORT_NAME,
            optional = true
    )
    protected double concentrationPriorAlpha = 1.;

    @Argument(
            doc = "Beta hyperparameter for Gamma-distribution prior on concentration parameter.",
            fullName = CONCENTRATION_PRIOR_BETA_LONG_NAME,
            shortName = CONCENTRATION_PRIOR_BETA_SHORT_NAME,
            optional = true
    )
    protected double concentrationPriorBeta = 1E8;

    @Argument(
            doc = "Alpha hyperparameter for Beta-distribution prior on variant-segment-fraction parameter.",
            fullName = VARIANT_SEGMENT_FRACTION_PRIOR_ALPHA_LONG_NAME,
            shortName = VARIANT_SEGMENT_FRACTION_PRIOR_ALPHA_SHORT_NAME,
            optional = true
    )
    protected double variantSegmentFractionPriorAlpha = 5.;

    @Argument(
            doc = "Beta hyperparameter for Beta-distribution prior on variant-segment-fraction parameter.",
            fullName = VARIANT_SEGMENT_FRACTION_PRIOR_BETA_LONG_NAME,
            shortName = VARIANT_SEGMENT_FRACTION_PRIOR_BETA_SHORT_NAME,
            optional = true
    )
    protected double variantSegmentFractionPriorBeta = 5.;

    @Argument(
            doc = "Penalty for complete allele deletion in ploidy-state prior.",
            fullName = PLOIDY_STATE_PRIOR_COMPLETE_DELETION_PENALTY_LONG_NAME,
            shortName = PLOIDY_STATE_PRIOR_COMPLETE_DELETION_PENALTY_SHORT_NAME,
            optional = true
    )
    protected double ploidyStatePriorCompleteDeletionPenalty = 1E-3;

    @Argument(
            doc = "Penalty for copy change in ploidy-state prior.",
            fullName = PLOIDY_STATE_PRIOR_CHANGE_PENALTY_LONG_NAME,
            shortName = PLOIDY_STATE_PRIOR_CHANGE_PENALTY_SHORT_NAME,
            optional = true
    )
    protected double ploidyStatePriorChangePenalty = 1E-3;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        validateArguments();

        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

        final File filteredSegmentsFile = new File(outputPrefix + FILTERED_SEGMENTS_FILE_SUFFIX);
        final List<ACNVModeledSegment> allSegments = SegmentUtils.readACNVModeledSegmentFile(allelicCNVFile);
        final List<ACNVModeledSegment> segments = doFilterSegments ?
                filterAndOutputSegments(allSegments, filteredSegmentsFile, logger, lengthPercentile, log2CrCredibleIntervalPercentile, mafCredibleIntervalPercentile) :
                allSegments;

        final TumorHeterogeneityData data = new TumorHeterogeneityData(segments);
        final File mafCrFile = new File(outputPrefix + MAF_CR_FILE_SUFFIX);
        outputMafCrFile(mafCrFile, data);

        final PloidyStatePrior variantPloidyStatePriorClonal = calculatePloidyStatePrior(ploidyStatePriorCompleteDeletionPenalty, ploidyStatePriorChangePenalty, maxAllelicCopyNumberClonal);
        final TumorHeterogeneityPriorCollection priorsClonal = new TumorHeterogeneityPriorCollection(NORMAL_PLOIDY_STATE, variantPloidyStatePriorClonal,
                concentrationPriorAlpha, concentrationPriorBeta, variantSegmentFractionPriorAlpha, variantSegmentFractionPriorBeta);

        final PloidyStatePrior variantPloidyStatePrior = calculatePloidyStatePrior(ploidyStatePriorCompleteDeletionPenalty, ploidyStatePriorChangePenalty, maxAllelicCopyNumber);
        final TumorHeterogeneityPriorCollection priors = new TumorHeterogeneityPriorCollection(NORMAL_PLOIDY_STATE, variantPloidyStatePrior,
                concentrationPriorAlpha, concentrationPriorBeta, variantSegmentFractionPriorAlpha, variantSegmentFractionPriorBeta);

        final File resultClonalFile = new File(outputPrefix + CLONAL_RESULT_FILE_SUFFIX);
        final TumorHeterogeneityModeller clonalModeller = new TumorHeterogeneityModeller(data, priorsClonal, NUM_POPULATIONS_CLONAL, numCells, swapIterationDivisorClonal, rng);
        clonalModeller.fitMCMC(numSamplesClonal, numBurnInClonal);
        clonalModeller.output(resultClonalFile, CREDIBLE_INTERVAL_ALPHA, ctx);

        final TumorHeterogeneityState initialState = initializeStateFromClonalResult(data, priors, clonalModeller, maxNumPopulations, numCells);
        final File resultFile = new File(outputPrefix + RESULT_FILE_SUFFIX);
        final TumorHeterogeneityModeller modeller = new TumorHeterogeneityModeller(data, initialState, swapIterationDivisor, rng);
        modeller.fitMCMC(numSamples, numBurnIn);
        modeller.output(resultFile, CREDIBLE_INTERVAL_ALPHA, ctx);

        logger.info("SUCCESS: Tumor heterogeneity run complete and result output to " + resultFile + ".");
    }

    private static void output(final FileWriter writer, final Logger logger, final String output) {
        try {
            writer.write(output);
//            logger.debug(output);
        } catch(final IOException e) {
            throw new GATKException("Cannot output.");
        }
    }

    private static void outputMafCrFile(final File mafCrFile,
                                        final TumorHeterogeneityData data) {
        try (final FileWriter mafCrWriter = new FileWriter(mafCrFile)) {
            mafCrFile.createNewFile();
            for (double f = 0.; f <= 0.5; f += 0.01) {
                for (double c = 1E-10; c <= 5; c += 0.05) {
                    final double maf = f;
                    final double cr = c;
                    final double density = IntStream.range(0, data.numSegments()).mapToDouble(i -> Math.exp(data.logDensity(i, cr, maf))).sum();
                    mafCrWriter.write(maf + "\t" + cr + "\t" + density + System.getProperty("line.separator"));
                }
            }
        } catch (final IOException e) {
            throw new GATKException("Error writing MAF-CR file.");
        }
    }

    private static PloidyStatePrior calculatePloidyStatePrior(final double ploidyStatePriorCompleteDeletionPenalty,
                                                              final double ploidyStatePriorChangePenalty,
                                                              final int maxAllelicCopyNumber) {
        final Function<PloidyState, Double> ploidyLogPDF = ps -> Math.log(ploidyStatePriorCompleteDeletionPenalty) * ((ps.m() == 0 ? 1 : 0) + (ps.n() == 0 ? 1 : 0))
                + Math.log(ploidyStatePriorChangePenalty) * (Math.abs(NORMAL_PLOIDY_STATE.m() - ps.m()) + Math.abs(NORMAL_PLOIDY_STATE.n() - ps.n()));
        final Map<PloidyState, Double> unnormalizedLogProbabilityMassFunctionMap = new LinkedHashMap<>();
        for (int n = 0; n <= maxAllelicCopyNumber; n++) {
            for (int m = 0; m <= n; m++) {
                if (!(m == 1 && n == 1)) {
                    unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(m, n), ploidyLogPDF.apply(new PloidyState(m, n)));
                }
            }
        }
        return new PloidyStatePrior(unnormalizedLogProbabilityMassFunctionMap);
    }

    private static List<ACNVModeledSegment> filterAndOutputSegments(final List<ACNVModeledSegment> allSegments, final File outputFile, final Logger logger,
                                                                    final double lengthPercentile, final double log2CrCredibleIntervalPercentile, final double mafCredibleIntervalPercentile) {
        try (final FileWriter writer = new FileWriter(outputFile)) {
            outputFile.createNewFile();
            final Percentile percentile = new Percentile();

            final double[] lengths = allSegments.stream().mapToDouble(s -> (double) s.getInterval().size()).toArray();
            final int lengthThreshold = (int) percentile.evaluate(lengths, lengthPercentile);

            final double[] log2CrCredibleIntervalSizes = allSegments.stream()
                    .mapToDouble(s -> s.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_90) - s.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_10))
                    .toArray();
            final double log2crCredibleIntervalThreshold = percentile.evaluate(log2CrCredibleIntervalSizes, log2CrCredibleIntervalPercentile);

            final double[] mafCredibleIntervalSizes = allSegments.stream()
                    .filter(s -> !Double.isNaN(s.getMinorAlleleFractionPosteriorSummary().getCenter()))
                    .mapToDouble(s -> s.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_90) - s.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_10))
                    .toArray();
            final double mafCredibleIntervalThreshold = percentile.evaluate(mafCredibleIntervalSizes, mafCredibleIntervalPercentile);

            final List<ACNVModeledSegment> segments = allSegments.stream()
                    .filter(s -> s.getInterval().size() > lengthThreshold)
                    .filter(s -> s.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_90) - s.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_10) < log2crCredibleIntervalThreshold)
                    .filter(s -> Double.isNaN(s.getMinorAlleleFractionPosteriorSummary().getCenter()) || s.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_90) - s.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_10) < mafCredibleIntervalThreshold)
                    .collect(Collectors.toList());

            output(writer, logger, "#num segments all: " + allSegments.size());
            output(writer, logger, System.getProperty("line.separator"));
            output(writer, logger, "#num segments after filtering: " + segments.size());
            output(writer, logger, System.getProperty("line.separator"));
            output(writer, logger, "#length threshold: " + lengthThreshold);
            output(writer, logger, System.getProperty("line.separator"));
            output(writer, logger, "#log2CR credible-interval threshold: " + log2crCredibleIntervalThreshold);
            output(writer, logger, System.getProperty("line.separator"));
            output(writer, logger, "#MAF credible-interval threshold: " + mafCredibleIntervalThreshold);
            output(writer, logger, System.getProperty("line.separator"));

            output(writer, logger, SegmentTableColumn.CONTIG + "\t" + SegmentTableColumn.START + "\t" + SegmentTableColumn.END + "\t" +
                    SegmentTableColumn.SEGMENT_MEAN_POSTERIOR_MODE + "\t" + SegmentTableColumn.SEGMENT_MEAN_POSTERIOR_LOWER + "\t" + SegmentTableColumn.SEGMENT_MEAN_POSTERIOR_UPPER + "\t" +
                    SegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_MODE + "\t" + SegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_LOWER + "\t" + SegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_UPPER);
            allSegments.stream().filter(s -> !segments.contains(s))
                    .forEach(s -> {
                        output(writer, logger, s.getContig() + "\t" + s.getStart() + "\t" + s.getEnd() + "\t" +
                                s.getSegmentMeanPosteriorSummary().getCenter() + "\t" + s.getSegmentMeanPosteriorSummary().getLower() + "\t" + s.getSegmentMeanPosteriorSummary().getUpper() + "\t" +
                                s.getMinorAlleleFractionPosteriorSummary().getCenter() + "\t" + s.getMinorAlleleFractionPosteriorSummary().getLower() + "\t" + s.getMinorAlleleFractionPosteriorSummary().getUpper());
                        output(writer, logger, System.getProperty("line.separator"));
                    });

            return segments;
        } catch (final IOException e) {
            throw new GATKException("Error writing filtered segments.");
        }
    }

    private static TumorHeterogeneityState initializeStateFromClonalResult(final TumorHeterogeneityData data,
                                                                           final TumorHeterogeneityPriorCollection priors,
                                                                           final TumorHeterogeneityModeller clonalModeller,
                                                                           final int maxNumPopulations,
                                                                           final int numCells) {
        final double clonalConcentration = Iterables.getLast(clonalModeller.getConcentrationSamples());
        final double clonalNormalFraction = Iterables.getLast(Iterables.getLast(clonalModeller.getPopulationFractionsSamples()));
        final TumorHeterogeneityState.PopulationIndicators initialPopulationIndicators = new TumorHeterogeneityState.PopulationIndicators(Iterables.getLast(clonalModeller.getPopulationIndicatorsSamples()));
        IntStream.range(0, numCells).filter(i -> initialPopulationIndicators.get(i) == 1).forEach(i -> initialPopulationIndicators.set(i, maxNumPopulations - 1));
        final TumorHeterogeneityState.VariantProfileCollection clonalVariantProfileCollection = Iterables.getLast(clonalModeller.getVariantProfileCollectionSamples());
        final List<Double> initialFractions = new ArrayList<>();
        initialFractions.add(1. - clonalNormalFraction);
        final List<TumorHeterogeneityState.VariantProfile> initialVariantProfiles = new ArrayList<>();
        initialVariantProfiles.addAll(clonalVariantProfileCollection);
        for (int i = 0; i < maxNumPopulations - NUM_POPULATIONS_CLONAL; i++) {
            initialFractions.add(1, 0.);
            initialVariantProfiles.add(1, TumorHeterogeneityModeller.initializeNormalProfile(data.numSegments()));
        }
        initialFractions.add(clonalNormalFraction);
        final TumorHeterogeneityState.PopulationFractions initialPopulationFractions =
                new TumorHeterogeneityState.PopulationFractions(initialFractions);
        final TumorHeterogeneityState.VariantProfileCollection initialVariantProfileCollection =
                new TumorHeterogeneityState.VariantProfileCollection(initialVariantProfiles);
        return new TumorHeterogeneityState(clonalConcentration, initialPopulationFractions, initialPopulationIndicators, initialVariantProfileCollection, priors);
    }

    //validate CLI arguments
    private void validateArguments() {
        Utils.validateArg(0. <= lengthPercentile && lengthPercentile <= 100., LENGTH_PERCENTILE_LONG_NAME + " must be in [0, 100].");
        Utils.validateArg(0. <= log2CrCredibleIntervalPercentile && log2CrCredibleIntervalPercentile <= 100., LOG2_COPY_RATIO_CREDIBLE_INTERVAL_PERCENTILE_LONG_NAME + " must be in [0, 100].");
        Utils.validateArg(0. <= mafCredibleIntervalPercentile && mafCredibleIntervalPercentile <= 100., MINOR_ALLELE_FRACTION_CREDIBLE_INTERVAL_PERCENTILE_LONG_NAME + " must be in [0, 100].");
        Utils.validateArg(maxAllelicCopyNumberClonal > 0, MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME + " must be positive.");
        Utils.validateArg(maxAllelicCopyNumber > 0, MAX_ALLELIC_COPY_NUMBER_LONG_NAME + " must be positive.");
        Utils.validateArg(maxNumPopulations >= 2, MAX_NUM_POPULATIONS_LONG_NAME + " must be greater than or equal to 2.");
        Utils.validateArg(numCells > 0, NUM_CELLS_LONG_NAME + " must be positive.");
        Utils.validateArg(swapIterationDivisorClonal > 0, SWAP_ITERATION_DIVISOR_CLONAL_LONG_NAME + " must be positive.");
        Utils.validateArg(swapIterationDivisor > 0, SWAP_ITERATION_DIVISOR_LONG_NAME + " must be positive.");
        Utils.validateArg(numSamplesClonal > 0, NUM_SAMPLES_CLONAL_LONG_NAME + " must be positive.");
        Utils.validateArg(numBurnInClonal > 0 && numBurnInClonal <= numSamplesClonal, NUM_BURN_IN_CLONAL_LONG_NAME + " must be positive and less than or equal to " + NUM_SAMPLES_CLONAL_LONG_NAME);
        Utils.validateArg(numSamples > 0, NUM_SAMPLES_LONG_NAME + " must be positive.");
        Utils.validateArg(numBurnIn > 0 && numBurnIn <= numSamples, NUM_BURN_IN_LONG_NAME + " must be positive and less than or equal to " + NUM_SAMPLES_LONG_NAME);
        Utils.validateArg(concentrationPriorAlpha > 0, CONCENTRATION_PRIOR_ALPHA_LONG_NAME + " must be positive.");
        Utils.validateArg(concentrationPriorBeta > 0, CONCENTRATION_PRIOR_BETA_LONG_NAME + " must be positive.");
        Utils.validateArg(TumorHeterogeneityModeller.CONCENTRATION_MIN < concentrationPriorAlpha / concentrationPriorBeta &&
                concentrationPriorAlpha / concentrationPriorBeta < TumorHeterogeneityModeller.CONCENTRATION_MAX,
                CONCENTRATION_PRIOR_ALPHA_LONG_NAME + " / " + CONCENTRATION_PRIOR_BETA_LONG_NAME + " must be in (" +
                        TumorHeterogeneityModeller.CONCENTRATION_MIN + ", " + TumorHeterogeneityModeller.CONCENTRATION_MAX + ").");
        Utils.validateArg(variantSegmentFractionPriorAlpha > 0, VARIANT_SEGMENT_FRACTION_PRIOR_ALPHA_LONG_NAME + " must be positive.");
        Utils.validateArg(variantSegmentFractionPriorBeta > 0, VARIANT_SEGMENT_FRACTION_PRIOR_BETA_LONG_NAME + " must be positive.");
        Utils.validateArg(ploidyStatePriorCompleteDeletionPenalty > 0, PLOIDY_STATE_PRIOR_COMPLETE_DELETION_PENALTY_LONG_NAME + " must be positive.");
        Utils.validateArg(ploidyStatePriorChangePenalty > 0, PLOIDY_STATE_PRIOR_CHANGE_PENALTY_LONG_NAME + " must be positive.");
    }
}