package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountTableColumn;
import org.broadinstitute.hellbender.tools.exome.pulldown.HetPulldownCalculator;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadConstants;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

@CommandLineProgramProperties(
        summary = "EXPERIMENTAL ... Output ref/alt counts at heterozygous SNPs in normal sample (and at same sites in tumor sample, if specified).\n" +
                "\t\tNOTE:  This method assumes that the first input bam file is the normal (control) and the second is the tumor (case)." +
                "\n\t\tNOTE:  This method requires the SM tag from the normal bam file in order to work.  ",
        oneLineSummary = "Output ref/alt counts at heterozygous SNPs in normal sample (and at same sites in tumor sample, if specified)",
        programGroup = CopyNumberProgramGroup.class
)
public final class GetHetCoverageLocusWalker extends LocusWalker {

    protected static final String PVALUE_THRESHOLD_FULL_NAME = "pvalueThreshold";
    protected static final String PVALUE_THRESHOLD_SHORT_NAME = "p";

    protected static final String MINIMUM_READ_COUNT_SHORT_NAME = "minRC";
    protected static final String MINIMUM_READ_COUNT_FULL_NAME = "minimumReadCount";

    @Argument(
            doc = "Output file for normal-sample ref/alt read counts (at heterozygous SNPs).",
            fullName = ExomeStandardArgumentDefinitions.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.NORMAL_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File normalHetOutputFile;

    @Argument(
            doc = "Output file for tumor-sample ref/alt read counts (at heterozygous SNPs in normal sample).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = true
    )
    protected File tumorHetOutputFile;

    @Argument(
            doc = "p-value threshold for binomial test for heterozygous SNPs in normal sample (must be in [0, 1]).",
            fullName = PVALUE_THRESHOLD_FULL_NAME,
            shortName = PVALUE_THRESHOLD_SHORT_NAME,
            optional = false
    )
    protected double pvalThreshold = 0.05;

    @Argument(
            doc = "Validation stringency for all BAM files read by this program.  Setting stringency to SILENT " +
                    "can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) " +
                    "do not otherwise need to be decoded.",
            common=true)
    protected ValidationStringency VALIDATION_STRINGENCY = ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY;

    @Argument(
            doc = "Minimum raw number of reads that must be present at a site for it to be considered for het compatibility.",
            shortName = MINIMUM_READ_COUNT_SHORT_NAME,
            fullName = MINIMUM_READ_COUNT_FULL_NAME,
            optional = true
    )
    protected int minimumRawReads = 15;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    private AllelicCountCollection controlHetPulldown = new AllelicCountCollection();

    private Map<String, AllelicCountCollection> caseSampleNameToHetPulldowns = new HashMap<>();

    private String controlSampleName = null;


    @Override
    public void onTraversalStart() {
        // TODO:  HACK:  get the normal (first -I) file ... this will not be needed once conventions are changed
        final List<Path> readArguments = this.readArguments.getReadPaths();

        if (readArguments.size() > 2) {
            throw new UserException.BadInput("Please specify only two bam files.  The first is the normal and the second is the tumor.  Changes to this convention are coming soon.");
        }

        // Get the normal sample name
        final ReadsDataSource normalReads = new ReadsDataSource(readArguments.get(0));
        controlSampleName = ReadUtils.getSamplesFromHeader(normalReads.getHeader()).iterator().next();

        // Get the tumor sample name
        final List<Path> casePaths = readArguments.subList(1, readArguments.size());
        for (final Path casePath : casePaths) {
            final ReadsDataSource caseReads = new ReadsDataSource(casePath);
            final String sampleName = ReadUtils.getSamplesFromHeader(caseReads.getHeader()).iterator().next();

            caseSampleNameToHetPulldowns.put(sampleName, new AllelicCountCollection());
        }

        logger.info("Starting... N: " + controlSampleName + "   T: " + String.join(", ", caseSampleNameToHetPulldowns.keySet()));

    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {

        final Map<String, AlignmentContext> sampleContexts = alignmentContext.splitContextBySampleName(getHeaderForReads());

        // TODO: Check that no two samples have the same name.
        // TODO: edge case where my normal is an empty bam
        // TODO: Error checking
        // TODO: Do I need to filter the reads?

        // Get the reads from the normal
        final AlignmentContext alignmentContextForNormal = sampleContexts.get(controlSampleName);
        if (alignmentContextForNormal == null) {
            return;
        }
        final List<GATKRead> allNormalReads = alignmentContextForNormal.getBasePileup().getReads();

        if (allNormalReads.size() < minimumRawReads) {
            return;
        }
        final byte refAsByte = referenceContext.getBase();

        final List<Integer> allNormalOffsets = alignmentContextForNormal.getBasePileup().getOffsets();

        final Nucleotide.Counter normalBaseCounts = new Nucleotide.Counter();
        IntStream.range(0, allNormalReads.size()).forEach(i -> normalBaseCounts.add(allNormalReads.get(i).getBase(allNormalOffsets.get(i))));
        final int totalBaseCountInNormal = Arrays.stream(HetPulldownCalculator.BASES).mapToInt(b -> (int) normalBaseCounts.get(b)).sum();

        if (HetPulldownCalculator.isPileupHetCompatible(normalBaseCounts, totalBaseCountInNormal, pvalThreshold)) {

            final int refReadCountInNormal = (int) normalBaseCounts.get(Nucleotide.valueOf(refAsByte));
            final int altReadCountInNormal = totalBaseCountInNormal - refReadCountInNormal;
            controlHetPulldown.add(new AllelicCount(referenceContext.getInterval(), refReadCountInNormal, altReadCountInNormal));

            // Process case samples
            for (final String caseSampleName : caseSampleNameToHetPulldowns.keySet()) {
                final AlignmentContext alignmentContextForCase = sampleContexts.get(caseSampleName);
                if (alignmentContextForCase == null) {
                    continue;
                }
                final List<GATKRead> allCaseReads = alignmentContextForCase.getBasePileup().getReads();
                if (allCaseReads.size() < minimumRawReads) {
                    continue;
                }
                final List<Integer> allCaseOffsets = alignmentContextForCase.getBasePileup().getOffsets();
                final Nucleotide.Counter caseBaseCounts = new Nucleotide.Counter();
                IntStream.range(0, allCaseReads.size()).forEach(i -> caseBaseCounts.add(allCaseReads.get(i).getBase(allCaseOffsets.get(i))));
                final int totalBaseCountInCase = Arrays.stream(HetPulldownCalculator.BASES).mapToInt(b -> (int) caseBaseCounts.get(b)).sum();
                final int refReadCountInCase = (int) caseBaseCounts.get(Nucleotide.valueOf(refAsByte));
                final int altReadCountInCase = totalBaseCountInCase - refReadCountInCase;
                caseSampleNameToHetPulldowns.get(caseSampleName).add(new AllelicCount(referenceContext.getInterval(), refReadCountInCase, altReadCountInCase));
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Writing normal het file: " + normalHetOutputFile.getAbsolutePath());
        controlHetPulldown.write(normalHetOutputFile, AllelicCountTableColumn.AllelicCountTableVerbosity.BASIC);
        // TODO: Update to allow more than one case sample.  This next line assumes one case sample, since it writes to the same output file no matter what.
        for (final Map.Entry<String, AllelicCountCollection> caseEntry : caseSampleNameToHetPulldowns.entrySet()) {
            logger.info("Writing case het file: " + tumorHetOutputFile.getAbsolutePath());
            caseEntry.getValue().write(tumorHetOutputFile, AllelicCountTableColumn.AllelicCountTableVerbosity.BASIC);
        }
        return "SUCCESS";
    }
}

