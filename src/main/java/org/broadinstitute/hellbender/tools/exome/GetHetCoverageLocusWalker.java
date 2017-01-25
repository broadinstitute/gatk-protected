package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountTableColumn;
import org.broadinstitute.hellbender.tools.exome.pulldown.HetPulldownCalculator;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

@CommandLineProgramProperties(
        summary = "EXPERIMENTAL ... Output ref/alt counts at heterozygous SNPs in normal sample (and at same sites in tumor sample, if specified).",
        oneLineSummary = "Output ref/alt counts at heterozygous SNPs in normal sample (and at same sites in tumor sample, if specified)",
        programGroup = CopyNumberProgramGroup.class
)
public final class GetHetCoverageLocusWalker extends LocusWalker {

    protected static final String PVALUE_THRESHOLD_FULL_NAME = "pvalueThreshold";
    protected static final String PVALUE_THRESHOLD_SHORT_NAME = "p";

    protected static final String MINIMUM_MAPPING_QUALITY_SHORT_NAME = "minMQ";
    protected static final String MINIMUM_MAPPING_QUALITY_FULL_NAME = "minimumMappingQuality";

    protected static final String MINIMUM_BASE_QUALITY_SHORT_NAME = "minBQ";
    protected static final String MINIMUM_BASE_QUALITY_FULL_NAME = "minimumBaseQuality";

    protected static final String MINIMUM_READ_COUNT_SHORT_NAME = "minRC";
    protected static final String MINIMUM_READ_COUNT_FULL_NAME = "minimumReadCount";

    @Argument(
            doc = "BAM file for normal sample.",
            fullName = ExomeStandardArgumentDefinitions.NORMAL_BAM_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.NORMAL_BAM_FILE_SHORT_NAME,
            optional = false
    )
    protected File normalBAMFile;

    @Argument(
            doc = "BAM file for tumor sample.",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_SHORT_NAME,
            optional = true
    )
    protected File tumorBAMFile;

//    @Argument(
//            doc = "Interval-list file of common SNPs.",
//            fullName = ExomeStandardArgumentDefinitions.SNP_FILE_LONG_NAME,
//            shortName = ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME,
//            optional = false
//    )
//    protected File snpFile;

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
            doc = "Minimum mapping quality; reads with lower quality will be filtered out of pileup.",
            shortName = MINIMUM_MAPPING_QUALITY_SHORT_NAME,
            fullName  = MINIMUM_MAPPING_QUALITY_FULL_NAME,
            optional = true
    )
    protected int minimumMappingQuality = 30;

    @Argument(
            doc = "Minimum base quality; base calls with lower quality will be filtered out of pileup.",
            shortName = MINIMUM_BASE_QUALITY_SHORT_NAME,
            fullName = MINIMUM_BASE_QUALITY_FULL_NAME,
            optional = true
    )
    protected int minimumBaseQuality = 20;

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

    private AllelicCountCollection hetPulldown = new AllelicCountCollection();

    @Override
    public void onTraversalStart() {
        logger.info("Starting...");
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {

        // TODO: Error checking
        // TODO: Do I need to filter the reads?
        final List<GATKRead> allReads = alignmentContext.getBasePileup().getReads();

        if (allReads.size() < minimumRawReads) {
            return;
        }
        final byte refAsByte = referenceContext.getBase();

        final List<Integer> allOffsets = alignmentContext.getBasePileup().getOffsets();

        final Nucleotide.Counter baseCounts = new Nucleotide.Counter();
        IntStream.range(0, allReads.size()).forEach(i-> baseCounts.add(allReads.get(i).getBase(allOffsets.get(i))));
        final int totalBaseCount = Arrays.stream(HetPulldownCalculator.BASES).mapToInt(b -> (int) baseCounts.get(b)).sum();

        // TODO: Make sure we are working with the normal
        if (!HetPulldownCalculator.isPileupHetCompatible(baseCounts, totalBaseCount, pvalThreshold)) {
            return;
        } else {

            final int refReadCount = (int) baseCounts.get(Nucleotide.valueOf(refAsByte));
            final int altReadCount = totalBaseCount - refReadCount;
            hetPulldown.add(new AllelicCount(referenceContext.getInterval(), refReadCount, altReadCount));

        }
    }

    @Override
    public Object onTraversalSuccess() {
        hetPulldown.write(normalHetOutputFile, AllelicCountTableColumn.AllelicCountTableVerbosity.BASIC);
        return "SUCCESS";
    }
}

