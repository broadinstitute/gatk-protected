package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IntervalList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import java.io.File;

/**
 * Calls heterozygous sites and output the ref and alt base counts at the called sites.
 * A list of SNP sites must be provided.
 *
 * The tool has three modes. The mode is automatically selected based on the provided arguments:
 *
 *   NORMAL_ONLY: this mode is chosen if just the normal BAM file is provided. The BALANCED prior
 *                is used for estimating the likelihood of het sites.
 *
 *   TUMOR_ONLY: this mode is chosen if just the tumor BAM file is provided. The HETEROGENEOUS prior
 *               is used for estimating the likelihood of het sites.
 *
 *   NORMAL_TUMOR: this mode is chosen if both normal and tumor BAM files are provided. The BALANCED prior is used
 *                 to estimate the likelihood of het sites on the normal reads. The base counts on the same het sites
 *                 are collected and reported from the tumor reads.
 *
 * @author Mehrtash Babadi &lt;mehrtash@brgoadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Calls heterozygous sites and output the ref and alt base counts at the called sites.",
        oneLineSummary = "Calls heterozygous sites and output the ref and alt base counts at the called sites.",
        programGroup = CopyNumberProgramGroup.class
)
public final class GetBayesianHetCoverage extends CommandLineProgram {

    private final Logger logger = LogManager.getLogger(GetBayesianHetCoverage.class);

    protected static final String READ_DEPTH_THRESHOLD_FULL_NAME = "readDepthThreshold";
    protected static final String READ_DEPTH_THRESHOLD_SHORT_NAME = "readThresh";

    protected static final String MINIMUM_MAPPING_QUALITY_FULL_NAME = "minimumMappingQuality";
    protected static final String MINIMUM_MAPPING_QUALITY_SHORT_NAME = "minMQ";

    protected static final String MINIMUM_BASE_QUALITY_FULL_NAME = "minimumBaseQuality";
    protected static final String MINIMUM_BASE_QUALITY_SHORT_NAME = "minBQ";

    protected static final String HET_CALLING_STRINGENCY_FULL_NAME = "hetCallingStringency";
    protected static final String HET_CALLING_STRINGENCY_SHORT_NAME = "hetS";

    protected static final String QUADRATURE_ORDER_FULL_NAME = "quadratureOrder";
    protected static final String QUADRATURE_ORDER_SHORT_NAME = "quad";

    protected static final String MINIMUM_ABNORMAL_FRACTION_FULL_NAME = "minimumAbnormalFraction";
    protected static final String MINIMUM_ABNORMAL_FRACTION_SHORT_NAME = "minAF";

    protected static final String MAXIMUM_ABNORMAL_FRACTION_FULL_NAME = "maximumAbnormalFraction";
    protected static final String MAXIMUM_ABNORMAL_FRACTION_SHORT_NAME = "maxAF";

    protected static final String MAXIMUM_COPY_NUMBER_FULL_NAME = "maximumCopyNumber";
    protected static final String MAXIMUM_COPY_NUMBER_SHORT_NAME = "maxCN";

    protected static final String ERROR_PROBABILITY_ADJUSTMENT_FACTOR_FULL_NAME = "errorAdjustmentFactor";
    protected static final String ERROR_PROBABILITY_ADJUSTMENT_FACTOR_SHORT_NAME = "errFact";

    @ArgumentCollection
    protected static final ReferenceInputArgumentCollection REFERENCE_ARGUMENTS =
            new RequiredReferenceInputArgumentCollection();

    @Argument(
            doc = "Normal BAM file (for normal-only and paired normal-tumor cases).",
            fullName = ExomeStandardArgumentDefinitions.NORMAL_BAM_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.NORMAL_BAM_FILE_SHORT_NAME,
            optional = true
    )
    protected File normalBamFile;

    @Argument(
            doc = "Tumor BAM file (for tumor-only and paired normal-tumor cases).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_SHORT_NAME,
            optional = true
    )
    protected File tumorBamFile;

    @Argument(
            doc = "IntervalList file of common SNPs.",
            fullName = ExomeStandardArgumentDefinitions.SNP_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpFile;

    @Argument(
            doc = "Output file for base counts at called heterozygous SNP sites (normal).",
            fullName = ExomeStandardArgumentDefinitions.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.NORMAL_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = true
    )
    protected File normalHetOutputFile;

    @Argument(
            doc = "Output file for base counts at called heterozygous SNP sites (tumor).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = true
    )
    protected File tumorHetOutputFile;

    @Argument(
            doc = "Minimum read depth; SNP sites with lower read depths will be ignored.",
            fullName = READ_DEPTH_THRESHOLD_FULL_NAME,
            shortName = READ_DEPTH_THRESHOLD_SHORT_NAME,
            optional = true
    )
    protected int readDepthThreshold = 15;

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
            doc = "Heterozygosous SNP calling stringency.",
            fullName = HET_CALLING_STRINGENCY_FULL_NAME,
            shortName = HET_CALLING_STRINGENCY_SHORT_NAME,
            optional = false
    )
    protected double hetCallingStringency = 5;

    @Argument(
            doc = "Estimated minimum abnormal cell fraction (for building the allele fraction prior; to be used for " +
                    "tumor-only cases)",
            fullName = MINIMUM_ABNORMAL_FRACTION_FULL_NAME,
            shortName = MINIMUM_ABNORMAL_FRACTION_SHORT_NAME,
            optional = true
    )
    protected double minimumAbnormalFraction = 0.5;

    @Argument(
            doc = "Estimated maximum abnormal cell fraction (for building the allele fraction prior; to be used for " +
                    "tumor-only cases)",
            fullName = MAXIMUM_ABNORMAL_FRACTION_FULL_NAME,
            shortName = MAXIMUM_ABNORMAL_FRACTION_SHORT_NAME,
            optional = true
    )
    protected double maximumAbnormalFraction = 0.8;

    @Argument(
            doc = "Estimated maximum copy number for either homologs in the abnormal portion of the sample (for " +
                    "building the allele fraction prior; to be used for tumor-only cases)",
            fullName = MAXIMUM_COPY_NUMBER_FULL_NAME,
            shortName = MAXIMUM_COPY_NUMBER_SHORT_NAME,
            optional = true
    )
    protected int maximumCopyNumber = 2;

    @Argument(
            doc = "Integration quadrature order; a good choice is the typical read depth (to be used for tumor-only cases)",
            fullName = QUADRATURE_ORDER_FULL_NAME,
            shortName = QUADRATURE_ORDER_SHORT_NAME,
            optional = true
    )
    protected int quadratureOrder = 200;

    @Argument(
            doc = "(Experimental) error probability adjustment factor (could be any positive value)",
            fullName = ERROR_PROBABILITY_ADJUSTMENT_FACTOR_FULL_NAME,
            shortName = ERROR_PROBABILITY_ADJUSTMENT_FACTOR_SHORT_NAME,
            optional = true
    )
    protected double errorProbabilityAdjustmentFactor = 1.0;

    private enum HetCoverageJobType {
        NORMAL_ONLY, TUMOR_ONLY, NORMAL_TUMOR, NONE
    }

    @Override
    protected Object doWork() {

        /* check for illegal arguments */
        if (normalBamFile == null && tumorBamFile == null) {
            throw new UserException.CommandLineException("At least one of normal and tumor BAM files must be provided. Stopping.");
        }
        if (normalBamFile != null && normalHetOutputFile == null) {
            throw new UserException.CommandLineException("If the normal BAM file is provided, the normal Het output file must also be provided. Stopping.");
        }
        if (tumorBamFile != null && tumorHetOutputFile == null) {
            throw new UserException.CommandLineException("If the tumor BAM file is provided, the tumor Het output file must also be provided. Stopping.");
        }

        /* determine the job type */
        HetCoverageJobType jobType = HetCoverageJobType.NONE;
        if (normalBamFile != null && tumorBamFile != null) {
            jobType = HetCoverageJobType.NORMAL_TUMOR;
            logger.info("NORMAL_TUMOR mode selected.");
        }
        if (normalBamFile != null && tumorBamFile == null) {
            jobType = HetCoverageJobType.NORMAL_ONLY;
            logger.info("NORMAL_ONLY mode selected.");
        }
        if (normalBamFile == null && tumorBamFile != null) {
            jobType = HetCoverageJobType.TUMOR_ONLY;
            logger.info("TUMOR_ONLY mode selected.");
        }

        BayesianHetPulldownCalculator normalHetPulldownCalculator, tumorHetPulldownCalculator;
        Pulldown normalHetPulldown, tumorHetPulldown;

        /* run the job */
        switch (jobType) {

            case NORMAL_ONLY:

                normalHetPulldownCalculator = new BayesianHetPulldownCalculator(REFERENCE_ARGUMENTS.getReferenceFile(),
                        IntervalList.fromFile(snpFile), minimumMappingQuality, minimumBaseQuality, readDepthThreshold,
                        VALIDATION_STRINGENCY, errorProbabilityAdjustmentFactor);

                logger.info("Calculating the Het pulldown from the normal BAM file using the BALANCED prior...");
                normalHetPulldown = normalHetPulldownCalculator.getHetPulldown(normalBamFile, hetCallingStringency);

                logger.info("Writing Het pulldown from normal reads to " + normalHetOutputFile.toString());
                normalHetPulldown.write(normalHetOutputFile);

                break;

            case TUMOR_ONLY:

                tumorHetPulldownCalculator = new BayesianHetPulldownCalculator(REFERENCE_ARGUMENTS.getReferenceFile(),
                        IntervalList.fromFile(snpFile), minimumMappingQuality, minimumBaseQuality, readDepthThreshold,
                        VALIDATION_STRINGENCY, errorProbabilityAdjustmentFactor);
                tumorHetPulldownCalculator.useHeterogeneousHetPrior(minimumAbnormalFraction, maximumAbnormalFraction, maximumCopyNumber,
                        quadratureOrder);

                logger.info("Calculating the Het pulldown from the tumor BAM file using the HETEROGENEOUS prior...");
                tumorHetPulldown = tumorHetPulldownCalculator.getHetPulldown(tumorBamFile, hetCallingStringency);

                logger.info("Writing Het pulldown from tumor reads to " + tumorHetOutputFile.toString());
                tumorHetPulldown.write(tumorHetOutputFile);

                break;

            case NORMAL_TUMOR:

                normalHetPulldownCalculator = new BayesianHetPulldownCalculator(REFERENCE_ARGUMENTS.getReferenceFile(),
                        IntervalList.fromFile(snpFile), minimumMappingQuality, minimumBaseQuality, readDepthThreshold,
                        VALIDATION_STRINGENCY, errorProbabilityAdjustmentFactor);

                logger.info("Calculating the Het pulldown from the normal BAM file using the BALANCED prior...");
                normalHetPulldown = normalHetPulldownCalculator.getHetPulldown(normalBamFile, hetCallingStringency);

                logger.info("Writing Het pulldown from normal reads to " + normalHetOutputFile.toString());
                normalHetPulldown.write(normalHetOutputFile);

                tumorHetPulldownCalculator = new BayesianHetPulldownCalculator(REFERENCE_ARGUMENTS.getReferenceFile(),
                        normalHetPulldown.getIntervals(), minimumMappingQuality, minimumBaseQuality, readDepthThreshold,
                        VALIDATION_STRINGENCY, errorProbabilityAdjustmentFactor);

                logger.info("Calculating the Het pulldown from the tumor BAM file on Hets detected in the normal BAM file...");
                tumorHetPulldown = tumorHetPulldownCalculator.getTumorHetPulldownFromNormalPulldown(tumorBamFile,
                        normalHetPulldown);

                logger.info("Writing Het pulldown from tumor reads to " + tumorHetOutputFile.toString());
                tumorHetPulldown.write(tumorHetOutputFile);

                break;

            case NONE:
                break;

            default:
                break;

        }

        return "SUCCESS";
    }
}