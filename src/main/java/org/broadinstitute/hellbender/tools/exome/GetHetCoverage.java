package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;

/**
 * Outputs reference/alternate read counts for a tumor sample at heterozygous SNP sites present in a normal sample.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Output ref/alt counts for tumor sample at heterozygous SNPs in normal sample.",
        oneLineSummary = "Output ref/alt counts for tumor sample at heterozygous SNPs in normal sample.",
        programGroup = CopyNumberProgramGroup.class
)
public final class GetHetCoverage extends CommandLineProgram {

    @ArgumentCollection
    protected static final ReferenceInputArgumentCollection REF_ARGUMENTS =
            new RequiredReferenceInputArgumentCollection();

    protected static final String ERROR_RATE_FULL_NAME = "errorRate";
    protected static final String ERROR_RATE_SHORT_NAME = "ER";

    protected static final String LIKELIHOOD_RATIO_THRESHOLD_FULL_NAME = "likelihoodRatioThreshold";
    protected static final String LIKELIHOOD_RATIO_THRESHOLD_SHORT_NAME = "LR";

    private static final double DEFAULT_LIKELIHOOD_RATIO_THRESHOLD = 1000;
    private static final double DEFAULT_ERROR_RATE = 0.01;

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
            optional = false
    )
    protected File tumorBAMFile;

    @Argument(
            doc = "Interval-list file of common SNPs.",
            fullName = ExomeStandardArgumentDefinitions.SNP_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpFile;

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
            optional = false
    )
    protected File tumorHetOutputFile;

    @Argument(
            doc = "Estimated sequencing error rate.",
            fullName = ERROR_RATE_FULL_NAME,
            shortName = ERROR_RATE_SHORT_NAME,
            optional = true
    )
    protected double errorRate = DEFAULT_ERROR_RATE;

    @Argument(
            doc = "het:hom likelihood ratio threshold for calling heterozygous SNPs in normal sample.",
            fullName = LIKELIHOOD_RATIO_THRESHOLD_FULL_NAME,
            shortName = LIKELIHOOD_RATIO_THRESHOLD_SHORT_NAME,
            optional = true
    )
    protected double likelihoodRatioThreshold = DEFAULT_LIKELIHOOD_RATIO_THRESHOLD;

    @Override
    protected Object doWork() {
        if (likelihoodRatioThreshold < 0) {
            throw new UserException.BadArgumentValue(LIKELIHOOD_RATIO_THRESHOLD_FULL_NAME,
                    Double.toString(likelihoodRatioThreshold),
                    "Likelihood ratio threshold may not be negative.");
        }

        if (errorRate < 0 || errorRate > 1) {
            throw new UserException.BadArgumentValue(ERROR_RATE_FULL_NAME,
                    Double.toString(errorRate),
                    "Error rate threshold must be between 0 and 1.");
        }

        final HetPulldownCalculator hetPulldown = new HetPulldownCalculator(REF_ARGUMENTS.getReferenceFile(), snpFile);

        logger.info("Getting normal het pulldown...");
        final Pulldown normalHetPulldown = hetPulldown.getNormal(normalBAMFile, errorRate, likelihoodRatioThreshold);
        normalHetPulldown.write(normalHetOutputFile);
        logger.info("Normal het pulldown written to " + normalHetOutputFile.toString());

        final IntervalList normalHetIntervals = normalHetPulldown.getIntervals();

        logger.info("Getting tumor het pulldown...");
        final Pulldown tumorHetPulldown = hetPulldown.getTumor(tumorBAMFile, normalHetIntervals);
        tumorHetPulldown.write(tumorHetOutputFile);
        logger.info("Tumor het pulldown written to " + tumorHetOutputFile.toString());

        return "SUCCESS";
    }
}