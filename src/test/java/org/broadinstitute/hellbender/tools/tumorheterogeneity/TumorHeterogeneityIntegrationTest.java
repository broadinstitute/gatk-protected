package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.*;

/**
 * Integration tests for {@link TumorHeterogeneity}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityIntegrationTest extends CommandLineProgramTest {
//    private static final List<File> ACNV_SEGMENT_FILES = Arrays.asList(
//            new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/purity-0.8/total_segments-log2cr_sd-0.05-maf_sd-0.005.acnv.seg")
//    );
    private static final List<File> ACNV_SEGMENT_FILES = Arrays.asList(
            new File("/home/slee/working/ipython/purity-ploidy/integration-test/2_clone/purity-1.0/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg")
    );
//    private static final List<File> ACNV_SEGMENT_FILES = Arrays.asList(
//            new File("/home/slee/working/ipython/purity-ploidy/purity-series/0-0-SM-74NEG-sim-final.seg"),
//            new File("/home/slee/working/ipython/purity-ploidy/purity-series/1-10-SM-74P2T-sim-final.seg"),
//            new File("/home/slee/working/ipython/purity-ploidy/purity-series/2-30-SM-74P35-sim-final.seg"),
//            new File("/home/slee/working/ipython/purity-ploidy/purity-series/3-40-SM-74P3J-sim-final.seg"),
//            new File("/home/slee/working/ipython/purity-ploidy/purity-series/4-60-SM-74P3M-sim-final.seg"),
//            new File("/home/slee/working/ipython/purity-ploidy/purity-series/5-70-SM-74P3K-sim-final.seg"),
//            new File("/home/slee/working/ipython/purity-ploidy/purity-series/6-80-SM-74P51-sim-final.seg"),
//            new File("/home/slee/working/ipython/purity-ploidy/purity-series/7-90-SM-74P56-sim-final.seg"),
//            new File("/home/slee/working/ipython/purity-ploidy/purity-series/8-100-SM-74P4M-sim-final.seg")
//    );

    @Test
    public void testTumorHeterogeneity() {
        for (final File ACNV_SEGMENT_FILE : ACNV_SEGMENT_FILES) {
            final String OUTPUT_PREFIX = ACNV_SEGMENT_FILE.getAbsolutePath().replace(".acnv.seg", "");
            final String[] arguments = {
                    "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, ACNV_SEGMENT_FILE.getAbsolutePath(),
                    "--" + TumorHeterogeneity.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX,
//                    "--" + TumorHeterogeneity.DO_FILTER_SEGMENTS_LONG_NAME,
//                    "--" + TumorHeterogeneity.LENGTH_PERCENTILE_LONG_NAME, "5",
//                    "--" + TumorHeterogeneity.LOG2_COPY_RATIO_CREDIBLE_INTERVAL_PERCENTILE_LONG_NAME, "95",
//                    "--" + TumorHeterogeneity.MINOR_ALLELE_FRACTION_CREDIBLE_INTERVAL_PERCENTILE_LONG_NAME, "95",
                    "--" + TumorHeterogeneity.NUM_SAMPLES_CLONAL_LONG_NAME, "2000",
                    "--" + TumorHeterogeneity.NUM_BURN_IN_CLONAL_LONG_NAME, "500",
                    "--" + TumorHeterogeneity.NUM_SAMPLES_LONG_NAME, "10000",
                    "--" + TumorHeterogeneity.NUM_BURN_IN_LONG_NAME, "2000",
                    "--" + TumorHeterogeneity.MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME, "6",
                    "--" + TumorHeterogeneity.MAX_ALLELIC_COPY_NUMBER_LONG_NAME, "6",
                    "--" + TumorHeterogeneity.CONCENTRATION_PRIOR_ALPHA_CLONAL_LONG_NAME, "1E1",
                    "--" + TumorHeterogeneity.CONCENTRATION_PRIOR_BETA_CLONAL_LONG_NAME, "1E1",
                    "--" + TumorHeterogeneity.CONCENTRATION_PRIOR_BETA_LONG_NAME, "1E1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NOISE_FACTOR_PRIOR_ALPHA_LONG_NAME, "25",
                    "--" + TumorHeterogeneity.COPY_RATIO_NOISE_FACTOR_PRIOR_BETA_LONG_NAME, "0.5",
                    "--" + TumorHeterogeneity.MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_ALPHA_LONG_NAME, "25",
                    "--" + TumorHeterogeneity.MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_BETA_LONG_NAME, "0.5",
                    "--" + TumorHeterogeneity.PLOIDY_STATE_PRIOR_COMPLETE_DELETION_PENALTY_LONG_NAME, "1E-1",
                    "--" + TumorHeterogeneity.PLOIDY_STATE_PRIOR_CHANGE_PENALTY_LONG_NAME, "1E-1",
                    "--verbosity", "DEBUG"
            };
            runCommandLine(arguments);
        }
    }
}