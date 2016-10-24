package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.testng.annotations.Test;

import java.io.File;

import static org.testng.Assert.*;

/**
 * Integration tests for {@link TumorHeterogeneity}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityIntegrationTest extends CommandLineProgramTest {
//    private static final File ACNV_SEGMENT_FILE = new File("/home/slee/working/ipython/purity-ploidy/integration-test/purity-1.0/total_segments-log2cr_sd-0.05-maf_sd-0.005.acnv.seg");
    private static final File ACNV_SEGMENT_FILE = new File("/home/slee/working/ipython/purity-ploidy/integration-test/purity-1.0/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg");
    private static final String OUTPUT_PREFIX = ACNV_SEGMENT_FILE.getAbsolutePath().replace(".acnv.seg", "");

    @Test
    public void testTumorHeterogeneity() {
        final String[] arguments = {
                "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, ACNV_SEGMENT_FILE.getAbsolutePath(),
                "--" + TumorHeterogeneity.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX,
                "--" + TumorHeterogeneity.NUM_SAMPLES_CLONAL_LONG_NAME, "200",
                "--" + TumorHeterogeneity.NUM_BURN_IN_CLONAL_LONG_NAME, "150",
                "--" + TumorHeterogeneity.NUM_SAMPLES_LONG_NAME, "500",
                "--" + TumorHeterogeneity.NUM_BURN_IN_LONG_NAME, "400",
                "--" + TumorHeterogeneity.NUM_CELLS_LONG_NAME, "50",
                "--" + TumorHeterogeneity.METROPOLIS_ITERATION_FRACTION_CLONAL_LONG_NAME, "0.25",
                "--" + TumorHeterogeneity.METROPOLIS_ITERATION_FRACTION_LONG_NAME, "0.25",
                "--" + TumorHeterogeneity.CONCENTRATION_PRIOR_BETA_LONG_NAME, "1E1",
                "--" + TumorHeterogeneity.PLOIDY_STATE_PRIOR_COMPLETE_DELETION_PENALTY_LONG_NAME, "1E-3",
                "--" + TumorHeterogeneity.PLOIDY_STATE_PRIOR_CHANGE_PENALTY_LONG_NAME, "1E-3",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }
}