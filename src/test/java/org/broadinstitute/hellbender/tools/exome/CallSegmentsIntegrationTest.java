package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

/**
 * Integration test for {@link CallSegments}.
 */
public final class CallSegmentsIntegrationTest extends CommandLineProgramTest{
    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/tools/exome/caller");
    private static final File TEST_TARGETS = new File(TEST_DIR,"targets.tsv");
    private static final File TEST_TARGETS_LEGACY = new File(TEST_DIR,"targets_legacy.tsv");
    private static final File TEST_SEGMENTS = new File(TEST_DIR,"segments.tsv");
    private static final String SAMPLE_NAME = "sample";

    @Test
    public void testCallSegmentsLegacy() {
        final File outputFile = createTempFile("test_call_legacy",".txt");

        final String[] arguments = {
                "-" + CallSegments.SEGFILE_SHORT_NAME, TEST_SEGMENTS.getAbsolutePath(),
                "-" + CallSegments.TARGET_FILE_SHORT_NAME, TEST_TARGETS_LEGACY.getAbsolutePath(),
                "-" + CallSegments.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "--" + CallSegments.SAMPLE_LONG_NAME, SAMPLE_NAME,
                "-" + CallSegments.Z_THRESHOLD_SHORT_NAME, Double.toString(ReCapSegCaller.DEFAULT_Z_SCORE_THRESHOLD),
                "-" + CallSegments.LEGACY_FORMAT_SHORT_NAME, "True"
        };
        runCommandLine(arguments);

        final TargetCollection<TargetCoverage> targets = TargetCoverageUtils.readModeledTargetFileIntoTargetCollection(TEST_TARGETS_LEGACY);

        final List<ModeledSegment> calls = SegmentUtils.readModeledSegmentsFromSegfile(outputFile);

        Assert.assertEquals(calls.get(0).getCall(), "+");
        Assert.assertEquals(calls.get(1).getCall(), "-");
        Assert.assertEquals(calls.get(2).getCall(), "0");
        Assert.assertEquals(calls.get(3).getCall(), "0");
    }

    @Test
    public void testCallSegments() {
        final File outputFile = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + CallSegments.SEGFILE_SHORT_NAME, TEST_SEGMENTS.getAbsolutePath(),
                "-" + CallSegments.TARGET_FILE_SHORT_NAME, TEST_TARGETS.getAbsolutePath(),
                "-" + CallSegments.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "--" + CallSegments.SAMPLE_LONG_NAME, SAMPLE_NAME,
                "-" + CallSegments.Z_THRESHOLD_SHORT_NAME, Double.toString(ReCapSegCaller.DEFAULT_Z_SCORE_THRESHOLD)
        };
        runCommandLine(arguments);

        final TargetCollection<TargetCoverage> targets = TargetCoverageUtils.readReadCountCollectionIntoTargetCollection(TEST_TARGETS);

        final List<ModeledSegment> calls = SegmentUtils.readModeledSegmentsFromSegfile(outputFile);

        Assert.assertEquals(calls.get(0).getCall(), "+");
        Assert.assertEquals(calls.get(1).getCall(), "-");
        Assert.assertEquals(calls.get(2).getCall(), "0");
        Assert.assertEquals(calls.get(3).getCall(), "0");
    }
}
