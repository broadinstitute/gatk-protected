package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Integration test for {@link CallSegments}.
 *
 * @author David Benjamin
 */
public class CallSegmentsTest extends CommandLineProgramTest{
    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/tools/exome/caller");
    private static final File TEST_TARGETS = new File(TEST_DIR,"targets.tsv");
    private static final File TEST_SEGMENTS = new File(TEST_DIR,"segments.tsv");

    @Test
    public void testCallSegments() throws IOException {
        final File outputFile = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + CallSegments.SEGFILE_SHORT_NAME, TEST_SEGMENTS.getAbsolutePath(),
                "-" + CallSegments.TARGET_FILE_SHORT_NAME, TEST_TARGETS.getAbsolutePath(),
                "-" + CallSegments.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
        };
        runCommandLine(arguments);

        HashedListExonCollection<TargetCoverage> targets =
                new HashedListExonCollection<TargetCoverage>(TargetCoverageUtils.readTargetsWithCoverage(TEST_TARGETS));

        List<Segment> segments = SegmentUtils.readCalledSegments(outputFile, targets);

        Assert.assertEquals(segments.get(0).getCall(), "+");
        Assert.assertEquals(segments.get(1).getCall(), "-");
        Assert.assertEquals(segments.get(2).getCall(), "0");
        Assert.assertEquals(segments.get(3).getCall(), "0");
    }
}
