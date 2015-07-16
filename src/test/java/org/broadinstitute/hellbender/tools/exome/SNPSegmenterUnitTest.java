package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link SNPSegmenter}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SNPSegmenterUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/tools/exome/";

    //Tests that segments are correctly determined using allelic counts from SNP sites.
    @Test
    public void testAllelicFractionBasedSegmentation() throws IOException {
        final double allelicFractionSkew = 1.;
        final float minLogValue = -10.f;
        final String sampleName = "TCGA-02-0001-01C-01D-0182-01";

        final File snpFile = new File(TEST_SUB_DIR + "snps.simplified_for_allelic_fraction_segmentation.tsv");
        final List<AllelicCount> snpCounts = new Pulldown(snpFile, hg19Header).asList();

        final File resultFile = createTempFile("snp-segmenter-test-result", ".seg");
        final File segmentFile = SNPSegmenter.findSegments(snpCounts, allelicFractionSkew, sampleName,
                resultFile, minLogValue);

        final File expectedFile = new File(TEST_SUB_DIR + "snp-segmenter-test-expected.seg");

        Assert.assertTrue(segmentFile.exists(), "SNPSegmenterTest output was not written to temp file: " + resultFile);
        assertEqualSegments(resultFile, expectedFile);
    }

    //Tests that allelic counts are correctly transformed to target coverages.
    @Test
    public void testTransformAllelicCounts() throws IOException {
        final double allelicFractionSkew = 0.96;

        final File snpFile = new File(TEST_SUB_DIR + "snps.simplified_for_allelic_fraction_transformation.tsv");
        final List<AllelicCount> snpCounts = new Pulldown(snpFile, hg19Header).asList();

        final List<TargetCoverage> resultTargets = snpCounts.stream()
                .map(count -> count.asTargetCoverage("snp-target", allelicFractionSkew)).collect(Collectors.toList());

        final List<TargetCoverage> expectedTargets = Arrays.asList(
            new TargetCoverage("snp-target", new SimpleInterval("1", 212360, 212360), 0.3559124087591240),
            new TargetCoverage("snp-target", new SimpleInterval("1", 241501, 241501), 0.1397938144329896),
            new TargetCoverage("snp-target", new SimpleInterval("1", 242173, 242173), 0.0372413793103448),
            new TargetCoverage("snp-target", new SimpleInterval("1", 256641, 256641), 0.3381818181818182),
            new TargetCoverage("snp-target", new SimpleInterval("1", 261164, 261164), 0.4868049792531120),
            new TargetCoverage("snp-target", new SimpleInterval("1", 267204, 267204), 0.1735483870967741),
            new TargetCoverage("snp-target", new SimpleInterval("1", 282282, 282282), 0.3209090909090909),
            new TargetCoverage("snp-target", new SimpleInterval("1", 291649, 291649), 0.1250955414012738),
            new TargetCoverage("snp-target", new SimpleInterval("1", 376402, 376402), 0.2981818181818181),
            new TargetCoverage("snp-target", new SimpleInterval("1", 408347, 408347), 0.2032704402515723),
            new TargetCoverage("snp-target", new SimpleInterval("1", 415813, 415813), 0.1721739130434782),
            new TargetCoverage("snp-target", new SimpleInterval("1", 426517, 426517), 0.4366666666666666),
            new TargetCoverage("snp-target", new SimpleInterval("1", 429357, 429357), 0.1670588235294118),
            new TargetCoverage("snp-target", new SimpleInterval("1", 455201, 455201), 0.0112499999999999),
            new TargetCoverage("snp-target", new SimpleInterval("1", 466369, 466369), 0.4287179487179487),
            new TargetCoverage("snp-target", new SimpleInterval("1", 545461, 545461), 0.1164912280701754),
            new TargetCoverage("snp-target", new SimpleInterval("1", 665716, 665716), 0.2191304347826086),
            new TargetCoverage("snp-target", new SimpleInterval("1", 679370, 679370), 0.0056115107913669),
            new TargetCoverage("snp-target", new SimpleInterval("1", 704935, 704935), 0.4597979797979797));

        Assert.assertEquals(resultTargets, expectedTargets);
    }

    /**
     * Compares the content of two segmenter output files.
     * @param actualOutput the actual segmenter output containing file.
     * @param expectedOutput the expected segmenter output containing file.
     * @throws NullPointerException if either {@code actualOutput} or {@code expectedOutput} is {@code null}.
     * @throws IOException if any was thrown when reading the input files.
     * @throws AssertionError if there are significant between both files.
     */
    private void assertEqualSegments(final File actualOutput, final File expectedOutput) throws IOException {
        try (final SegmentReader actual = new SegmentReader(actualOutput);
             final SegmentReader expected = new SegmentReader(expectedOutput)) {
            final List<SegmentMean> actualSegmentMeans = actual.stream().collect(Collectors.toList());
            final List<SegmentMean> expectedSegmentMeans = expected.stream().collect(Collectors.toList());
            Assert.assertEquals(actualSegmentMeans, expectedSegmentMeans);
        }
    }

    /**
     * Represent a Segment in the segmenter output.
     */
    private static class SegmentMean {
        public final double segmentMean;

        public SegmentMean(final double segmentMean) {
            this.segmentMean = segmentMean;
        }

        @Override
        public boolean equals(final Object other) {
            if (!(other instanceof SegmentMean)) {
                return false;
            }

            final SegmentMean otherSegmentMean = (SegmentMean) other;
            return Math.abs(otherSegmentMean.segmentMean - segmentMean) < 0.0000000000001;
        }

        @Override
        public int hashCode() {
            return Double.hashCode(segmentMean);
        }
    }

    /**
     * Tsv reader for the Segmenter output.
     */
    private static class SegmentReader extends TableReader<SegmentMean> {
        public SegmentReader(final File file) throws IOException {
            super(file);
        }

        @Override
        protected SegmentMean createRecord(DataLine dataLine) {
            return new SegmentMean(dataLine.getDouble("Segment_Mean"));
        }
    }
}
