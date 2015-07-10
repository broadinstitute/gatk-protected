package org.broadinstitute.hellbender.tools.exome;

//class SNPSegmenterTest(unittest.TestCase):

//
//        def test_transform_allelic_fractions(self):
//        """
//        Tests that the allelic fractions are transformed correctly.
//        """
//        snps_filename = os.path.join(self._testdata_directory, "snps.simplified_for_allelic_fraction_transformation.tsv")
//        dataframe = TestUtils.load_snps_dataframe(snps_filename)
//
//        reference_counts = dataframe["t_ref_count"]
//        alternate_counts = dataframe["t_alt_count"]
//        allelic_fraction_skew = 0.96
//        transformed_allelic_fractions = \
//        SNPSegmenter._transform_allelic_fractions(reference_counts=reference_counts,
//        alternate_counts=alternate_counts,
//        allelic_fraction_skew=allelic_fraction_skew)
//        gt_transformed_allelic_fractions = pandas.Series([0.355912408759124054835609740621293894946575164794921875,
//        0.13979381443298966036792307932046242058277130126953125,
//        0.03724137931034487980497260650736279785633087158203125,
//        0.33818181818181825004643314969143830239772796630859375,
//        0.48680497925311205786869095391011796891689300537109375,
//        0.173548387096774170412771809424157254397869110107421875,
//        0.32090909090909092160615045941085554659366607666015625,
//        0.12509554140127387977798889551195316016674041748046875,
//        0.298181818181818159008145130428601987659931182861328125,
//        0.203270440251572315215611297389841638505458831787109375,
//        0.1721739130434782882872468690038658678531646728515625,
//        0.43666666666666664742280090649728663265705108642578125,
//        0.16705882352941181512306911827181465923786163330078125,
//        0.011249999999999982236431605997495353221893310546875,
//        0.428717948717948715842140927634318359196186065673828125,
//        0.11649122807017542324814485255046747624874114990234375,
//        0.219130434782608685129190462248516269028186798095703125,
//        0.00561151079136690267290532574406825006008148193359375,
//        0.45979797979797976115179380940389819443225860595703125])
//
//        transformed_allelic_fractions_tolerance = 1e-10
//        for index, transformed_allelic_fraction in enumerate(transformed_allelic_fractions):
//        gt_transformed_allelic_fraction = gt_transformed_allelic_fractions[index]
//        self.assertLess(TestUtils.determine_percent_error(measured_value=transformed_allelic_fraction,
//        actual_value=gt_transformed_allelic_fraction),
//        transformed_allelic_fractions_tolerance,
//        "Transformed allelic fraction must be %s but it was %s."
//        % (gt_transformed_allelic_fractions[index], transformed_allelic_fraction))

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link SNPSegmenter}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SNPSegmenterUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/tools/exome/";

    //Tests that segments are correctly determined using allelic-count data from SNPs.
    @Test
    public void testAllelicFractionBasedSegmentation() throws IOException {
        final double allelicFractionSkew = 1.;
        final float minLogValue = -10.f;
        final String sampleName = "TCGA-02-0001-01C-01D-0182-01";
        final File snpFile = new File(TEST_SUB_DIR + "snps.simplified_for_allelic_fraction_transformation.tsv");
        final File expectedFile = new File(TEST_SUB_DIR + "snp-segmenter-test-expected.seg");

        final List<AllelicCount> snpCounts = new Pulldown(snpFile, hg19Header).asList();

        final File resultFile = createTempFile("snp-segmenter-test-result", ".seg");
        final File segmentFile = SNPSegmenter.findSegments(snpCounts, allelicFractionSkew, sampleName,
                resultFile, minLogValue);

        Assert.assertTrue(segmentFile.exists(), "SNPSegmenterTest output was not written to temp file: " + resultFile);
        assertEqualSegments(resultFile, expectedFile);
    }

    /**
     * Compares the content of two segmenter output files.  Borrowed from SegmenterTest, extract later?
     * @param actualOutput the actual segmenter output containing file.
     * @param expectedOutput the expected segmenter output containing file.
     * @throws NullPointerException if either {@code actualOutput} or {@code expectedOutput} is {@code null}.
     * @throws IOException if any was thrown when reading the input files.
     * @throws AssertionError if there are significant between both files.
     */
    private void assertEqualSegments(final File actualOutput, final File expectedOutput) throws IOException {
        try (final SegmentReader actual = new SegmentReader(actualOutput);
             final SegmentReader expected = new SegmentReader(expectedOutput)) {
            final List<Segment> actualSegments = actual.stream().collect(Collectors.toList());
            final List<Segment> expectedSegments = expected.stream().collect(Collectors.toList());
            Assert.assertEquals(actualSegments, expectedSegments);
        }
    }

    /**
     * Represent a Segment in the segmenter output.  Borrowed from SegmenterTest, temporary until structure of
     * ModelSegment class decided on.
     */
    private static class Segment {
        public final double segmentMean;

        public Segment(final double segmentMean) {
            this.segmentMean = segmentMean;
        }

        @Override
        public boolean equals(final Object other) {
            if (!(other instanceof Segment))
                return false;
            final Segment otherSegment = (Segment) other;
            return Math.abs(otherSegment.segmentMean - segmentMean) < 0.0000000000001;
        }

        @Override
        public int hashCode() {
            return Double.hashCode(segmentMean);
        }
    }

    /**
     * Tsv reader for the Segmenter output.  Borrowed from SegmenterTest, temporary until structure of
     * ModelSegment class decided on.
     */
    private static class SegmentReader extends TableReader<Segment> {
        public SegmentReader(final File file) throws IOException {
            super(file);
        }

        @Override
        protected Segment createRecord(DataLine dataLine) {
            return new Segment(dataLine.getDouble("Segment_Mean"));
        }
    }
}
