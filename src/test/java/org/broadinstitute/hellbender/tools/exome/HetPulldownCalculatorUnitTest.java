package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.EnumMap;
import java.util.Map;

/**
 * Unit tests for {@link HetPulldownCalculator}.  Uses BAM and SNP files generated from hg19mini using wgsim.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class HetPulldownCalculatorUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    private static final File NORMAL_BAM_FILE = new File(TEST_SUB_DIR + "normal.sorted.bam");
    private static final File NORMAL_UNSORTED_BAM_FILE = new File(TEST_SUB_DIR + "normal.unsorted.bam");
    private static final File TUMOR_BAM_FILE = new File(TEST_SUB_DIR + "tumor.sorted.bam");
    private static final File SNP_FILE = new File(TEST_SUB_DIR + "common_SNP.interval_list");
    private static final File REF_FILE = new File(hg19MiniReference);

    private static final HetPulldownCalculator calculator = new HetPulldownCalculator(REF_FILE, SNP_FILE);

    private static SAMFileHeader normalHeader;
    private static SAMFileHeader tumorHeader;

    @BeforeClass
    public void initHeaders() throws IOException {
        try (final SamReader normalBamReader = SamReaderFactory.makeDefault().open(NORMAL_BAM_FILE);
             final SamReader tumorBamReader = SamReaderFactory.makeDefault().open(TUMOR_BAM_FILE)) {
            normalHeader = normalBamReader.getFileHeader();
            tumorHeader = tumorBamReader.getFileHeader();
        }
    }

    private static EnumMap<Nucleotide, Integer> makeBaseCounts(final int aCount, final int cCount,
                                                               final int gCount, final int tCount) {
        final EnumMap<Nucleotide, Integer> baseCounts = new EnumMap<>(Nucleotide.class);
        baseCounts.put(Nucleotide.A, aCount);
        baseCounts.put(Nucleotide.C, cCount);
        baseCounts.put(Nucleotide.G, gCount);
        baseCounts.put(Nucleotide.T, tCount);
        return baseCounts;
    }

    private static void assertEqualBaseCounts(final Nucleotide.Counter actual, final EnumMap<Nucleotide, Integer> expected) {
        for (final Nucleotide base : HetPulldownCalculator.BASES) {
            Assert.assertEquals(actual.get(base), (long) expected.get(base));
        }
    }

    @DataProvider(name = "inputGetPileupBaseCount")
    public Object[][] inputGetPileupBaseCount() throws IOException {
        try (final SamReader bamReader = SamReaderFactory.makeDefault().open(NORMAL_BAM_FILE)) {
            final IntervalList intervals = new IntervalList(bamReader.getFileHeader());
            intervals.add(new Interval("1", 100, 100));
            intervals.add(new Interval("1", 11000, 11000));
            intervals.add(new Interval("1", 14000, 14000));
            intervals.add(new Interval("1", 14630, 14630));

            final SamLocusIterator locusIterator = new SamLocusIterator(bamReader, intervals);

            final EnumMap<Nucleotide, Integer> baseCounts1 = makeBaseCounts(0, 0, 0, 0);
            final EnumMap<Nucleotide, Integer> baseCounts2 = makeBaseCounts(0, 9, 0, 0);
            final EnumMap<Nucleotide, Integer> baseCounts3 = makeBaseCounts(12, 0, 0, 0);
            final EnumMap<Nucleotide, Integer> baseCounts4 = makeBaseCounts(0, 0, 8, 9);

            if (!locusIterator.hasNext()) {
                throw new SAMException("Can't get locus to start iteration. Check that " + NORMAL_BAM_FILE.toString()
                        + " contains 1:0-16000.");
            }
            final SamLocusIterator.LocusInfo locus1 = locusIterator.next();
            final SamLocusIterator.LocusInfo locus2 = locusIterator.next();
            final SamLocusIterator.LocusInfo locus3 = locusIterator.next();
            final SamLocusIterator.LocusInfo locus4 = locusIterator.next();
            locusIterator.close();

            return new Object[][]{
                    {locus1, baseCounts1},
                    {locus2, baseCounts2},
                    {locus3, baseCounts3},
                    {locus4, baseCounts4}
            };
        }
    }

    @Test(dataProvider = "inputGetPileupBaseCount")
    public void testGetPileupBaseCount(final SamLocusIterator.LocusInfo locus,
                                       EnumMap<Nucleotide, Integer> expected) {
        final Nucleotide.Counter result = HetPulldownCalculator.getPileupBaseCounts(locus);
        assertEqualBaseCounts(result, expected);
    }

    @DataProvider(name = "inputIsPileupHetCompatible")
    public Object[][] inputIsPileupHetCompatible() {
        final double normalErrorRate = 0.01;
        final double equalOddsThreshold = 1.0;

        //if likelihood ratio < likelihoodRatioThreshold, expected = false
        return new Object[][]{
                {50, 50, normalErrorRate, 1e50, true}, //likelihood ratio = 1.3 x 10^70
                {10, 10, normalErrorRate, 1e10, true}, //likelihood ratio = 1.0 x 10^14
                {99, 1, normalErrorRate, 1e-20, false},  // likelihood ratio = 2.1 x 10^-28
                {90, 10, normalErrorRate, equalOddsThreshold, false},  // likelihood ratio = 1.9 x 10^-10
                {86, 14, normalErrorRate, equalOddsThreshold, false},  // likelihood ratio = 0.019
                {85, 15, normalErrorRate, equalOddsThreshold, true},  // likelihood ratio = 1.85
                {10, 1, normalErrorRate, equalOddsThreshold, false},   // likelihood ratio = 0.054
                {10, 2, normalErrorRate, equalOddsThreshold, true},   // likelihood ratio = 5.40
                {10, 2, 0.1, equalOddsThreshold, false}   // likelihood ratio = 0.14
        };
    }

    @Test(dataProvider = "inputIsPileupHetCompatible")
    public void testIsHet(final int refBaseCount, final int altBaseCount,
                                          final double errorRate, final double likelihoodRatioThreshold,
                                          final boolean expected) {
        final boolean result = HetPulldownCalculator.isHet(refBaseCount, altBaseCount,
                errorRate, likelihoodRatioThreshold);
        Assert.assertEquals(result, expected);
    }

    @DataProvider(name = "inputGetNormalHetPulldown")
    public Object[][] inputGetNormalHetPulldown() {
        final Pulldown normalHetPulldown1 = new Pulldown(normalHeader);
        normalHetPulldown1.add(new SimpleInterval("1", 10736, 10736), 9, 2);
        normalHetPulldown1.add(new SimpleInterval("1", 11522, 11522), 7, 4);
        normalHetPulldown1.add(new SimpleInterval("1", 12098, 12098), 8, 6);
        normalHetPulldown1.add(new SimpleInterval("1", 14630, 14630), 9, 8);
        normalHetPulldown1.add(new SimpleInterval("2", 14689, 14689), 6, 9);
        normalHetPulldown1.add(new SimpleInterval("2", 14982, 14982), 6, 5);

        //changing error rate from 0.01 to 0.1 removes first het SNP
        final Pulldown normalHetPulldown2 = new Pulldown(normalHeader);
        normalHetPulldown2.add(new SimpleInterval("1", 11522, 11522), 7, 4);
        normalHetPulldown2.add(new SimpleInterval("1", 12098, 12098), 8, 6);
        normalHetPulldown2.add(new SimpleInterval("1", 14630, 14630), 9, 8);
        normalHetPulldown2.add(new SimpleInterval("2", 14689, 14689), 6, 9);
        normalHetPulldown2.add(new SimpleInterval("2", 14982, 14982), 6, 5);

        return new Object[][]{
                {0.01, 1.0, normalHetPulldown1},
                {0.1, 1.0, normalHetPulldown2}
        };
    }

    @Test(dataProvider = "inputGetNormalHetPulldown")
    public void testGetNormalHetPulldown(final double errorRate, final double likelihoodRatioThreshold,
                                         final Pulldown expected) {
        final Pulldown result = calculator.getNormal(NORMAL_BAM_FILE, errorRate, likelihoodRatioThreshold);
        Assert.assertEquals(result, expected);
    }

    @Test(expectedExceptions = UserException.class)
    public void testGetHetPulldownWithUnsortedBAMFile() {
        final Pulldown result = calculator.getNormal(NORMAL_UNSORTED_BAM_FILE, -1, -1);
    }

    @DataProvider(name = "inputGetTumorHetPulldown")
    public Object[][] inputGetTumorHetPulldown() {
        //first het SNP in normalHetPulldown1 has <= 10 reads in tumor BAM and hence does not pass read depth filter
        final Pulldown tumorHetPulldown = new Pulldown(normalHeader);
        tumorHetPulldown.add(new SimpleInterval("1", 11522, 11522), 7, 4);
        tumorHetPulldown.add(new SimpleInterval("1", 12098, 12098), 8, 6);
        tumorHetPulldown.add(new SimpleInterval("1", 14630, 14630), 9, 8);
        tumorHetPulldown.add(new SimpleInterval("2", 14689, 14689), 6, 9);
        tumorHetPulldown.add(new SimpleInterval("2", 14982, 14982), 6, 5);

        final IntervalList normalHetIntervals = new IntervalList(tumorHeader);
        normalHetIntervals.add(new Interval("1", 11522, 11522));
        normalHetIntervals.add(new Interval("1", 12098, 12098));
        normalHetIntervals.add(new Interval("1", 14630, 14630));
        normalHetIntervals.add(new Interval("2", 14689, 14689));
        normalHetIntervals.add(new Interval("2", 14982, 14982));

        return new Object[][]{
                {normalHetIntervals, tumorHetPulldown}
        };
    }

    @Test(dataProvider = "inputGetTumorHetPulldown")
    public void testGetTumorHetPulldown(final IntervalList normalHetIntervals,
                                        final Pulldown expected) {
        final Pulldown result = calculator.getTumor(TUMOR_BAM_FILE, normalHetIntervals);
        Assert.assertEquals(result, expected);
    }
}
