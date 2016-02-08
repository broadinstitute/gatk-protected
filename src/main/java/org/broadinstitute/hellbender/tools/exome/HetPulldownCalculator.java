package org.broadinstitute.hellbender.tools.exome;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.NotPrimaryAlignmentFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Gets heterozygous SNP pulldown for normal and tumor samples.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class HetPulldownCalculator {
    private final Logger logger = LogManager.getLogger(HetPulldownCalculator.class);

    private final File refFile;
    private final IntervalList snpIntervals;

    /** Set quality and read-depth thresholds for pulldown, interval threshold for indexing for SamLocusIterator. */
    private static final int MIN_QUALITY = 0;
    private static final int READ_DEPTH_THRESHOLD = 10;
    private static final int MAX_INTERVALS_FOR_INDEX = 25000;

    @VisibleForTesting
    static final Nucleotide[] BASES = {Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.T};
    
    public HetPulldownCalculator(final File refFile, final File snpFile) {
        this.refFile = refFile;
        this.snpIntervals = IntervalList.fromFile(snpFile);
    }

    /**
     * Provide flags for running getHetPulldown based on sample type (normal or tumor).
     */
    private enum SampleType {
        NORMAL, TUMOR
    }

    /**
     * Returns base-pair counts at a given locus.
     * @param locus locus
     * @return      map of base-pair counts
     */
    @VisibleForTesting
    static Nucleotide.Counter getPileupBaseCounts(final SamLocusIterator.LocusInfo locus) {
        final Nucleotide.Counter result = new Nucleotide.Counter();
        for (final SamLocusIterator.RecordAndOffset rec : locus.getRecordAndPositions()) {
            result.add(rec.getReadBase());
        }
        return result;
    }

    private static Pair<Integer, Integer> refAndAltReadCounts(final Nucleotide.Counter baseCounts,
                                                                 final Nucleotide refBase) {
        final long refReadCount = baseCounts.get(refBase);
        final long altReadCount = Arrays.stream(BASES)
                .filter(base -> base != refBase)
                .mapToLong(baseCounts::get).max().getAsLong();
        return new ImmutablePair<>((int) refReadCount, (int) altReadCount);
    }

    /**
     * Returns true if the distribution of major and other base-pair counts from a pileup at a locus is compatible with
     * a given heterozygous allele fraction.
     *
     * <p>
     *     Compatibility is defined by a likelihood ratio threshold between two hypotheses.  First, that the site is
     *     heterozygous in which case the likelihood is nCa * (1/2)^n.  Second, that the site is homozygous and minor
     *     allele reads are sequencing errors, with likelihood nCa * errorRate^a * (1 - errorRate)^(n-a).
     *     Here n and a are total and alt read counts.
     * </p>
     * @param refReadCount      number of ref reads at this site
     * @param altReadCount      number of reads of most common alt at this site
     * @param errorRate         estimated substitution error rate -- result is not very sensitive to this
     * @param likelihoodRatioThreshold    ratio of het to hom likelihood required to call a het
     * @return                  boolean compatibility with heterozygous allele fraction
     */
    @VisibleForTesting
    static boolean isHet(final int refReadCount, final int altReadCount,
                                final double errorRate, final double likelihoodRatioThreshold) {
        final int minorReadCount = Math.min(refReadCount, altReadCount);
        final int majorReadCount = Math.max(refReadCount, altReadCount);

        //work in log space to avoid underflow
        // ignore nCa combinatorial factors common to het and hom because this cancels in the ratio
        final double hetLogLikelihood = (minorReadCount + majorReadCount) * Math.log(0.5);
        final double homLogLikelihood = minorReadCount * Math.log(errorRate) + majorReadCount * Math.log(1 - errorRate);
        final double logLikelihoodRatio = hetLogLikelihood - homLogLikelihood;
        return logLikelihoodRatio > Math.log(likelihoodRatioThreshold);
    }

    /**
     * Calls {@link HetPulldownCalculator#getHetPulldown} with flags set for a normal sample.
     */
    public Pulldown getNormal(final File normalBAMFile, final double errorRate, final double likelihoodRatioThreshold) {
        return getHetPulldown(normalBAMFile, this.snpIntervals, SampleType.NORMAL, errorRate, likelihoodRatioThreshold);
    }

    /**
     * Calls {@link HetPulldownCalculator#getHetPulldown} with flags set for a tumor sample.
     */
    public Pulldown getTumor(final File tumorBAMFile, final IntervalList normalHetIntervals) {
        return getHetPulldown(tumorBAMFile, normalHetIntervals, SampleType.TUMOR, -1, -1);
    }

    /**
     * For a normal or tumor sample, returns a data structure giving (intervals, reference counts, alternate counts),
     * where intervals give positions of likely heterozygous SNP sites.
     *
     * <p>
     *     For a normal sample:
     *     <ul>
     *         The IntervalList snpIntervals gives common SNP sites in 1-based format.
     *     </ul>
     *     <ul>
     *         The estimated sequencing error rate and likelihood ratio threshold define a likelihood
     *         test for heterozygous SNP sites, given the sample.  Only these sites are output.
     *     </ul>
     * </p>
     * <p>
     *     For a tumor sample:
     *     <ul>
     *         The IntervalList snpIntervals gives heterozygous SNP sites likely to be present in the normal sample.
     *         This should be from {@link HetPulldownCalculator#getNormal} in 1-based format.
     *         Only these sites are output.
     *     </ul>
     * </p>
     * @param bamFile           sorted BAM file for sample
     * @param snpIntervals      IntervalList of SNP sites
     * @param sampleType        flag indicating type of sample (SampleType.NORMAL or SampleType.TUMOR)
     *                          (determines whether to perform binomial test)
     * @param errorRate         estimated substitution error rate of sequencing
     * @param likelihoodRatioThreshold     het:hom likelihood ratio threshold for calling a het
     * @return                  Pulldown of heterozygous SNP sites in 1-based format
     */
    private Pulldown getHetPulldown(final File bamFile, final IntervalList snpIntervals, SampleType sampleType,
                                    final double errorRate, final double likelihoodRatioThreshold) {
        try (final SamReader bamReader = SamReaderFactory.makeDefault().open(bamFile);
             final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(this.refFile)) {
            if (bamReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new UserException.BadInput("BAM file " + bamFile.toString() + " must be coordinate sorted.");
            }

            final Pulldown hetPulldown = new Pulldown(bamReader.getFileHeader());

            final int totalNumberOfSNPs = snpIntervals.size();
            final SamLocusIterator locusIterator = new SamLocusIterator(bamReader, snpIntervals,
                    totalNumberOfSNPs < MAX_INTERVALS_FOR_INDEX);

            //set read and locus filters [note: read counts match IGV, but off by a few from pysam.mpileup]
            final List<SamRecordFilter> samFilters = Arrays.asList(new NotPrimaryAlignmentFilter(),
                    new DuplicateReadFilter());
            locusIterator.setSamFilters(samFilters);
            locusIterator.setEmitUncoveredLoci(false);
            locusIterator.setQualityScoreCutoff(MIN_QUALITY);
            locusIterator.setIncludeNonPfReads(false);

            logger.info("Examining " + totalNumberOfSNPs + " sites...");
            final int iterationsPerStatus = Math.max((int) Math.floor(totalNumberOfSNPs / 20.), 1);
            int locusCount = 1;
            for (final SamLocusIterator.LocusInfo locus : locusIterator) {
                if (locusCount % iterationsPerStatus == 0) {
                    logger.info("Examined " + locusCount + " out of " + totalNumberOfSNPs + " sites.");
                }
                locusCount++;

                //include N, etc. reads here
                final int totalReadCount = locus.getRecordAndPositions().size();
                if (totalReadCount <= READ_DEPTH_THRESHOLD) {
                    continue;
                }

                final Nucleotide.Counter baseCounts = getPileupBaseCounts(locus);
                final Nucleotide refBase = Nucleotide.valueOf(refWalker.get(locus.getSequenceIndex()).getBases()[locus.getPosition() - 1]);
                final Pair<Integer, Integer> refAltCounts = refAndAltReadCounts(baseCounts, refBase);
                final int refReadCount = refAltCounts.getLeft();
                final int altReadCount = refAltCounts.getRight();

                if (sampleType == SampleType.TUMOR || isHet(refReadCount, altReadCount, errorRate, likelihoodRatioThreshold)) {
                    final SimpleInterval interval = new SimpleInterval(locus.getSequenceName(), locus.getPosition(), locus.getPosition());
                    hetPulldown.add(interval, refReadCount, altReadCount);
                }
            }
            logger.info("Examined " + totalNumberOfSNPs + " sites.");
            return hetPulldown;
        } catch (final IOException e) {
            throw new UserException(e.getMessage());
        }
    }
}
