package org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount;

import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Collects reference/alternate allele counts at sites.  The alt count is defined as the total count minus the ref count,
 * and the alt nucleotide is defined as the non-ref base with the highest count, with ties broken by the order of the
 * bases in {@link AllelicCountCollector#BASES}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicCountCollector {

    private static final Logger logger = LogManager.getLogger(AllelicCountCollector.class);

    public static final List<Nucleotide> BASES = Collections.unmodifiableList(Arrays.asList(Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.T));

    private final AllelicCountCollection allelicCounts = new AllelicCountCollection();

    public AllelicCountCollector() {
    }

    //TODO: Populate docs
    public void collectAtLocus(final Nucleotide refBase, final ReadPileup pileup, final Locatable locus, final int minBaseQuality) {

        if (!BASES.contains(refBase)) {
            logger.warn(String.format("The reference position at %s has an unknown base call (value: %s). Skipping...",
                    locus, refBase.toString()));
            return;
        }

        final Nucleotide.Counter nucleotideCounter = new Nucleotide.Counter();

        Utils.stream(pileup.iterator())
                .filter(r -> !r.isDeletion())
                .filter(r -> r.getQual() >= minBaseQuality)
                .forEach(r -> nucleotideCounter.add(r.getBase()));

        final int totalBaseCount = BASES.stream().mapToInt(b -> (int) nucleotideCounter.get(b)).sum(); //only include total ACGT counts in binomial test (exclude N, etc.)
        final int refReadCount = (int) nucleotideCounter.get(refBase);
        final int altReadCount = totalBaseCount - refReadCount;                                 //we take alt = total - ref instead of the actual alt count
        final Nucleotide altBase = inferAltFromPileupBaseCounts(nucleotideCounter, refBase);

        allelicCounts.add(new AllelicCount(
                new SimpleInterval(locus.getContig(), locus.getStart(), locus.getEnd()),
                refReadCount, altReadCount, refBase, altBase));
    }

    /**
     * Get the allelic counts gathered so far.
     *
     * @return a *reference* to the AllelicCountCollection
     */
    public AllelicCountCollection getAllelicCounts() {
        return allelicCounts;
    }

    /**
     * Returns the non-ref base with highest count (if there is a tie, the first base in the order given in
     * {@link AllelicCountCollector#BASES} will be returned).
     */
    private static Nucleotide inferAltFromPileupBaseCounts(final Nucleotide.Counter baseCounts,
                                                           final Nucleotide refNucleotide) {
        return BASES.stream()
                .filter(b -> b != refNucleotide)
                .sorted((b1, b2) -> Long.compare(baseCounts.get(b1), baseCounts.get(b2)))
                .findFirst().get();
    }
}
