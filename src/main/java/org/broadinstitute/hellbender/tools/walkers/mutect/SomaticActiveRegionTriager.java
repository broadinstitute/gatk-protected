package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.List;

/**
 * Created by David Benjamin on 5/23/17.
 */
public class SomaticActiveRegionTriager {
    private static final AltEvidence tumorEvidence = new AltEvidence();
    private static final AltEvidence normalEvidence = new AltEvidence();
    private final double tumorLogOddsThreshold;
    private final double normalLogOddsThreshold;

    public SomaticActiveRegionTriager(final M2ArgumentCollection MTAC) {
        tumorLogOddsThreshold = MathUtils.log10ToLog(MTAC.initialTumorLog10OddsThreshold);
        normalLogOddsThreshold = MathUtils.log10ToLog(MTAC.initialNormalLog10OddsThreshold);
    }

    public boolean isActive(final ReadPileup tumorPileup, final byte refBase) {
        tumorEvidence.findQuals(tumorPileup, refBase);
        return AltEvidence.isEvidenceOfSomaticVariant(tumorEvidence, tumorLogOddsThreshold, tumorPileup.size());
    }

    public boolean isActive(final ReadPileup tumorPileup, final ReadPileup normalPileup, final byte refBase) {
        tumorEvidence.findQuals(tumorPileup, refBase);
        normalEvidence.findQuals(normalPileup, refBase);
        return AltEvidence.isEvidenceOfSomaticVariant(tumorEvidence, normalEvidence, tumorLogOddsThreshold, normalLogOddsThreshold, tumorPileup.size(), normalPileup.size());
    }


    private static class AltEvidence {
        private static final byte DEFAULT_INSERTION_QUAL = (byte) 40;
        private static final byte DEFAULT_SOFT_CLIP_QUAL = (byte) 30;
        private static final byte DEFAULT_DELETION_QUAL = (byte) 40;
        private static final byte MIN_BASE_QUAL = (byte) 15;
        private static final List<Nucleotide> regularBases = Arrays.asList(Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.T);

        private final EnumMap<Nucleotide, List<Byte>> snvQuals = new EnumMap<>(Nucleotide.class);
        private final List<Byte> deletionQuals = new ArrayList<>();
        private final List<Byte> insertionQuals = new ArrayList<>();

        public AltEvidence() {
            Arrays.stream(Nucleotide.values()).forEach(nucleotide -> snvQuals.put(nucleotide, new ArrayList<>()));
        }

        public void findQuals(final ReadPileup pileup, final byte refBase) {
            clear();

            for (final PileupElement el : pileup) {
                if (el.isDeletion() || el.isBeforeDeletionStart() || el.isAfterDeletionEnd()) {
                    deletionQuals.add(DEFAULT_DELETION_QUAL);
                } else if (el.isNextToSoftClip()) {
                    deletionQuals.add(DEFAULT_SOFT_CLIP_QUAL);
                }

                if (el.isBeforeInsertion() || el.isAfterInsertion()) {
                    insertionQuals.add(DEFAULT_INSERTION_QUAL);
                }

                final byte base = el.getBase();
                if (base != refBase && el.getQual() > MIN_BASE_QUAL) {
                    snvQuals.get(Nucleotide.valueOf(base)).add(el.getQual());
                }
            }

        }

        public static boolean isEvidenceOfSomaticVariant(final AltEvidence tumorEvidence, final double logOddsThreshold, final int totalDepth) {
            for (final Nucleotide nucleotide : regularBases) {
                if (logOdds(tumorEvidence.snvQuals.get(nucleotide), totalDepth) > logOddsThreshold) {
                    return true;
                }
            }

            return logOdds(tumorEvidence.deletionQuals, totalDepth) > logOddsThreshold || logOdds(tumorEvidence.insertionQuals, totalDepth) > logOddsThreshold;
        }

        public static boolean isEvidenceOfSomaticVariant(final AltEvidence tumorEvidence, final AltEvidence normalEvidence,
                                                         final double tumorLogOddsThreshold, final double normalLogOddsThreshold,
                                                         final int tumorDepth, final int normalDepth) {
            for (final Nucleotide nucleotide : regularBases) {
                if (logOdds(tumorEvidence.snvQuals.get(nucleotide), tumorDepth) > tumorLogOddsThreshold && logOdds(normalEvidence.snvQuals.get(nucleotide), normalDepth) < normalLogOddsThreshold) {
                    return true;
                }
            }

            return (logOdds(tumorEvidence.deletionQuals, tumorDepth) > tumorLogOddsThreshold && logOdds(normalEvidence.deletionQuals, normalDepth) < normalLogOddsThreshold)
                    || (logOdds(tumorEvidence.insertionQuals, tumorDepth) > tumorLogOddsThreshold && logOdds(normalEvidence.insertionQuals, normalDepth) < normalLogOddsThreshold);
        }

        private void clear() {
            Arrays.stream(Nucleotide.values()).forEach(nucleotide -> snvQuals.get(nucleotide).clear());
            deletionQuals.clear();
            insertionQuals.clear();
        }

        private static double entropy(final double f) {
            return f < 0.0001 || f > 0.999 ? 0 : f * FastMath.log(f) + (1 - f) * FastMath.log(1 - f);
        }

        private static double logOdds(final List<Byte> quals, final int totalDepth) {
            final double f = (double) quals.size() / totalDepth;
            return totalDepth * entropy(f) - MathUtils.log10ToLog(quals.stream().mapToDouble(QualityUtils::qualToErrorProbLog10).sum());

        }
    }
}
