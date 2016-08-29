package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link DataCollection} for the tumor-heterogeneity model that allows the calculation of log posterior probabilities
 * for (copy ratio, minor-allele fraction) for each {@link ACNVModeledSegment}.  Given a list of
 * {@link ACNVModeledSegment}, a map from (copy ratio, minor-allele fraction) to discretized log posterior probability
 * is constructed for each segment.  The discretized log posterior probabilities are determined by the copy-ratio and
 * minor-allele-fraction posterior deciles for each segment.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityData implements DataCollection {
    private final List<ACNVSegmentPosterior> segmentPosteriors;

    public TumorHeterogeneityData(final List<ACNVModeledSegment> segments) {
        Utils.nonNull(segments);
        segmentPosteriors = segments.stream().map(ACNVSegmentPosterior::new).collect(Collectors.toList());
    }

    public double logProbability(final int segmentIndex, final double copyRatio, final double minorAlleleFraction) {
        Utils.validateArg(0 <= segmentIndex && segmentIndex < segmentPosteriors.size(), "Segment index is not in valid range.");
        Utils.validateArg(copyRatio >= 0, "Copy ratio must be non-negative.");
        Utils.validateArg(0 <= minorAlleleFraction && minorAlleleFraction <= 0.5, "Minor-allele fraction must be in [0, 0.5].");
        return segmentPosteriors.get(segmentIndex).logProbability(copyRatio, minorAlleleFraction);
    }
    
    private final class ACNVSegmentPosterior {
        private static final int NUM_BINS_1D = DecileCollection.NUM_DECILES - 1;
        private static final int NUM_BINS = NUM_BINS_1D * NUM_BINS_1D;
        private static final double BIN_PROBABILITY = 1. / NUM_BINS;
        private static final double EPSILON_BIN_AREA = 1E-10;
        private static final double EPSILON_LOG_POSTERIOR = -100;

        private final boolean isMinorAlleleFractionNaN;
        private final List<Double> copyRatioDecileValues;
        private final List<Double> minorAlleleFractionDecileValues;
        private final Map<Pair, Double> binToLogPosteriorMap = new HashMap<>(NUM_BINS); //keys are (lower CR bin edge, lower MAF bin edge)
        
        ACNVSegmentPosterior(final ACNVModeledSegment segment) {
            final List<Double> log2CopyRatioDecileValues = segment.getSegmentMeanPosteriorSummary().getDeciles().getAll();
            copyRatioDecileValues = log2CopyRatioDecileValues.stream().map(log2cr -> Math.pow(2, log2cr)).collect(Collectors.toList());
            minorAlleleFractionDecileValues = segment.getMinorAlleleFractionPosteriorSummary().getDeciles().getAll();
            isMinorAlleleFractionNaN = Double.isNaN(segment.getMinorAlleleFractionPosteriorSummary().getCenter());
            initializeBinToLogPosteriorMap(binToLogPosteriorMap, copyRatioDecileValues, minorAlleleFractionDecileValues, isMinorAlleleFractionNaN);
        }

        double logProbability(final double copyRatio, final double minorAlleleFraction) {
            final int crIndex = Collections.binarySearch(copyRatioDecileValues, copyRatio);
            final int mafIndex = Collections.binarySearch(minorAlleleFractionDecileValues, minorAlleleFraction);
            return binToLogPosteriorMap.getOrDefault(Pair.of(crIndex, mafIndex), EPSILON_LOG_POSTERIOR);
        }
        
        private void initializeBinToLogPosteriorMap(final Map<Pair, Double> binToLogPosteriorMap,
                                                    final List<Double> copyRatioDecileValues,
                                                    final List<Double> minorAlleleFractionDecileValues,
                                                    final boolean isMinorAlleleFractionNaN) {
            final List<Double> copyRatioBinSizes = calculateBinSizesFromDecileValues(copyRatioDecileValues);
            final List<Double> minorAlleleFractionBinSizes = isMinorAlleleFractionNaN ?
                    Collections.nCopies(NUM_BINS_1D, 0.5 / NUM_BINS_1D) : calculateBinSizesFromDecileValues(minorAlleleFractionDecileValues);
            for (int crIndex = 0; crIndex < DecileCollection.NUM_DECILES; crIndex++) {
                for (int mafIndex = 0; mafIndex < DecileCollection.NUM_DECILES; mafIndex++) {
                    final double copyRatioDecile = copyRatioDecileValues.get(crIndex);
                    final double minorAlleleFractionDecile = minorAlleleFractionDecileValues.get(mafIndex);
                    final Pair bin = Pair.of(copyRatioDecile, minorAlleleFractionDecile);
                    final double binArea = copyRatioBinSizes.get(crIndex) * minorAlleleFractionBinSizes.get(mafIndex);
                    final double binProbabilityDensity = BIN_PROBABILITY / (binArea + EPSILON_BIN_AREA);
                    final double logPosterior = Math.log(binProbabilityDensity);
                    binToLogPosteriorMap.put(bin, logPosterior);
                }
            }
        }
        
        private List<Double> calculateBinSizesFromDecileValues(final List<Double> decileValues) {
            return IntStream.range(1, DecileCollection.NUM_DECILES).boxed()
                    .map(i -> decileValues.get(i) - decileValues.get(i - 1))
                    .collect(Collectors.toList());
        }
    }
}