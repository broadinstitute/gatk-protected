package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;
import org.broadinstitute.hellbender.utils.mcmc.Decile;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link DataCollection} for the tumor-heterogeneity model.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityData implements DataCollection {
    private final List<ACNVModeledSegment> segments;
    private final int numBinsCopyRatio;
    private final int numBinsMinorAlleleFraction;
    private final List<ACNVSegmentPosterior> segmentPosteriors;

    public TumorHeterogeneityData(final List<ACNVModeledSegment> segments,
                                  final int numBinsCopyRatio,
                                  final int numBinsMinorAlleleFraction) {
        Utils.nonNull(segments);
        Utils.validateArg(numBinsCopyRatio > 0, "Number of copy-ratio bins must be positive.");
        Utils.validateArg(numBinsMinorAlleleFraction > 0, "Number of minor-allele-fraction bins must be positive.");
        this.segments = new ArrayList<>(segments);
        this.numBinsCopyRatio = numBinsCopyRatio;
        this.numBinsMinorAlleleFraction = numBinsMinorAlleleFraction;
        segmentPosteriors = segments.stream().map(ACNVSegmentPosterior::new).collect(Collectors.toList());
    }

    public double logProbability(final int segmentIndex, final double copyRatio, final double minorAlleleFraction) {
        return segmentPosteriors.get(segmentIndex).logProbability(copyRatio, minorAlleleFraction);
    }
    
    private final class ACNVSegmentPosterior {
        private static final int NUM_BINS = (DecileCollection.NUM_DECILES - 1) * (DecileCollection.NUM_DECILES - 1);
        private static final double BIN_PROBABILITY = 1. / NUM_BINS;
        private static final double EPSILON_BIN_AREA = 1E-10;
        private static final double EPSILON_LOG_POSTERIOR = -100;

        private final List<Double> copyRatioDecileValues;
        private final List<Double> minorAlleleFractionDecileValues;
        private final Map<Pair, Double> binToLogPosteriorMap = new HashMap<>(NUM_BINS); //keys are (lower CR bin edge, lower MAF bin edge)
        
        ACNVSegmentPosterior(final ACNVModeledSegment segment) {
            final List<Double> log2CopyRatioDecileValues = segment.getSegmentMeanPosteriorSummary().getDeciles().getAll();
            copyRatioDecileValues = log2CopyRatioDecileValues.stream().map(log2cr -> Math.pow(2, log2cr)).collect(Collectors.toList());
            minorAlleleFractionDecileValues = segment.getMinorAlleleFractionPosteriorSummary().getDeciles().getAll();
            initializeBinToLogPosteriorMap(binToLogPosteriorMap, copyRatioDecileValues, minorAlleleFractionDecileValues);
        }

        double logProbability(final double copyRatio, final double minorAlleleFraction) {
            final int crIndex = Collections.binarySearch(copyRatioDecileValues, copyRatio);
            //take care of case where no hets separately
        }
        
        private void initializeBinToLogPosteriorMap(final Map<Pair, Double> binToLogPosteriorMap,
                                                    final List<Double> copyRatioDecileValues,
                                                    final List<Double> minorAlleleFractionDecileValues) {
            final List<Double> copyRatioBinSizes = calculateBinSizesFromDecileValues(copyRatioDecileValues);
            final List<Double> minorAlleleFractionBinSizes = calculateBinSizesFromDecileValues(minorAlleleFractionDecileValues);
            //take care of case where no hets separately
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