package org.broadinstitute.hellbender.tools.exome.allelefractionrevision;

import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.hellbender.utils.mcmc.ParameterSampler;
import org.broadinstitute.hellbender.utils.mcmc.SliceSampler;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Sampler classes for the allele-fraction model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class AlleleFractionSamplers {
    private AlleleFractionSamplers() {}

    // sample sample depth
    protected static final class SampleDepthSampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> {
        private static final double MIN_SAMPLE_DEPTH = 0.;

        private final double maxSampleDepth;
        private final double sampleDepthSliceSamplingWidth;

        public SampleDepthSampler(final double maxSampleDepth, final double sampleDepthSliceSamplingWidth) {
            this.maxSampleDepth = maxSampleDepth;
            this.sampleDepthSliceSamplingWidth = sampleDepthSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            return new SliceSampler(rng, x -> AlleleFractionLikelihoods.logLikelihood(
                    state.globalParameters().copyWithNewSampleDepth(x),  state.minorFractions(), data),
                    MIN_SAMPLE_DEPTH, maxSampleDepth, sampleDepthSliceSamplingWidth)
                    .sample(state.sampleDepth());
        }
    }

    // sample sample bias
    protected static final class SampleBiasSampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> {
        private static final double MIN_SAMPLE_BIAS = 1E-10;

        private final double maxSampleBias;
        private final double sampleBiasSliceSamplingWidth;

        public SampleBiasSampler(final double maxSampleBias, final double sampleBiasSliceSamplingWidth) {
            this.maxSampleBias = maxSampleBias;
            this.sampleBiasSliceSamplingWidth = sampleBiasSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            return new SliceSampler(rng, x -> AlleleFractionLikelihoods.logLikelihood(
                    state.globalParameters().copyWithNewSampleBias(x), state.minorFractions(), data),
                    MIN_SAMPLE_BIAS, maxSampleBias, sampleBiasSliceSamplingWidth)
                    .sample(state.sampleBias());
        }
    }

    // sample outlier probability
    protected static final class OutlierProbabilitySampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> {
        private static final double MIN_OUTLIER_PROBABILITY = 0.;
        private static final double MAX_OUTLIER_PROBABILITY = 1.;

        private final double outlierProbabilitySliceSamplingWidth;

        public OutlierProbabilitySampler(final double outlierProbabilitySliceSamplingWidth) {
            this.outlierProbabilitySliceSamplingWidth = outlierProbabilitySliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            return new SliceSampler(rng, x -> AlleleFractionLikelihoods.logLikelihood(
                    state.globalParameters().copyWithNewOutlierProbability(x), state.minorFractions(), data),
                    MIN_OUTLIER_PROBABILITY, MAX_OUTLIER_PROBABILITY, outlierProbabilitySliceSamplingWidth)
                    .sample(state.outlierProbability());
        }
    }

    // sample minor fraction of a single segment
    private static final class PerSegmentMinorFractionSampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> {
        private static double MIN_MINOR_FRACTION = 0.;
        private static double MAX_MINOR_FRACTION = 0.5;

        private final int segmentIndex;
        private final double sliceSamplingWidth;

        public PerSegmentMinorFractionSampler(final int segmentIndex, final double sliceSamplingWidth) {
            this.segmentIndex = segmentIndex;
            this.sliceSamplingWidth = sliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            if (data.getNumHetsInSegment(segmentIndex) == 0) {
                return Double.NaN;
            }
            return new SliceSampler(rng, f -> AlleleFractionLikelihoods.segmentLogLikelihood(
                    state.globalParameters(), f, data.getCountsInSegment(segmentIndex), data.getPoN()),
                    MIN_MINOR_FRACTION, MAX_MINOR_FRACTION, sliceSamplingWidth)
                    .sample(state.segmentMinorFraction(segmentIndex));
        }
    }

    // sample minor fractions of all segments
    protected static final class MinorFractionsSampler implements ParameterSampler<AlleleFractionState.MinorFractions, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> {
        private final List<PerSegmentMinorFractionSampler> perSegmentSamplers = new ArrayList<>();

        public MinorFractionsSampler(final List<Double> sliceSamplingWidths) {
            final int numSegments = sliceSamplingWidths.size();
            for (int segment = 0; segment < numSegments; segment++) {
                perSegmentSamplers.add(new PerSegmentMinorFractionSampler(segment, sliceSamplingWidths.get(segment)));
            }
        }

        @Override
        public AlleleFractionState.MinorFractions sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            return new AlleleFractionState.MinorFractions(perSegmentSamplers.stream()
                    .map(sampler -> sampler.sample(rng, state, data)).collect(Collectors.toList()));
        }
    }
}
