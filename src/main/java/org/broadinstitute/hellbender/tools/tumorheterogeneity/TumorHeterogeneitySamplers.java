package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.primitives.Doubles;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.mcmc.ParameterSampler;
import org.broadinstitute.hellbender.utils.mcmc.SliceSampler;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class TumorHeterogeneitySamplers {
    private static final double EPSILON = 1E-10;

    private TumorHeterogeneitySamplers() {}

    protected static final class ConcentrationSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private final double concentrationMin;
        private final double concentrationMax;
        private final double concentrationSliceSamplingWidth;
        private final double concentrationPriorAlpha;
        private final double concentrationPriorBeta;

        public ConcentrationSampler(final double concentrationMin, final double concentrationMax, final double concentrationSliceSamplingWidth,
                                    final double concentrationPriorAlpha, final double concentrationPriorBeta) {
            this.concentrationMin = concentrationMin;
            this.concentrationMax = concentrationMax;
            this.concentrationSliceSamplingWidth = concentrationSliceSamplingWidth;
            this.concentrationPriorAlpha = concentrationPriorAlpha;
            this.concentrationPriorBeta = concentrationPriorBeta;
        }

        @Override
        public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            final int numPopulations = state.numPopulations();
            final Function<Double, Double> logConditionalPDF = newConcentration -> {
                final double populationFractionsTerm = IntStream.range(0, numPopulations)
                        .mapToDouble(i -> (newConcentration - 1) * Math.log(state.populationFraction(i) + EPSILON)).sum();
                return (concentrationPriorAlpha - 1.) * Math.log(newConcentration) - concentrationPriorBeta * newConcentration +
                        Gamma.logGamma(newConcentration * numPopulations) - numPopulations * Gamma.logGamma(newConcentration) + populationFractionsTerm;
            };
            return new SliceSampler(rng, logConditionalPDF, concentrationMin, concentrationMax, concentrationSliceSamplingWidth).sample(state.concentration());
        }
    }

    protected static final class PopulationFractionsSampler implements ParameterSampler<TumorHeterogeneityState.PopulationFractions, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        public PopulationFractionsSampler() {
        }

        @Override
        public TumorHeterogeneityState.PopulationFractions sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            final List<MutableInt> populationCounts = sumPopulationCounts(state);
            final List<Double> populationFractions = samplePopulationFractions(rng, state, populationCounts);
            return new TumorHeterogeneityState.PopulationFractions(populationFractions);
        }

        private List<MutableInt> sumPopulationCounts(final TumorHeterogeneityState state) {
            final int numPopulations = state.numPopulations();
            final List<MutableInt> populationCounts = IntStream.range(0, numPopulations).boxed()
                    .map(j -> new MutableInt(0)).collect(Collectors.toList());
            for (int cellIndex = 0; cellIndex < state.numCells(); cellIndex++) {
                final int populationIndex = state.populationIndex(cellIndex);
                populationCounts.get(populationIndex).increment();
            }
            return populationCounts;
        }

        private List<Double> samplePopulationFractions(final RandomGenerator rng, final TumorHeterogeneityState state, final List<MutableInt> populationCounts) {
            //sampling from Dirichlet(alpha_vec) is equivalent to sampling from individual Gamma(alpha_vec_i, 1) distributions and normalizing
            final double[] unnormalizedPopulationFractions = populationCounts.stream()
                    .mapToDouble(c -> new GammaDistribution(rng, state.concentration() + c.doubleValue(), 1.).sample()).toArray();
            return Doubles.asList(MathUtils.normalizeFromRealSpace(unnormalizedPopulationFractions));
        }
    }

    protected static final class PopulationIndicatorsSampler implements ParameterSampler<TumorHeterogeneityState.PopulationIndicators, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        public PopulationIndicatorsSampler() {}

        @Override
        public TumorHeterogeneityState.PopulationIndicators sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            final int numPopulations = state.numPopulations();
            final List<Integer> populationIndicators = new ArrayList<>(numPopulations);
            final List<Integer> populationIndices = IntStream.range(0, numPopulations).boxed().collect(Collectors.toList());
            final int numPoints = dataCollection.numPoints();
            final double inverseDenominator = 1. / (2. * state.variance());
            final double log10Prefactor = 0.5 * Math.log10(inverseDenominator / Math.PI);
            for (int dataIndex = 0; dataIndex < numPoints; dataIndex++) {
                final double[] unnormalizedPopulationLog10Probabilities = new double[numPopulations];
                for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
                    final double meanDifference = dataCollection.getPoint(dataIndex) - state.mean(populationIndex);
                    unnormalizedPopulationLog10Probabilities[populationIndex] =
                            Math.log10(state.populationFraction(populationIndex) + EPSILON) + log10Prefactor
                                    - meanDifference * meanDifference * inverseDenominator * MathUtils.LOG10_OF_E;
                }
                final double[] normalizedPopulationProbabilities = MathUtils.normalizeFromLog10(unnormalizedPopulationLog10Probabilities);
                final Function<Integer, Double> probabilityFunction = j -> normalizedPopulationProbabilities[j];
                final int populationIndex = GATKProtectedMathUtils.randomSelect(populationIndices, probabilityFunction, rng);
                populationIndicators.add(populationIndex);
            }
            return new TumorHeterogeneityState.PopulationIndicators(populationIndicators);
        }
    }

    protected static final class VariantSegmentFractionHyperparametersSampler implements ParameterSampler<TumorHeterogeneityState.HyperparameterValues, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        public VariantSegmentFractionHyperparametersSampler() {}

        @Override
        public TumorHeterogeneityState.HyperparameterValues sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            return new TumorHeterogeneityState.HyperparameterValues(1., 1.);
        }
    }

    protected static final class PopulationCollectionSampler implements ParameterSampler<TumorHeterogeneityState.PopulationCollection, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private int numPopulations;
        private final List<PopulationStateSampler> populationStateSamplers;

        public PopulationCollectionSampler(final int numPopulations) {
            this.numPopulations = numPopulations;
            populationStateSamplers = IntStream.range(0, numPopulations).boxed().map(PopulationStateSampler::new).collect(Collectors.toList());
        }

        public TumorHeterogeneityState.PopulationCollection sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
        }
    }


    protected static final class PopulationStateSampler implements ParameterSampler<TumorHeterogeneityState.PopulationState, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private int populationIndex;
        private final VariantSegmentFractionSampler variantSegmentFractionSampler;

        public PopulationStateSampler(final int populationIndex) {
            this.populationIndex = populationIndex;
            variantSegmentFractionSampler = new VariantSegmentFractionSampler(populationIndex);
        }

        public TumorHeterogeneityState.PopulationState sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            final double variantSegmentFraction = variantSegmentFractionSampler.sample(rng, state, dataCollection);
            return new TumorHeterogeneityState.PopulationState(variantSegmentFraction, );
        }

        private static final class VariantSegmentFractionSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
            private final int populationIndex;

            public VariantSegmentFractionSampler(final int populationIndex) {
                this.populationIndex = populationIndex;
            }

            @Override
            public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
                final int numSegmentsVariant = (int) IntStream.range(0, dataCollection.numSegments()).filter(i -> state.isVariant(populationIndex, i)).count();
                return new BetaDistribution(rng,
                        state.variantSegmentFractionPriorAlpha() + numSegmentsVariant,
                        state.variantSegmentFractionPriorBeta() + dataCollection.numSegments() - numSegmentsVariant).sample();
            }
        }
    }
}
