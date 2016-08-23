package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.primitives.Doubles;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.mcmc.ParameterSampler;
import org.broadinstitute.hellbender.utils.mcmc.SliceSampler;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
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

    protected static final class VarianceSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private final double varianceMin;
        private final double varianceMax;
        private final double varianceSliceSamplingWidth;

        public VarianceSampler(final double varianceMin, final double varianceMax, final double varianceSliceSamplingWidth) {
            this.varianceMin = varianceMin;
            this.varianceMax = varianceMax;
            this.varianceSliceSamplingWidth = varianceSliceSamplingWidth;
        }

        public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            final Function<Double, Double> logConditionalPDF = newVariance -> {
                double ll = -0.5 * dataCollection.numPoints() * Math.log(newVariance);
                final double inverseDenominator = 1. / (2. * newVariance);
                for (int dataIndex = 0; dataIndex < dataCollection.numPoints(); dataIndex++) {
                    for (int populationIndex = 0; populationIndex < state.numPopulations(); populationIndex++) {
                        if (state.isInPopulation(dataIndex, populationIndex)) {
                            final double meanDifference = dataCollection.getPoint(dataIndex) - state.mean(populationIndex);
                            ll -= meanDifference * meanDifference * inverseDenominator;
                            break;
                        }
                    }
                }
                return ll;
            };
            return new SliceSampler(rng, logConditionalPDF, varianceMin, varianceMax, varianceSliceSamplingWidth).sample(state.variance());
        }
    }

    protected static final class PopulationFractionsSampler implements ParameterSampler<TumorHeterogeneityState.PopulationFractions, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        public PopulationFractionsSampler() {
        }

        public TumorHeterogeneityState.PopulationFractions sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            final int numPopulations = state.numPopulations();
            final List<MutableInt> populationCounts = IntStream.range(0, numPopulations).boxed().map(j -> new MutableInt(0)).collect(Collectors.toList());
            for (int dataIndex = 0; dataIndex < dataCollection.numPoints(); dataIndex++) {
                for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
                    if (state.isInPopulation(dataIndex, populationIndex)) {
                        populationCounts.get(populationIndex).increment();
                        break;
                    }
                }
            }
            final List<Double> unnormalizedPopulationFractions = populationCounts.stream()
                    .map(c -> new GammaDistribution(rng, state.concentration() + c.doubleValue(), 1.).sample()).collect(Collectors.toList());
            final List<Double> populationFractions = Doubles.asList(MathUtils.normalizeFromRealSpace(Doubles.toArray(unnormalizedPopulationFractions)));
            return new TumorHeterogeneityState.PopulationFractions(populationFractions);
        }
    }

    protected static final class MeansSampler implements ParameterSampler<TumorHeterogeneityState.Means, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private final double meanMin;
        private final double meanMax;
        private final double meanSliceSamplingWidth;

        public MeansSampler(final double meanMin, final double meanMax, final double meanSliceSamplingWidth) {
            this.meanMin = meanMin;
            this.meanMax = meanMax;
            this.meanSliceSamplingWidth = meanSliceSamplingWidth;
        }

        public TumorHeterogeneityState.Means sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            final int numPopulations = state.numPopulations();
            final List<Double> newMeans = new ArrayList<>(numPopulations);
            for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
                final int j = populationIndex;
                final Function<Double, Double> logConditionalPDF = newMean -> {
                    final double inverseDenominator = 1. / (2. * state.variance());
                    return IntStream.range(0, dataCollection.numPoints())
                            .filter(i -> state.isInPopulation(i, j))
                            .mapToDouble(i -> -(dataCollection.getPoint(i) - newMean) * (dataCollection.getPoint(i) - newMean) * inverseDenominator)
                            .sum();
                };
                newMeans.add(new SliceSampler(rng, logConditionalPDF, meanMin, meanMax, meanSliceSamplingWidth).sample(state.mean(populationIndex)));
            }
            return new TumorHeterogeneityState.Means(newMeans);
        }
    }

    protected static final class PopulationIndicatorsSampler implements ParameterSampler<TumorHeterogeneityState.PopulationIndicators, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        public PopulationIndicatorsSampler() {}

        public TumorHeterogeneityState.PopulationIndicators sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            final int numPopulations = state.numPopulations();
            final List<Integer> populationIndicators = new ArrayList<>(numPopulations);
            final List<Integer> populationIndices = IntStream.range(0, numPopulations).boxed().collect(Collectors.toList());
            final int numPoints = dataCollection.numPoints();
            final double inverseDenominator = 1. / (2. * state.variance());
            final double logPrefactor = 0.5 * Math.log(inverseDenominator / Math.PI);
            for (int dataIndex = 0; dataIndex < numPoints; dataIndex++) {
                final double point = dataCollection.getPoint(dataIndex);
                final List<Double> unnormalizedPopulationLog10Probabilities = IntStream.range(0, numPopulations).boxed()
                        .map(j -> Math.log(state.populationFraction(j) + EPSILON) + logPrefactor - (point - state.mean(j)) * (point - state.mean(j)) * inverseDenominator)
                        .map(MathUtils::logToLog10)
                        .collect(Collectors.toList());
                final double[] normalizedPopulationProbabilities = MathUtils.normalizeFromLog10(Doubles.toArray(unnormalizedPopulationLog10Probabilities));
                final Function<Integer, Double> probabilityFunction = j -> normalizedPopulationProbabilities[j];
                final int populationIndex = GATKProtectedMathUtils.randomSelect(populationIndices, probabilityFunction, rng);
                populationIndicators.add(populationIndex);
            }
            return new TumorHeterogeneityState.PopulationIndicators(populationIndicators);
        }
    }
}
