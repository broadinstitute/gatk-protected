package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.primitives.Doubles;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
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
        private static final double INITIAL_HIDDEN_POINT_STEP_SIZE = 0.01;
        private static final double DEFAULT_OPTIMAL_ACCEPTANCE_RATE = 0.4;
        private static final double DEFAULT_TIME_SCALE = 100;
        private static final double DEFAULT_ADJUSTMENT_RATE = 1.;

        private final List<Double> hiddenPoints;
        private double hiddenPointStepSize;
        private int iteration;

        public PopulationFractionsSampler(final RandomGenerator rng, final int numPopulations) {
            hiddenPoints = IntStream.range(0, numPopulations).boxed().map(i -> new ExponentialDistribution(rng, 1.).sample()).collect(Collectors.toList());
            hiddenPointStepSize = INITIAL_HIDDEN_POINT_STEP_SIZE;
            iteration = 0;
        }

        public TumorHeterogeneityState.PopulationFractions sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            final List<MutableInt> populationCounts = sumPopulationCounts(state, dataCollection);
            final List<Double> populationFractions = samplePopulationFractions(rng, state, populationCounts);
            return new TumorHeterogeneityState.PopulationFractions(populationFractions);
        }

        private List<MutableInt> sumPopulationCounts(final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
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
            return populationCounts;
        }

        private List<Double> samplePopulationFractions(final RandomGenerator rng, final TumorHeterogeneityState state, final List<MutableInt> populationCounts) {
            final double[] unnormalizedPopulationFractions = populationCounts.stream()
                    .mapToDouble(c -> new GammaDistribution(rng, state.concentration() + c.doubleValue(), 1.).sample()).toArray();
            return Doubles.asList(MathUtils.normalizeFromRealSpace(unnormalizedPopulationFractions));
        }

        private List<Double> samplePopulationFractionsMH(final RandomGenerator rng, final TumorHeterogeneityState state, final List<MutableInt> populationCounts) {
            final int numPopulations = state.numPopulations();
            final List<Double> currentPopulationFractions = IntStream.range(0, numPopulations).boxed().map(state::populationFraction).collect(Collectors.toList());
            final MutableDouble sum = new MutableDouble(0.);
            for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
                final double currentPoint = hiddenPoints.get(populationIndex);
                final double nextPoint = sampleHiddenPoint(rng, currentPoint, iteration);
                hiddenPoints.set(populationIndex, nextPoint);
                sum.add(nextPoint);
            }
            final List<Double> nextPopulationFractions = hiddenPoints.stream().map(p -> p / sum.doubleValue()).collect(Collectors.toList());
            final double metropolisLogRatio = IntStream.range(0, numPopulations).boxed()
                    .mapToDouble(j -> (state.concentration() - 1. + populationCounts.get(j).doubleValue()) *
                            (Math.log(nextPopulationFractions.get(j) + EPSILON) - Math.log(currentPopulationFractions.get(j) + EPSILON)))
                    .sum();
            final double acceptanceProbability = Math.min(1., Math.exp(metropolisLogRatio));
            iteration++;
            return rng.nextDouble() < acceptanceProbability ? nextPopulationFractions : new ArrayList<>(currentPopulationFractions);
        }

        private double sampleHiddenPoint(final RandomGenerator rng, final double currentPoint, final int iteration) {
            final double normalSample = new NormalDistribution(rng, 0., hiddenPointStepSize).sample();
            final double nextPoint = currentPoint * Math.exp(normalSample);
            final double metropolisLogRatio = -nextPoint + currentPoint;
            final double hastingsLogRatio = Math.log(nextPoint + EPSILON) - Math.log(currentPoint + EPSILON);
            final double acceptanceProbability = Math.min(1., Math.exp(metropolisLogRatio + hastingsLogRatio));
            final double correctionFactor = (acceptanceProbability - DEFAULT_OPTIMAL_ACCEPTANCE_RATE) * DEFAULT_ADJUSTMENT_RATE *
                    (DEFAULT_TIME_SCALE / (DEFAULT_TIME_SCALE + iteration));
            hiddenPointStepSize *= Math.exp(correctionFactor);
            return rng.nextDouble() < acceptanceProbability ? nextPoint : currentPoint;
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
}
