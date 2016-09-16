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

        public ConcentrationSampler(final double concentrationMin, final double concentrationMax, final double concentrationSliceSamplingWidth) {
            this.concentrationMin = concentrationMin;
            this.concentrationMax = concentrationMax;
            this.concentrationSliceSamplingWidth = concentrationSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            final int numPopulations = state.numPopulations();
            final Function<Double, Double> logConditionalPDF = newConcentration -> {
                final double populationFractionsTerm = IntStream.range(0, numPopulations)
                        .mapToDouble(i -> (newConcentration - 1) * Math.log(state.populationFraction(i) + EPSILON)).sum();
                return (dataCollection.concentrationPriorAlpha() - 1.) * Math.log(newConcentration) - dataCollection.concentrationPriorBeta() * newConcentration +
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
        private final List<Integer> populationIndices;

        public PopulationIndicatorsSampler(final int numPopulations) {
            populationIndices = Collections.unmodifiableList(IntStream.range(0, numPopulations).boxed().collect(Collectors.toList()));
        }

        @Override
        public TumorHeterogeneityState.PopulationIndicators sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            //TODO
            final List<Integer> populationIndicators = IntStream.range(0, state.numCells()).boxed()
                    .map(i -> GATKProtectedMathUtils.randomSelect(populationIndices, state::populationFraction, rng))
                    .collect(Collectors.toList());
            return new TumorHeterogeneityState.PopulationIndicators(populationIndicators);
        }
    }

    protected static final class PopulationStateCollectionSampler implements ParameterSampler<TumorHeterogeneityState.PopulationStateCollection, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private final List<PopulationStateSampler> populationStateSamplers;

        public PopulationStateCollectionSampler(final int numPopulations) {
            populationStateSamplers = IntStream.range(0, numPopulations).boxed().map(PopulationStateSampler::new).collect(Collectors.toList());
        }

        public TumorHeterogeneityState.PopulationStateCollection sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            final List<TumorHeterogeneityState.PopulationState> populationStates = populationStateSamplers.stream().map(sampler -> sampler.sample(rng, state, dataCollection)).collect(Collectors.toList());
            return new TumorHeterogeneityState.PopulationStateCollection(populationStates);
        }
    }

    protected static final class PopulationStateSampler implements ParameterSampler<TumorHeterogeneityState.PopulationState, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private int populationIndex;
        private final VariantSegmentFractionSampler variantSegmentFractionSampler;
        private final VariantIndicatorsSampler variantIndicatorsSampler;
        private final VariantPloidyStateIndicatorsSampler variantPloidyStateIndicatorsSampler;

        public PopulationStateSampler(final int populationIndex) {
            this.populationIndex = populationIndex;
            variantSegmentFractionSampler = new VariantSegmentFractionSampler();
            variantIndicatorsSampler = new VariantIndicatorsSampler();
            variantPloidyStateIndicatorsSampler = new VariantPloidyStateIndicatorsSampler();
        }

        public TumorHeterogeneityState.PopulationState sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
            final double variantSegmentFraction = variantSegmentFractionSampler.sample(rng, state, dataCollection);
            final TumorHeterogeneityState.PopulationState.VariantIndicators variantIndicators = variantIndicatorsSampler.sample(rng, state, dataCollection);
            final TumorHeterogeneityState.PopulationState.VariantPloidyStateIndicators variantPloidyStateIndicators = variantPloidyStateIndicatorsSampler.sample(rng, state, dataCollection);
            return new TumorHeterogeneityState.PopulationState(variantSegmentFraction, variantIndicators, variantPloidyStateIndicators);
        }

        private final class VariantSegmentFractionSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
            private VariantSegmentFractionSampler() {}

            @Override
            public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
                final int numSegmentsVariant = (int) IntStream.range(0, dataCollection.numSegments()).filter(i -> state.isVariant(populationIndex, i)).count();
                return new BetaDistribution(rng,
                        dataCollection.variantSegmentFractionPriorAlpha() + numSegmentsVariant,
                        dataCollection.variantSegmentFractionPriorBeta() + dataCollection.numSegments() - numSegmentsVariant).sample();
            }
        }

        private final class VariantIndicatorsSampler implements ParameterSampler<TumorHeterogeneityState.PopulationState.VariantIndicators, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
            private VariantIndicatorsSampler() {}

            @Override
            public TumorHeterogeneityState.PopulationState.VariantIndicators sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
                //TODO
                return new TumorHeterogeneityState.PopulationState.VariantIndicators(new ArrayList<>());
            }
        }

        private final class VariantPloidyStateIndicatorsSampler implements ParameterSampler<TumorHeterogeneityState.PopulationState.VariantPloidyStateIndicators, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
            private VariantPloidyStateIndicatorsSampler() {}

            @Override
            public TumorHeterogeneityState.PopulationState.VariantPloidyStateIndicators sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData dataCollection) {
                //TODO
                return new TumorHeterogeneityState.PopulationState.VariantPloidyStateIndicators(new ArrayList<>());
            }
        }
    }
}
