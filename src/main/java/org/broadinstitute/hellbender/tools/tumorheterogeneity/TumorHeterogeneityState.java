package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityState extends ParameterizedState<TumorHeterogeneityParameter> {
    private static final double POPULATION_FRACTION_NORMALIZATION_EPSILON = 1E-3;

    private final int numPopulations;

    public static final class PopulationFractions extends ArrayList<Double> {
        //list of doubles, size = number of populations, i-th element = fraction of population i by cell number
        private static final long serialVersionUID = 79454L;
        public PopulationFractions(final List<Double> populationFractions) {
            super(populationFractions);
        }
    }

    public static final class PopulationStates extends ArrayList<PopulationState> {
        //list of PopulationStates, size = number of populations, i-th element = PopulationState for population i
        private static final long serialVersionUID = 76498L;
        public PopulationStates(final List<PopulationState> populationStates) {
            super(populationStates);
        }
    }

    public static final class PopulationIndicators extends ArrayList<Integer> {
        //list of integers, size = number of cells, i-th element = population index for cell i
        private static final long serialVersionUID = 81915L;
        public PopulationIndicators(final List<Integer> populationIndicators) {
            super(populationIndicators);
        }
    }

    public TumorHeterogeneityState(final double concentration,
                                   final PopulationFractions populationFractions,
                                   final PopulationIndicators populationIndicators,
                                   final HyperparameterValues variantSegmentFractionHyperparameters,
                                   final PopulationStates populationStates) {
        super(Arrays.asList(
                new Parameter<>(TumorHeterogeneityParameter.CONCENTRATION, concentration),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_FRACTIONS, populationFractions),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_INDICATORS, populationIndicators),
                new Parameter<>(TumorHeterogeneityParameter.VARIANT_SEGMENT_FRACTION_HYPERPARAMETERS, variantSegmentFractionHyperparameters),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_STATES, populationStates)));
        Utils.validateArg(populationFractions.size() > 0,
                "Number of populations must be positive.");
        final double populationFractionNormalization = populationFractions.stream().mapToDouble(Double::doubleValue).sum();
        Utils.validateArg(Math.abs(1. - populationFractionNormalization) <= POPULATION_FRACTION_NORMALIZATION_EPSILON,
                "Population fractions must sum to unity.");
        Utils.validateArg(populationIndicators.size() > 0,
                "Number of cells must be positive.");
        Utils.validateArg(populationFractions.size() == populationStates.size(),
                "Number of populations must be same for fractions and states.");
        Utils.validateArg(Collections.max(populationIndicators) <= populationFractions.size(),
                "Number of populations must be same for population fractions and indicators.");
        numPopulations = populationStates.size();
    }

    public int numPopulations() {
        return numPopulations;
    }

    public double concentration() {
        return get(TumorHeterogeneityParameter.CONCENTRATION, Double.class);
    }

    public double populationFraction(final int populationIndex) {
        return get(TumorHeterogeneityParameter.POPULATION_FRACTIONS, PopulationFractions.class).get(populationIndex);
    }

    public boolean isInPopulation(final int populationIndex, final int cellIndex) {
        return get(TumorHeterogeneityParameter.POPULATION_INDICATORS, PopulationIndicators.class).get(cellIndex).equals(populationIndex);
    }

    public double variantFractionAlpha() {
        return get(TumorHeterogeneityParameter.VARIANT_SEGMENT_FRACTION_HYPERPARAMETERS, HyperparameterValues.class).alpha;
    }

    public double variantFractionBeta() {
        return get(TumorHeterogeneityParameter.VARIANT_SEGMENT_FRACTION_HYPERPARAMETERS, HyperparameterValues.class).beta;
    }

    public double variantSegmentFraction(final int populationIndex) {
        return get(TumorHeterogeneityParameter.POPULATION_STATES, PopulationStates.class).get(populationIndex).variantFraction;
    }

    public boolean isVariant(final int populationIndex, final int segmentIndex) {
        return get(TumorHeterogeneityParameter.POPULATION_STATES, PopulationStates.class).get(populationIndex).variantIndicators.get(segmentIndex);
    }

    public boolean isInVariantPloidyState(final int populationIndex, final int segmentIndex, final int variantPloidyStateIndex) {
        return get(TumorHeterogeneityParameter.POPULATION_STATES, PopulationStates.class).get(populationIndex).variantPloidyStateIndicators.get(segmentIndex).equals(variantPloidyStateIndex);
    }

    static class HyperparameterValues {
        private final double alpha;
        private final double beta;

        HyperparameterValues(final double alpha, final double beta) {
            Utils.validateArg(alpha > 0, "Hyperparameter alpha must be positive.");
            Utils.validateArg(beta > 0, "Hyperparameter beta must be positive.");
            this.alpha = alpha;
            this.beta = beta;
        }
    }

    /**
     * For each population, represents variant-segment fraction and per-segment variant and ploidy indicators.
     */
    static class PopulationState {
        public static final class VariantIndicators extends ArrayList<Boolean> {
            //list of booleans, size = number of segments, i-th element = true if segment i is variant
            private static final long serialVersionUID = 35746L;
            public VariantIndicators(final List<Boolean> variantIndicators) {
                super(variantIndicators);
            }
        }

        public static final class VariantPloidyStateIndicators extends ArrayList<Integer> {
            //list of integers, size = number of segments, i-th element = variant-ploidy-state index of segment i
            private static final long serialVersionUID = 78476L;
            public VariantPloidyStateIndicators(final List<Integer> variantPloidyStateIndicators) {
                super(variantPloidyStateIndicators);
            }
        }

        private final double variantFraction;   //fraction of segments that are variant in this population
        private final VariantIndicators variantIndicators;
        private final VariantPloidyStateIndicators variantPloidyStateIndicators;

        PopulationState(final double variantSegmentFraction,
                        final VariantIndicators variantIndicators,
                        final VariantPloidyStateIndicators variantPloidyStateIndicators) {
            Utils.validateArg(0. <= variantSegmentFraction && variantSegmentFraction <= 1.,
                    "Variant-segment fraction must be in [0, 1].");
            Utils.validateArg(variantIndicators.size() == variantPloidyStateIndicators.size(),
                    "Number of segments must be same for variant and ploidy-state indicators.");
            this.variantFraction = variantSegmentFraction;
            this.variantIndicators = variantIndicators;
            this.variantPloidyStateIndicators = variantPloidyStateIndicators;
        }
    }
}
