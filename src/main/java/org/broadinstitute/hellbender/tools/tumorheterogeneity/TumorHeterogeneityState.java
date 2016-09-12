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
                                   final HyperparameterValues variantFractionHyperparameters,
                                   final PopulationStates populationStates) {
        super(Arrays.asList(
                new Parameter<>(TumorHeterogeneityParameter.CONCENTRATION, concentration),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_FRACTIONS, populationFractions),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_INDICATORS, populationIndicators),
                new Parameter<>(TumorHeterogeneityParameter.VARIANT_FRACTION_HYPERPARAMETERS, variantFractionHyperparameters),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_STATES, populationStates)));
        Utils.validateArg(populationFractions.size() > 0, "Number of populations must be positive.");
        Utils.validateArg(populationIndicators.size() > 0, "Number of cells must be positive.");
        Utils.validateArg(populationFractions.size() == populationStates.size(), "Number of populations must be same for fractions and states.");
        Utils.validateArg(Collections.max(populationIndicators) <= populationFractions.size(), "Number of populations must be same for fractions and indicators.");
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
        return get(TumorHeterogeneityParameter.VARIANT_FRACTION_HYPERPARAMETERS, HyperparameterValues.class).alpha;
    }

    public double variantFractionBeta() {
        return get(TumorHeterogeneityParameter.VARIANT_FRACTION_HYPERPARAMETERS, HyperparameterValues.class).beta;
    }

    public double variantFraction(final int populationIndex) {
        return get(TumorHeterogeneityParameter.POPULATION_STATES, PopulationStates.class).get(populationIndex).variantFraction;
    }

    public boolean isVariant(final int populationIndex, final int segmentIndex) {
        return get(TumorHeterogeneityParameter.POPULATION_STATES, PopulationStates.class).get(populationIndex).variantIndicators.get(segmentIndex);
    }

    public boolean isVariantPloidy(final int populationIndex, final int segmentIndex, final int variantPloidyIndex) {
        return get(TumorHeterogeneityParameter.POPULATION_STATES, PopulationStates.class).get(populationIndex).variantPloidyIndicators.get(segmentIndex).equals(variantPloidyIndex);
    }

    static class HyperparameterValues {
        private final double alpha;
        private final double beta;

        HyperparameterValues(final double alpha, final double beta) {
            this.alpha = alpha;
            this.beta = beta;
        }
    }

    /**
     * Represents variant fraction and per-segment variant and ploidy indicators for each population.
     */
    static class PopulationState {
        public static final class VariantIndicators extends ArrayList<Boolean> {
            //list of booleans, size = number of segments, i-th element = true if segment i is variant
            private static final long serialVersionUID = 35746L;
            public VariantIndicators(final List<Boolean> variantIndicators) {
                super(variantIndicators);
            }
        }

        public static final class VariantPloidyIndicators extends ArrayList<Integer> {
            //list of integers, size = number of segments, i-th element = variant-ploidy index of segment i
            private static final long serialVersionUID = 78476L;
            public VariantPloidyIndicators(final List<Integer> variantPloidyIndicators) {
                super(variantPloidyIndicators);
            }
        }

        private final double variantFraction;
        private final VariantIndicators variantIndicators;
        private final VariantPloidyIndicators variantPloidyIndicators;

        PopulationState(final double variantFraction,
                        final VariantIndicators variantIndicators,
                        final VariantPloidyIndicators variantPloidyIndicators) {
            Utils.validateArg(0. <= variantFraction && variantFraction <= 1., "Variant fraction must be in [0, 1].");
            Utils.validateArg(variantIndicators.size() == variantPloidyIndicators.size(), "Number of segments must be same for variant and ploidy indicators.");
            this.variantFraction = variantFraction;
            this.variantIndicators = variantIndicators;
            this.variantPloidyIndicators = variantPloidyIndicators;
        }
    }
}
