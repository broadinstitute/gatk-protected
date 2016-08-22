package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityState extends ParameterizedState<TumorHeterogeneityParameter> {
    private final int numPopulations;

    public static final class PopulationFractions extends ArrayList<Double> {
        private static final long serialVersionUID = 32905L;
        public PopulationFractions(final List<Double> populationFractions) {
            super(populationFractions);
        }
    }

    public static final class Means extends ArrayList<Double> {
        private static final long serialVersionUID = 23059L;
        public Means(final List<Double> means) {
            super(means);
        }
    }

    static final class PopulationIndicator extends ArrayList<Boolean> {
        private static final long serialVersionUID = 44561L;
        public PopulationIndicator(final List<Boolean> oneOfNumPopulationsArray) {
            super(oneOfNumPopulationsArray);
        }
    }

    public static final class PopulationIndicators extends ArrayList<ArrayList<Boolean>> {
        private static final long serialVersionUID = 44345L;
        public PopulationIndicators(final List<PopulationIndicator> populationIndicators) {
            super(populationIndicators);
        }
    }

    public TumorHeterogeneityState(final double concentration, final double variance,
                                   final PopulationFractions populationFractions, final Means means,
                                   final PopulationIndicators populationIndicators) {
        super(Arrays.asList(
                new Parameter<>(TumorHeterogeneityParameter.CONCENTRATION, concentration),
                new Parameter<>(TumorHeterogeneityParameter.VARIANCE, variance),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_FRACTIONS, populationFractions),
                new Parameter<>(TumorHeterogeneityParameter.MEANS, means),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_INDICATORS, populationIndicators)));
        Utils.validateArg(populationFractions.size() > 0, "Number of populations must be positive.");
        Utils.validateArg(populationIndicators.size() > 0, "Number of points must be positive.");
        Utils.validateArg(populationFractions.size() == means.size(), "Number of populations must be same for fractions and means.");
        Utils.validateArg(populationFractions.size() == populationIndicators.get(0).size(), "Number of populations must be same for fractions and indicators.");
        numPopulations = populationFractions.size();
    }

    public int numPopulations() {
        return numPopulations;
    }

    public double concentration() {
        return get(TumorHeterogeneityParameter.CONCENTRATION, Double.class);
    }

    public double variance() {
        return get(TumorHeterogeneityParameter.VARIANCE, Double.class);
    }

    public double populationFraction(final int populationIndex) {
        return get(TumorHeterogeneityParameter.POPULATION_FRACTIONS, PopulationFractions.class).get(populationIndex);
    }

    public double mean(final int index) {
        return get(TumorHeterogeneityParameter.MEANS, Means.class).get(index);
    }

    public boolean isInPopulation(final int dataIndex, final int populationIndex) {
        return get(TumorHeterogeneityParameter.POPULATION_INDICATORS, PopulationIndicators.class).get(dataIndex).get(populationIndex);
    }
}
