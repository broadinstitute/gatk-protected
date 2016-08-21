package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.exome.copyratio.*;
import org.broadinstitute.hellbender.tools.exome.copyratio.CopyRatioParameter;
import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorState extends ParameterizedState<TumorParameter> {
    public static final class PopulationFractions extends ArrayList<Double> {
        private static final long serialVersionUID = 32905L;
        public PopulationFractions(final List<Double> populationFractions) {
            super(new ArrayList<>(populationFractions));
        }
    }

    public static final class Means extends ArrayList<Double> {
        private static final long serialVersionUID = 23059L;
        public Means(final List<Double> means) {
            super(new ArrayList<>(means));
        }
    }

    public TumorState(final double concentration, final double variance,
                      final PopulationFractions populationFractions, final Means means) {
        super(Arrays.asList(
                new Parameter<>(TumorParameter.CONCENTRATION, concentration),
                new Parameter<>(TumorParameter.VARIANCE, variance),
                new Parameter<>(TumorParameter.POPULATION_FRACTIONS, populationFractions),
                new Parameter<>(TumorParameter.MEANS, means)));
    }

    public double concentration() {
        return get(TumorParameter.CONCENTRATION, Double.class);
    }

    public double variance() {
        return get(TumorParameter.VARIANCE, Double.class);
    }

    public double populationFraction(final int index) {
        return get(TumorParameter.POPULATION_FRACTIONS, PopulationFractions.class).get(index);
    }

    public double mean(final int index) {
        return get(TumorParameter.MEANS, Means.class).get(index);
    }
}
