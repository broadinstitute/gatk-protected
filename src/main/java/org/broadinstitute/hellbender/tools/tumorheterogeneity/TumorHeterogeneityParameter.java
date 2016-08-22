package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;


/**
 * Enumerates the parameters for {@link TumorHeterogeneityState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum TumorHeterogeneityParameter implements ParameterEnum {
    CONCENTRATION("TH_concentration"),
    VARIANCE("TH_variance"),
    POPULATION_FRACTIONS("TH_population_fractions"),
    MEANS("TH_means"),
    POPULATION_INDICATORS("TH_population_indicators");

    public final String name;

    TumorHeterogeneityParameter(final String name) {
        this.name = name;
    }
}

