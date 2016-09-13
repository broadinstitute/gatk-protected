package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;


/**
 * Enumerates the parameters for {@link TumorHeterogeneityState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum TumorHeterogeneityParameter implements ParameterEnum {
    CONCENTRATION("TH_concentration"),
    POPULATION_FRACTIONS("TH_population_fractions"),
    POPULATION_INDICATORS("TH_population_indicators"),
    VARIANT_SEGMENT_FRACTION_HYPERPARAMETERS("TH_variant_segment_fraction_hyperparameters"),
    POPULATION_STATES("TH_population_states");

    public final String name;

    TumorHeterogeneityParameter(final String name) {
        this.name = name;
    }
}

