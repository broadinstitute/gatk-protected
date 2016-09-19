package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;


/**
 * Enumerates the parameters for {@link TumorHeterogeneityState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum TumorHeterogeneityParameter implements ParameterEnum {
    CONCENTRATION("concentration"),
    POPULATION_FRACTIONS("population_fractions"),
    POPULATION_INDICATORS("population_indicators"),
    VARIANT_SEGMENT_FRACTION_HYPERPARAMETERS("variant_segment_fraction_hyperparameters"),
    VARIANT_PROFILES("variant_profiles");

    public final String name;

    TumorHeterogeneityParameter(final String name) {
        this.name = name;
    }
}

