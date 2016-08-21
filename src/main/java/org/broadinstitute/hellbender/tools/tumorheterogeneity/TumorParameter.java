package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;


/**
 * Enumerates the parameters for {@link TumorState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum TumorParameter implements ParameterEnum {
    CONCENTRATION("TH_concentration"),
    VARIANCE("TH_variance"),
    POPULATION_FRACTIONS("TH_population_fractions"),
    MEANS("TH_means");

    public final String name;

    TumorParameter(final String name) {
        this.name = name;
    }
}

