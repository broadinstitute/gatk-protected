package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;


/**
 * Enumerates the parameters for {@link TumorHeterogeneityState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum TumorHeterogeneityParameter implements ParameterEnum {
    CONCENTRATION("CONCENTRATION"),
    COPY_RATIO_NOISE_CONSTANT("CR_NOISE_CONSTANT"),
    COPY_RATIO_NOISE_FACTOR("CR_NOISE_FACTOR"),
    MINOR_ALLELE_FRACTION_NOISE_FACTOR("MAF_NOISE_FACTOR"),
    INITIAL_PLOIDY("INITIAL_PLOIDY"),
    PLOIDY("PLOIDY"),
    POPULATION_MIXTURE("POPULATION_MIXTURE");

    public final String name;

    TumorHeterogeneityParameter(final String name) {
        this.name = name;
    }
}

