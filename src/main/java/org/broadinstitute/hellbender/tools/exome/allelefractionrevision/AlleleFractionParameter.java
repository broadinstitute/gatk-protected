package org.broadinstitute.hellbender.tools.exome.allelefractionrevision;

import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;

/**
 * Enumerates the parameters for {@link AlleleFractionState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum AlleleFractionParameter implements ParameterEnum {
    SAMPLE_DEPTH("AF_sample_depth"),
    SAMPLE_BIAS("AF_sample_bias"),
    OUTLIER_PROBABILITY("AF_outlier_probability"),
    MINOR_ALLELE_FRACTIONS("AF_minor_allele_fractions");

    public final String name;

    AlleleFractionParameter(final String name) {
        this.name = name;
    }
}
