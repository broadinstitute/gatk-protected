package org.broadinstitute.hellbender.utils.mcmc;


import org.apache.commons.math3.random.RandomGenerator;

import java.util.List;

/**
 * Interface for sampling of a global parameter.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface GlobalSampler {
    /**
     * Generate a single sample of a global parameter from its conditional probability density function,
     * which is a univariate function that depends on global parameters, local (i.e., segment-level or site-level)
     * parameters, and data.  The probability density function is assumed to be independent across segments/sites
     * for each local parameter.
     * @param rng               random number generator
     * @param globalParameters  list of global-parameter values
     * @param localParameters   list of lists of local-parameter values at each segment or site
     * @param data              list of lists of data (e.g., a list of lists of site-level ref and alt counts)
     * @return                  single sample of a global parameter
     */
    double sample(final RandomGenerator rng,
                  final List<Double> globalParameters, final List<List<Double>> localParameters,
                  final List<List<Double>> data);
}