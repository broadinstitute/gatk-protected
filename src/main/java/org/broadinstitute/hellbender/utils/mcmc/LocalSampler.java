package org.broadinstitute.hellbender.utils.mcmc;


import org.apache.commons.math3.random.RandomGenerator;

import java.util.List;

/**
 * Interface for sampling of a local (i.e., segment-level or site-level) parameter.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface LocalSampler<T extends DataCollection> {
    /**
     * Generate a single sample of a local parameter (i.e., a list of segment-level or site-level samples)
     * from its conditional probability density function, which is a vector-valued function that depends on global
     * parameters, local (i.e., segment-level or site-level) parameters, and data.  The probability density function
     * is assumed to be independent across segments/sites for each local parameter.
     * @param rng               random number generator
     * @param globalParameters  list of global-parameter values
     * @param localParameters   list of lists of local-parameter values at each segment or site
     * @param dataCollection    collection of data sets (e.g., site-level ref and alt counts)
     * @return                  single sample of a local parameter (i.e., a list of segment-level or site-level samples)
     */
    List<Double> sample(final RandomGenerator rng,
                        final List<Double> globalParameters, final List<List<Double>> localParameters,
                        final T dataCollection);
}