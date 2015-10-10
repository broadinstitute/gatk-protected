package org.broadinstitute.hellbender.utils.mcmc;


import org.apache.commons.math3.random.RandomGenerator;

public interface Sampler<S, T extends ParameterizedState, U extends DataCollection> {
    S sample(final RandomGenerator rng, final T modelState, final U dataCollection);
}