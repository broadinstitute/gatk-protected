package org.broadinstitute.hellbender.utils.mcmc;


import org.apache.commons.math3.random.RandomGenerator;

public interface Sampler<T> {
    T sample(final RandomGenerator rng, final ParameterizedState state, final DataCollection dataCollection);
}