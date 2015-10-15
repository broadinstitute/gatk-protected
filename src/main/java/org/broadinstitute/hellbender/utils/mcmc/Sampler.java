package org.broadinstitute.hellbender.utils.mcmc;


import org.apache.commons.math3.random.RandomGenerator;

public interface Sampler<U, S extends AbstractParameterizedState, T extends DataCollection> {
    U sample(final RandomGenerator rng, final S state, final T dataCollection);
}