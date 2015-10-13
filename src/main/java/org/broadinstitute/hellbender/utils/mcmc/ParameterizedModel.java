package org.broadinstitute.hellbender.utils.mcmc;

import org.apache.commons.math3.random.RandomGenerator;

import java.util.HashMap;
import java.util.Map;

public final class ParameterizedModel {
    private enum UpdateMethod {
        GIBBS
    }

    private final ParameterizedState state;
    private final DataCollection dataCollection;
    private final Map<String, Sampler<?>> samplerMap;
    private final UpdateMethod updateMethod;

    public static final class GibbsBuilder {
        private final ParameterizedState state;
        private final DataCollection dataCollection;
        private final Map<String, Sampler<?>> samplerMap = new HashMap<>();

        public GibbsBuilder(final ParameterizedState state, final DataCollection dataCollection) {
            this.state = state;
            this.dataCollection = dataCollection;
        }

        public <T> GibbsBuilder addParameterSampler(final String parameterName, final Sampler<T> sampler) {
            if (samplerMap.containsKey(parameterName)) {
                throw new UnsupportedOperationException("Cannot add more than one sampler per parameter.");
            }
            if (!state.parameterNames().contains(parameterName)) {
                throw new UnsupportedOperationException("Cannot add sampler for parameter not specified in initial state.");
            }
            samplerMap.put(parameterName, sampler);
            return this;
        }

        public ParameterizedModel build() {
            if (!samplerMap.keySet().containsAll(state.parameterNames())) {
                throw new UnsupportedOperationException("Each parameter must have a corresponding sampler specified.");
            }
            return new ParameterizedModel(this);
        }
    }

    private ParameterizedModel(final GibbsBuilder builder) {
        this.state = builder.state;
        this.dataCollection = builder.dataCollection;
        this.samplerMap = builder.samplerMap;
        this.updateMethod = UpdateMethod.GIBBS;
    }

    public ParameterizedState state() {
        return state.copy();
    }

    public void update(final RandomGenerator rng) {
        if (updateMethod == UpdateMethod.GIBBS) {
            doGibbsUpdate(rng);
        }
    }

    private void doGibbsUpdate(final RandomGenerator rng) {
        for (final String parameterName : state.parameterNames()) {
            state.updateParameter(parameterName, samplerMap.get(parameterName).sample(rng, state, dataCollection));
        }
    }
}
