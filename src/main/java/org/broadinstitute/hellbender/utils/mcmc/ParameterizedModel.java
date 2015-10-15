package org.broadinstitute.hellbender.utils.mcmc;

import org.apache.commons.math3.random.RandomGenerator;

import java.lang.reflect.Type;
import java.util.HashMap;
import java.util.Map;

public final class ParameterizedModel<S1 extends ParameterizedState, T1 extends DataCollection> {
    private enum UpdateMethod {
        GIBBS
    }

    private final S1 state;
    private final T1 dataCollection;
    private final Map<String, Sampler<?, S1, T1>> samplerMap;
    private final UpdateMethod updateMethod;

    public static final class GibbsBuilder<S2 extends ParameterizedState, T2 extends DataCollection> {
        private final S2 state;
        private final T2 dataCollection;
        private final Map<String, Sampler<?, S2, T2>> samplerMap = new HashMap<>();

        public GibbsBuilder(final S2 state, final T2 dataCollection) {
            if (dataCollection.size() == 0) {
                throw new IllegalArgumentException("The collection of datasets cannot be empty.");
            }
            this.state = state;
            this.dataCollection = dataCollection;
        }

        public <U> GibbsBuilder<S2, T2> addParameterSampler(final String parameterName, final Sampler<U, S2, T2> sampler) {
            if (samplerMap.containsKey(parameterName)) {
                throw new UnsupportedOperationException("Cannot add more than one sampler per parameter.");
            }
            if (!state.parameterNames().contains(parameterName)) {
                throw new UnsupportedOperationException("Cannot add sampler for parameter not specified in initial state.");
            }
            samplerMap.put(parameterName, sampler);
            return this;
        }

        public ParameterizedModel<S2, T2> build() {
            if (!samplerMap.keySet().containsAll(state.parameterNames())) {
                throw new UnsupportedOperationException("Each parameter must have a corresponding sampler specified.");
            }
            return new ParameterizedModel<>(this);
        }
    }

    private ParameterizedModel(final GibbsBuilder<S1, T1> builder) {
        this.state = builder.state;
        this.dataCollection = builder.dataCollection;
        this.samplerMap = builder.samplerMap;
        this.updateMethod = UpdateMethod.GIBBS;
    }

    public S1 state(final Class<S1> stateClass) {
        return (S1) state.copy();
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
