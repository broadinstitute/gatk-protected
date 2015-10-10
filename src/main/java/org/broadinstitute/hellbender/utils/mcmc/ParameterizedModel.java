package org.broadinstitute.hellbender.utils.mcmc;

import org.apache.commons.math3.random.RandomGenerator;

import java.util.HashMap;
import java.util.Map;

public class ParameterizedModel<S1 extends ParameterizedState, T1 extends DataCollection> {
    private enum UpdateMethod {
        GIBBS
    }

    private final S1 state;
    private final T1 dataCollection;
    private final Map<String, Sampler<?, S1, T1>> samplerMap;
    private final UpdateMethod updateMethod;

    public static final class GibbsBuilder<S2 extends ParameterizedState, T2 extends DataCollection> {
        private final S2 modelState;
        private final T2 dataCollection;
        private final Map<String, Sampler<?, S2, T2>> samplerMap = new HashMap<>();

        public GibbsBuilder(final S2 modelState, final T2 dataCollection) {
            this.modelState = modelState;
            this.dataCollection = dataCollection;
        }

        public GibbsBuilder<S2, T2> addParameterSampler(final String parameterName, final Sampler<Parameter, S2, T2> sampler) {
            if (samplerMap.containsKey(parameterName)) {
                throw new UnsupportedOperationException("Cannot add more than one sampler per parameter.");
            }
            samplerMap.put(parameterName, sampler);
            return this;
        }

        public ParameterizedModel<S2, T2> build() {
            if (!samplerMap.keySet().containsAll(modelState.parameterNames())) {
                throw new IllegalArgumentException("Each parameter must have a corresponding sampler specified.");
            }
            return new ParameterizedModel<>(this);
        }
    }

    private ParameterizedModel(final GibbsBuilder<S1, T1> builder) {
        this.state = builder.modelState;
        this.dataCollection = builder.dataCollection;
        this.samplerMap = builder.samplerMap;
        this.updateMethod = UpdateMethod.GIBBS;
    }

    public S1 state() {
        return (S1) state.copy();
    }

    public void doGibbsUpdate(final RandomGenerator rng) {
        if (updateMethod != UpdateMethod.GIBBS) {
            throw new UnsupportedOperationException("Update method must be set to Gibbs.");
        }
        for (final String parameterName : state.parameterNames()) {
            state.updateParameter(parameterName, (Parameter) samplerMap.get(parameterName).sample(rng, state, dataCollection));
        }
    }
}
