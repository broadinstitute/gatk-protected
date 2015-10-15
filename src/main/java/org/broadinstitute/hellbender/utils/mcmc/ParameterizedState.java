package org.broadinstitute.hellbender.utils.mcmc;

import java.util.List;

public final class ParameterizedState extends AbstractParameterizedState {
    public ParameterizedState(final List<Parameter<?>> parameters) {
        super(parameters);
    }

    public ParameterizedState(final ParameterizedState state) {
        super(state);
    }

    @Override
    protected <S extends AbstractParameterizedState> S copy(final Class<S> stateClass) {
        return stateClass.cast(new ParameterizedState(this));
    }
}
