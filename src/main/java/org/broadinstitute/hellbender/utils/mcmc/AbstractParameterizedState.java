package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

public abstract class AbstractParameterizedState {
    private final Map<String, Parameter<?>> parameterMap = new HashMap<>();

    public AbstractParameterizedState(final List<Parameter<?>> parameters) {
        Utils.nonNull(parameters, "List of parameters cannot be null.");
        for (final Parameter<?> parameter : parameters) {
            if (parameterMap.containsKey(parameter.name())) {
                throw new IllegalArgumentException("List of parameters cannot contain duplicate parameter names.");
            }
            parameterMap.put(parameter.name(), parameter);
        }
    }

    public AbstractParameterizedState(final AbstractParameterizedState state) {
        this(new ArrayList<>(state.parameterMap.values()));
    }

    public <T> AbstractParameterizedState(final String parameterNamePrefix, final List<T> parameterValues) {
        this(initializeParameters(parameterNamePrefix, parameterValues));
    }

    public <T> T get(final String parameterName, final Class<T> parameterValueClass) {
        try {
            return parameterValueClass.cast(parameterMap.get(parameterName).value());
        } catch (final NullPointerException | ClassCastException e) {
            if (e instanceof NullPointerException) {
                throw new IllegalArgumentException("Can only get pre-existing parameters; check parameter name.");
            }
            throw new IllegalArgumentException("Type of parameter specified in getter does not match pre-existing type.");
        }
    }

    protected abstract <S extends AbstractParameterizedState> S copy(final Class<S> stateClass);

    protected Set<String> parameterNames() {
        return parameterMap.keySet();
    }

    protected <T> void updateParameter(final String parameterName, final T value) {
        try {
            parameterMap.put(parameterName, new Parameter<>(parameterName, value));
        } catch (final NullPointerException e) {
            throw new UnsupportedOperationException("Can only update pre-existing parameters; check parameter name.");
        }
    }

    private static <T> List<Parameter<?>> initializeParameters(final String parameterNamePrefix,
                                                               final List<T> parameterValues) {
        final List<Parameter<?>> initialParameters = new ArrayList<>();
        for (int i = 0; i < parameterValues.size(); i++) {
            initialParameters.add(new Parameter<>(parameterNamePrefix + i, parameterValues.get(i)));
        }
        return initialParameters;
    }
}
