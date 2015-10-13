package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

public final class ParameterizedState {
    private final Map<String, Parameter<?>> parameterMap = new HashMap<>();

    public ParameterizedState(final List<Parameter<?>> parameters) {
        Utils.nonNull(parameters, "List of parameters cannot be null.");
        for (final Parameter<?> parameter : parameters) {
            if (parameterMap.containsKey(parameter.name())) {
                throw new IllegalArgumentException("List of parameters cannot contain duplicate parameter names.");
            }
            parameterMap.put(parameter.name(), parameter);
        }
    }

    public ParameterizedState(final ParameterizedState state) {
        this(new ArrayList<>(state.parameterMap.values()));
    }

    public <T> T get(final String parameterName, final Class<T> type) {
        try {
            return type.cast(parameterMap.get(parameterName).value());
        } catch (final NullPointerException | ClassCastException e) {
            if (e instanceof NullPointerException) {
                throw new IllegalArgumentException("Can only get pre-existing parameters; check parameter name.");
            }
            throw new IllegalArgumentException("Specified type of parameter does not match pre-existing type.");
        }
    }

    protected ParameterizedState copy() {
        return new ParameterizedState(this);
    }

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
}
