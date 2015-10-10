package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

public class ParameterizedState {
    private final Map<String, Parameter> parameterMap = new HashMap<>();

    public ParameterizedState(final List<Parameter> parameters) {
        Utils.nonNull(parameters, "List of parameters cannot be null.");
        for (final Parameter parameter : parameters) {
            if (parameterMap.containsKey(parameter.name())) {
                throw new IllegalArgumentException("List of parameters cannot contain duplicate parameter names.");
            }
            parameterMap.put(parameter.name(), parameter);
        }
    }

    public Parameter getParameter(final String parameterName) {
        if (!parameterMap.containsKey(parameterName)) {
            throw new UnsupportedOperationException("Can only get pre-existing parameters.");
        }
        return parameterMap.get(parameterName);
    }

    public ParameterizedState copy() {
        return new ParameterizedState(new ArrayList<>(parameterMap.values()));
    }

    public Set<String> parameterNames() {
        return parameterMap.keySet();
    }

    public void updateParameter(final String parameterName, final Parameter parameter) {
        if (!parameterMap.containsKey(parameterName)) {
            throw new UnsupportedOperationException("Can only update pre-existing parameters.");
        }
        if (parameterName != parameter.name()) {
            throw new UnsupportedOperationException("Parameter names of sampler and sample must match.");
        }
        parameterMap.put(parameterName, parameter);
    }
}
