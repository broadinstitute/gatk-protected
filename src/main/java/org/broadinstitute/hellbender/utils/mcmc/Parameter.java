package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.Utils;

public final class Parameter<T> {
    private final String name;
    private final T value;

    public Parameter(final String name, final T value) {
        Utils.nonNull(name, "The parameter name cannot be null.");
        Utils.nonNull(value, "The parameter value cannot be null.");
        this.name = name;
        this.value = value;
    }

    public String name() {
        return name;
    }

    public T value() {
        return value;
    }
}
