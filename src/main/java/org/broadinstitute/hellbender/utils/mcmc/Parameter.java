package org.broadinstitute.hellbender.utils.mcmc;

public class Parameter<T> {
    private final String name;
    private final T value;

    public Parameter(final String name, final T value) {
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
