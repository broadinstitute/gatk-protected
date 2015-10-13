package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public final class DataCollection {
    private final Map<String, Data<?>> datasetMap = new HashMap<>();

    public DataCollection(final List<Data<?>> datasets) {
        Utils.nonNull(datasets, "The list of datasets cannot be null.");
        for (final Data<?> dataset : datasets) {
            if (datasetMap.containsKey(dataset.name())) {
                throw new IllegalArgumentException("The list of datasets cannot contain duplicate parameter names.");
            }
            datasetMap.put(dataset.name(), dataset);
        }
    }

    @SuppressWarnings("unchecked")
    public <T> List<T> get(final String datasetName) {
        try {
            return (List<T>) datasetMap.get(datasetName).values();
        } catch (final NullPointerException | ClassCastException e) {
            if (e instanceof NullPointerException) {
                throw new IllegalArgumentException("Can only get pre-existing datasets; check dataset name.");
            }
            throw new UnsupportedOperationException("Type specified in method type parameter does not match pre-existing type of dataset.");
        }
    }
}
