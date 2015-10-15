package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

public class DataCollection {
    private final Map<String, Data<?>> datasetMap = new HashMap<>();

    public DataCollection(final Collection<Data<?>> datasets) {
        Utils.nonNull(datasets, "The collection of datasets cannot be null.");
        for (final Data<?> dataset : datasets) {
            if (datasetMap.containsKey(dataset.name())) {
                throw new IllegalArgumentException("Each dataset in the collection must have a unique name.");
            }
            datasetMap.put(dataset.name(), dataset);
        }
    }

    public DataCollection() {
        this(Collections.emptyList());
    }

    protected int size() {
        return datasetMap.values().size();
    }

    protected void add(final Data<?> dataset) {
        if (datasetMap.containsKey(dataset.name())) {
            throw new IllegalArgumentException("Cannot add a dataset with the same name as another dataset in the collection.");
        }
        datasetMap.put(dataset.name(), dataset);
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
