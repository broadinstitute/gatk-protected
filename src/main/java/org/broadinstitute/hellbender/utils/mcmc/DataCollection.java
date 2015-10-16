package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Represents a collection of named datasets, each consisting of numeric datapoints and represented
 * by a {@link Data} object.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class DataCollection {
    private final Map<String, Data<?>> datasetMap = new HashMap<>();

    /**
     * Constructs a DataCollection object given a Collection of datasets represented by {@link Data} objects,
     * which should be well-formed and uniquely named.
     * @param datasets  Collection of Data objects
     * @throws IllegalArgumentException if {@code datasets} is {@code null},
     *                                  or if not of all of the datasets are uniquely named
     */
    public DataCollection(final Collection<Data<?>> datasets) {
        Utils.nonNull(datasets, "The collection of datasets cannot be null.");
        for (final Data<?> dataset : datasets) {
            if (datasetMap.containsKey(dataset.name())) {
                throw new IllegalArgumentException("Each dataset in the collection must have a unique name.");
            }
            datasetMap.put(dataset.name(), dataset);
        }
    }

    /**
     * Constructs an empty DataCollection.
     */
    public DataCollection() {
        this(Collections.emptyList());
    }

    /**
     * Returns the size of the DataCollection.
     * @return  size of the DataCollection
     */
    protected int size() {
        return datasetMap.values().size();
    }

    /**
     * Adds a dataset represented by a {@link Data} object to the DataCollection.
     * @param dataset   Data object to be added
     */
    protected void add(final Data<?> dataset) {
        if (datasetMap.containsKey(dataset.name())) {
            throw new IllegalArgumentException("Cannot add a dataset with the same name as another dataset in the collection.");
        }
        datasetMap.put(dataset.name(), dataset);
    }

    /**
     * Returns an unmodifiable view of the List of the values of a dataset in the DataCollection, given the name
     * of the dataset.
     * @param datasetName   name of the dataset to be retrieved
     * @param <T>           type of the datapoints to be retrieved
     * @return              unmodifiable view of the List of the datapoints of the named dataset
     * @throws IllegalArgumentException if no dataset with name given {@code datasetName} is present in the DataCollection
     * @throws UnsupportedOperationException if the type {@code <T>} does not match the type of the datapoints
     *                                       in the named dataset
     */
    @SuppressWarnings("unchecked")
    public <T> List<T> get(final String datasetName) {
        try {
            return (List<T>) Collections.unmodifiableList(datasetMap.get(datasetName).values());
        } catch (final NullPointerException | ClassCastException e) {
            if (e instanceof NullPointerException) {
                throw new IllegalArgumentException("Can only get pre-existing datasets; check dataset name.");
            }
            throw new UnsupportedOperationException("Type of dataset specified in getter does not match pre-existing type.");
        }
    }
}
