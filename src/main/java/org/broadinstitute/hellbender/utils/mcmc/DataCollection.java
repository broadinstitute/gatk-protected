package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class DataCollection {
    private final List<Data> collection;

    public DataCollection(final List<Data> collection) {
        Utils.nonNull(collection, "The data collection cannot be null.");
        this.collection = new ArrayList<>(collection);
    }

    public List<Data> get() {
        return Collections.unmodifiableList(collection);
    }

    public Data get(final int index) {
        return collection.get(index);
    }
}
