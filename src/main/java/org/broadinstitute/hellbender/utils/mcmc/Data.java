package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

public final class Data<N extends Number> {
    private final String name;
    private final List<N> data;

    public Data(final String name, final List<N> data) {
        Utils.nonNull(name, "The name of the dataset cannot be null.");
        Utils.nonNull(data, "The dataset cannot be null.");
        if (data.size() == 0) {
            throw new IllegalArgumentException("The dataset cannot be empty.");
        }
        this.name = name;
        this.data = new ArrayList<>(data);
    }

    public Data(final String name, final File file, final Function<String, N> parse) {
        this(name, loadData(file, parse));
    }

    protected String name() {
        return name;
    }

    protected List<N> values() {
        return Collections.unmodifiableList(data);
    }

    protected N value(final int index) {
        return data.get(index);
    }

    protected int size() {
        return data.size();
    }

    private static <T> List<T> loadData(final File file, final Function<String, T> parse) {
        final List<T> list = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = br.readLine()) != null) {
                list.add(parse.apply(line));
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(file, e);
        }
        return list;
    }
}
