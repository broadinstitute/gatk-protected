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

/**
 * Created by slee on 09/10/15.
 */
public class Data<T extends Number> {
    private final List<T> data;

    public Data(final List<T> data) {
        Utils.nonNull(data, "The data cannot be null.");
        this.data = new ArrayList<>(data);
    }

    public int size() {
        return data.size();
    }

    public List<T> values() {
        return Collections.unmodifiableList(data);
    }

    public T value(final int index) {
        return data.get(index);
    }

    public static <T extends Number> Data<T> loadData(final File file, final Function<String, T> parse) {
        final List<T> list = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = br.readLine()) != null) {
                list.add(parse.apply(line));
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(file, e);
        }
        return new Data<>(list);
    }
}
