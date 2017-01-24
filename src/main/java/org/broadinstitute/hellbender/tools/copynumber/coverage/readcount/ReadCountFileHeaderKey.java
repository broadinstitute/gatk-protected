package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.barclay.utils.Utils;

import org.apache.commons.lang3.tuple.Pair;

import java.util.Arrays;
import java.util.Optional;

/**
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public enum ReadCountFileHeaderKey {
    SAMPLE_NAME("sampleName"),
    READ_COUNT_TYPE("readCountType");

    private final String headerKeyName;

    // this regex expression is used to extract the value of the corresponding key from the read counts file header
    private final static String headerKeyValuePattern = "^[# ]*(\\S*)[\\s\t]*=[\\s\t]*(.*\\S)[\\s\t]*$";

    ReadCountFileHeaderKey(final String name) {
        this.headerKeyName = Utils.nonNull(name);
    }

    public String getHeaderKeyName() {
        return headerKeyName;
    }


    /**
     * Parse the header line and return its key value pair if the key matches one of the {@link ReadCountFileHeaderKey}
     * instances
     * @param headerLine a single comment line
     * @return a pair object
     */
    public static Pair<ReadCountFileHeaderKey, String> getHeaderValueForKey(String headerLine) {
        String key = headerLine.replaceAll(headerKeyValuePattern, "$1");
        String value = headerLine.replaceAll(headerKeyValuePattern, "$2");

        Optional<ReadCountFileHeaderKey> optionalHeaderKey = Arrays.asList(values()).
                stream().filter(t -> t.getHeaderKeyName().equals(key)).findAny();
        return optionalHeaderKey.isPresent() ? new ImmutablePair<>(optionalHeaderKey.get(), value) : null;
    }
}
