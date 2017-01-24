package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

import org.broadinstitute.barclay.utils.Utils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.DataLine;

import java.util.*;
import java.util.stream.Collector;

/**
 * Stores raw read count for a single interval
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public class RawReadCountData extends ReadCountData {

    private Target target;
    private int count;
    private final String RAW_COLUMN = "RAW";

    public RawReadCountData(Target target) {
        this.target = Utils.nonNull(target, "Target cannot be null. ");
        count = 0;
    }

    @Override
    public Target getTarget() {
        return target;
    }

    @Override
    public void updateReadCount(GATKRead read) {
        count++;
    }

    @Override
    public void appendCountsTo(DataLine dataLine) {
        dataLine.append(count);
    }

    @Override
    public List<String> getReadCountDataColumns() {
        return Arrays.asList(RAW_COLUMN);
    }
}
