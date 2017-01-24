package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.DataLine;

import java.util.List;

/**
 * A container that stores and updates read count data pertaining to a single interval
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */

public abstract class ReadCountData {

    /**
     * Updates information about the interval's read count data given a new read
     * @param read GATK read
     */
    public abstract void updateReadCount(GATKRead read);

    /**
     * Appends this instance read counts to a data-line object.
     *
     * @param dataLine the destination data-line.
     * @throws IllegalArgumentException if {@code dataLine is {@code null}}.
     * @throws IllegalStateException    if there is not enough room in {@code dataLine} from its current appending position
     *                                  to add all the record counts.
     */
    public abstract void appendCountsTo(final DataLine dataLine);

    /**
     * Get the target of this record
     *
     * @return never {@code null}.
     */
    public abstract Target getTarget();

    /**
     * Get the columns of the read count record.
     * Classes that extend {@link ReadCountData} should implement a concrete version of this method
     * @return columns of the read count data
     */
    public List<String> getReadCountDataColumns() { return null; }

}
