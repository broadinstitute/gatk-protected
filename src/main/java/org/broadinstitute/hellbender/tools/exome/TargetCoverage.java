package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Exome analysis target with coverage
 *
 * <p>A name entity with genomic interval and coverage</p>
 *
 * @author David Benjamin
 */

public final class TargetCoverage extends Target {
    private double coverage;

    /**
     * Construct a new TargetCoverage given all its properties.
     * @param interval the interval.
     * @param name the name of the interval.
     * @param coverage the coverage (read count or normalized coverage) of this interval
     *
     * @throws IllegalArgumentException if either {@code interval} or {@code start}
     *    and {@code end} do not represent a valid interval
     */
    public TargetCoverage(final String name, final SimpleInterval interval, final double coverage) {
        super(name, interval);
        this.coverage = coverage;
    }

    /**
     * Returns the coverage.
     */
    public double getCoverage() {
        return coverage;
    }

    /**
     * Sets the coverage.
     */
    public void setCoverage(final double coverage) {
        this.coverage = coverage;
    }

    /**
     * write a list of targets with coverage to file
     */
    public static void writeTargetsWithCoverage(final File outFile, final String sampleName,
                                                List<TargetCoverage> targets) throws IOException {
        try (final TableWriter<TargetCoverage> writer = TableUtils.writer(outFile,
                new TableColumnCollection("name", "contig", "start", "stop", sampleName),

                //lambda for filling an initially empty DataLine
                (target, dataLine) -> {
                    final String name = target.getName();
                    final SimpleInterval interval = target.getInterval();
                    final String contig = interval.getContig();
                    final int start = interval.getStart();
                    final int end = interval.getEnd();
                    final double coverage = target.getCoverage();
                    dataLine.append(name).append(contig).append(start, end).append(coverage);
                })) {
            for (final TargetCoverage target : targets) {
                writer.writeRecord(target);
            }
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof TargetCoverage)) {
            return false;
        }
        if (!super.equals(o)) {
            return false;
        }

        TargetCoverage that = (TargetCoverage) o;
        return getInterval().equals(that.getInterval()) && Math.abs(coverage - that.coverage) < 0.0000000000001;
    }
}
