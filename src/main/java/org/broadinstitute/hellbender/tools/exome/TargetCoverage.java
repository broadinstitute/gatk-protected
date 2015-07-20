package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

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

<<<<<<< HEAD
=======
    /**
     * read a list of targets with coverage from a file
     */
    public static List<TargetCoverage> readTargetsWithCoverage(final File file) throws IOException {
        try (final TableReader<TargetCoverage> reader = TableUtils.reader(file,
                (columns, formatExceptionFactory) -> {
                    if (!columns.matchesAll(0, "CONTIG", "START", "END", "NAME"))   //coverage is fifth column w/ header = <sample name>
                        throw formatExceptionFactory.apply("Bad header");

                    // return the lambda to translate dataLines into targets.
                    return (dataLine) -> new TargetCoverage(dataLine.get(3),
                            new SimpleInterval(dataLine.get(0), dataLine.getInt(1), dataLine.getInt(2)),
                            dataLine.getDouble(4));
                })) {
            return reader.stream().collect(Collectors.toList());
        } catch (UncheckedIOException e) {
            throw e.getCause();
        }
    }
>>>>>>> 5a9f284... Port of ReCapSeg Caller
}
