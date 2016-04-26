package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.*;
import org.fusesource.leveldbjni.All;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Simple data structure to pass and read/write a List of {@link AllelicCount} objects.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class AllelicCountCollection {
    private final List<AllelicCount> counts;

    public AllelicCountCollection() {
        counts = new ArrayList<>();
    }

    /**
     * Constructor from from file. Checks whether the input file has just the basic columns, or full columns.
     * @param inputFile file to read from
     */
    public AllelicCountCollection(final File inputFile) {
        Utils.nonNull(inputFile);
        Utils.regularReadableUserFile(inputFile);

        try (final TableReader<AllelicCount> reader = TableUtils.reader(inputFile,
                (columns, formatExceptionFactory) -> {
                    /* in the order of decrasing verbosity */
                    if (columns.containsAll(AllelicCountTableColumns.FULL_COLUMN_NAME_ARRAY)) {
                        return AllelicCount.fullParser;
                    }
                    if (columns.containsAll(AllelicCountTableColumns.INTERMEDIATE_COLUMN_NAME_ARRAY)) {
                        return AllelicCount.intermediateParser;
                    }
                    if (columns.containsAll(AllelicCountTableColumns.BASIC_COLUMN_NAME_ARRAY)) {
                        return AllelicCount.basicParser;
                    }
                    /* if we are here, the table is malformed */
                    throw formatExceptionFactory.apply("The AllelicCountCollection input file does not have the" +
                            " mandatory colums. Required: " + StringUtils.join(
                            AllelicCountTableColumns.BASIC_COLUMN_NAME_ARRAY, ", "));
                })) {
            counts = reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        }
    }

    /**
     * Adds a new {@link AllelicCount} to counts.
     */
    public void add(final AllelicCount allelicCount) {
        counts.add(Utils.nonNull(allelicCount));
    }

    /** Returns an unmodifiable view of the list of AllelicCounts.   */
    public List<AllelicCount> getCounts() {
        return Collections.unmodifiableList(counts);
    }

    /**
     * @return a map from SimpleIntervals in counts to their integer index
     */
    public Map<SimpleInterval, Integer> getSimpleIntervalToIndexMap() {
        return IntStream.range(0, counts.size()).boxed()
                .collect(Collectors.toMap(i -> counts.get(i).getInterval(), i -> i));
    }

    /**
     * Writes out basic pulldown data (sequence, position, reference count, alternate count) to specified file.
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     */
    public void writeBasicCollection(final File outputFile) {
        try (final TableWriter<AllelicCount> writer = TableUtils.writer(outputFile,
                new TableColumnCollection(AllelicCountTableColumns.BASIC_COLUMN_NAME_ARRAY),
                (count, dataLine) -> {
                    final SimpleInterval interval = count.getInterval();
                    final int refReadCount = count.getRefReadCount();
                    final int altReadCount = count.getAltReadCount();
                    dataLine.append(interval.getContig())
                            .append(interval.getEnd())
                            .append(refReadCount)
                            .append(altReadCount);
                })) {
            writer.writeAllRecords(counts);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    /**
     * Writes out intermediate pulldown data (sequence, position, reference count, alternate count, reference nucleotide,
     * alternate nucleotide) to specified file.
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     */
    public void writeIntermediateCollection(final File outputFile) {
        try (final TableWriter<AllelicCount> writer = TableUtils.writer(outputFile,
                new TableColumnCollection(AllelicCountTableColumns.INTERMEDIATE_COLUMN_NAME_ARRAY),
                (count, dataLine) -> {
                    final SimpleInterval interval = count.getInterval();
                    final int refReadCount = count.getRefReadCount();
                    final int altReadCount = count.getAltReadCount();
                    final Nucleotide refNucleotide = count.getRefNucleotide();
                    final Nucleotide altNucleotide = count.getAltNucleotide();
                    final int readDepth = count.getReadDepth();
                    dataLine.append(interval.getContig())
                            .append(interval.getEnd())
                            .append(refReadCount)
                            .append(altReadCount)
                            .append(refNucleotide.name())
                            .append(altNucleotide.name())
                            .append(readDepth);
                })) {
            writer.writeAllRecords(counts);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    /**
     * Writes out full pulldown data (sequence, position, reference count, alternate count, reference nucleotide,
     * alternate nucleotide, het log odds) to specified file.
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     */
    public void writeFullCollection(final File outputFile) {
        try (final TableWriter<AllelicCount> writer = TableUtils.writer(outputFile,
                new TableColumnCollection(AllelicCountTableColumns.FULL_COLUMN_NAME_ARRAY),
                (count, dataLine) -> {
                    final SimpleInterval interval = count.getInterval();
                    final int refReadCount = count.getRefReadCount();
                    final int altReadCount = count.getAltReadCount();
                    final Nucleotide refNucleotide = count.getRefNucleotide();
                    final Nucleotide altNucleotide = count.getAltNucleotide();
                    final int readDepth = count.getReadDepth();
                    final double hetLogOdds = count.getHetLogOdds();
                    dataLine.append(interval.getContig())
                            .append(interval.getEnd())
                            .append(refReadCount)
                            .append(altReadCount)
                            .append(refNucleotide.name())
                            .append(altNucleotide.name())
                            .append(readDepth)
                            .append(hetLogOdds);
                })) {
            writer.writeAllRecords(counts);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof AllelicCountCollection)) {
            return false;
        }

        final AllelicCountCollection allelicCountCollection = (AllelicCountCollection) o;
        return counts.equals(allelicCountCollection.counts);
    }

    @Override
    public int hashCode() {
        return counts.hashCode();
    }
}
