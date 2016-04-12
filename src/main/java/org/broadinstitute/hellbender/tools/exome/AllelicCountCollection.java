package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.collections.CollectionUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Simple data structure to pass and read/write a List of {@link AllelicCount} objects.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AllelicCountCollection {
    private final List<AllelicCount> counts;

    public AllelicCountCollection() {
        counts = new ArrayList<>();
    }

    /**
     * Constructor that reads (sequence, position, reference count, alternate count) from the specified file.
     * @param inputFile     file to read from
     */
    public AllelicCountCollection(final File inputFile) {
        Utils.nonNull(inputFile);
        Utils.regularReadableUserFile(inputFile);
        try (final TableReader<AllelicCount> reader = TableUtils.reader(inputFile,
                (columns, formatExceptionFactory) -> {
                    TableUtils.checkMandatoryColumns(columns, AllelicCountTableColumns.COLUMN_NAME_ARRAY, formatExceptionFactory);

                    // return the lambda to translate dataLines into AllelicCounts.
                    return (dataLine) -> {
                        final int position = dataLine.getInt(AllelicCountTableColumns.POSITION.toString());
                        final SimpleInterval interval = new SimpleInterval(dataLine.get(AllelicCountTableColumns.CONTIG.toString()), position, position);
                        final int refReadCount = dataLine.getInt(AllelicCountTableColumns.REF_COUNT.toString());
                        final int altReadCount = dataLine.getInt(AllelicCountTableColumns.ALT_COUNT.toString());
                        return new AllelicCount(interval, refReadCount, altReadCount);
                    };
                })) {
            counts = reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        }
    }

    /**
     * Adds (interval, reference count, alternate count) to respective lists.
     * @param interval      site in 1-based format
     * @param refReadCount  number of reads at site matching the reference
     * @param altReadCount  number of reads at site different from the reference
     */
    public void add(final SimpleInterval interval, final int refReadCount, final int altReadCount) {
        counts.add(new AllelicCount(interval, refReadCount, altReadCount));
    }

    /**
     * Adds (interval, reference count, alternate count, reference nucleotide, alternate nucleotide) to respective lists.
     * @param interval      site in 1-based format
     * @param refReadCount  number of reads at site matching the reference
     * @param altReadCount  number of reads at site different from the reference
     * @param refNucleotide the ref nucleotide
     * @param altNucleotide the alt nucleotide
     */
    public void add(final SimpleInterval interval, final int refReadCount, final int altReadCount,
                    final Nucleotide refNucleotide, final Nucleotide altNucleotide) {
        counts.add(new AllelicCount(interval, refReadCount, altReadCount, refNucleotide, altNucleotide));
    }

    /** Returns an unmodifiable view of the list of AllelicCounts.   */
    public List<AllelicCount> getCounts() {
        return Collections.unmodifiableList(counts);
    }

    /**
     * @return a map from SimpleIntervals in counts to their integer index
     */
    public Map<SimpleInterval, Integer> getSimpleIntervalToIndexMap() {
        Map<SimpleInterval, Integer> intervalToIndexMap = new HashMap<>();
        IntStream.range(0, counts.size()).forEach(i -> intervalToIndexMap.put(counts.get(i).getInterval(), i));
        return intervalToIndexMap;
    }

    /**
     * Writes out (sequence, position, reference count, alternate count) to specified file.
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     */
    public void write(final File outputFile) {
        try (final TableWriter<AllelicCount> writer = TableUtils.writer(outputFile,
                new TableColumnCollection(AllelicCountTableColumns.COLUMN_NAME_ARRAY),
                //lambda for filling an initially empty DataLine
                (count, dataLine) -> {
                    final SimpleInterval interval = count.getInterval();
                    final int refReadCount = count.getRefReadCount();
                    final int altReadCount = count.getAltReadCount();
                    dataLine.append(interval.getContig()).append(interval.getEnd(), refReadCount, altReadCount);
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
