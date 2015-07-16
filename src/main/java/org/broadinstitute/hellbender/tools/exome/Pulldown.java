package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Simple data structure to pass and write pulldown results.  Should probably replace with a more generic class later.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class Pulldown {
    private final SAMFileHeader header;
    private final List<AllelicCount> counts;

    public Pulldown(final SAMFileHeader header) {
        Utils.nonNull(header, "SAMFileHeader must be supplied.");

        this.header = header;
        this.counts = new ArrayList<>();
    }

    /**
     * Constructor that reads (sequence, position, reference count, alternate count) from the specified file and
     * uses external SAMFile header to construct Pulldown.
     * @param inputFile     file to read from
     * @param header        SAMFile header for IntervalList
     * TODO remove dependency from IntervalList/SamLocusIterator on external header once LocusWalker implemented?
     */
    public Pulldown(final File inputFile, final SAMFileHeader header) {
        this(header);

        try (final TableReader<AllelicCount> reader = TableUtils.reader(inputFile,
                (columns, formatExceptionFactory) -> {
                    if (!columns.matchesExactly("SEQ", "POS", "REF_COUNT", "ALT_COUNT"))
                        throw formatExceptionFactory.apply("Bad header");

                    // return the lambda to translate dataLines into Pulldown rows.
                    return (dataLine) -> {
                        final Interval interval = new Interval(dataLine.get(0), dataLine.getInt(1), dataLine.getInt(1));
                        final int refReadCount = dataLine.getInt(2);
                        final int altReadCount = dataLine.getInt(3);
                        return new AllelicCount(interval, refReadCount, altReadCount);
                    };
                })) {
            for (final AllelicCount count : reader) {
                this.counts.add(count);
            }
        } catch (final Exception e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e.getMessage());
        }
    }

    /**
     * Adds (interval, reference count, alternate count) to respective lists.
     * @param interval      heterozygous SNP site in 1-based format
     * @param refReadCount  number of reads at SNP site matching the reference
     * @param altReadCount  number of reads at SNP site different from the reference
     */
    public void add(final Interval interval, final int refReadCount, final int altReadCount) {
        final AllelicCount count = new AllelicCount(interval, refReadCount, altReadCount);
        counts.add(count);
    }

    /** Returns the IntervalList of SNP sites.   */
    public IntervalList getIntervals() {
        final IntervalList intervals = new IntervalList(header);
        intervals.addall(counts.stream().map(count -> count.getInterval()).collect(Collectors.toList()));
        return intervals;
    }

    /** Returns a list of the allelic counts.   */
    public List<AllelicCount> asList() {
        return new ArrayList<>(counts);
    }

    /**
     * Writes out (sequence, position, reference count, alternate count) to specified file.
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     */
    public void write(final File outputFile) {
        try (final TableWriter<AllelicCount> writer = TableUtils.writer(outputFile,
                new TableColumnCollection("SEQ", "POS", "REF_COUNT", "ALT_COUNT"),
                //lambda for filling an initially empty DataLine
                (count, dataLine) -> {
                    final Interval interval = count.getInterval();
                    final int refReadCount = count.getRefReadCount();
                    final int altReadCount = count.getAltReadCount();
                    dataLine.append(interval.getContig()).append(interval.getEnd(), refReadCount, altReadCount);
                })) {
            for (final AllelicCount count : counts) {
                writer.writeRecord(count);
            }
        } catch (final Exception e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e.getMessage());
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof Pulldown)) {
            return false;
        }

        final Pulldown pulldown = (Pulldown) o;
        return header.equals(pulldown.header) && counts.equals(pulldown.counts);
    }

    @Override
    public int hashCode() {
        int result = header.hashCode();
        result = 31 * result + counts.hashCode();
        return result;
    }
}
