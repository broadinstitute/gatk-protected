package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by davidben on 7/17/15.
 */
public final class TargetCoverageUtils {
    private static final String TARGET_NAME_COLUMN = "name";
    private static final String CONTIG_COLUMN = "contig";
    private static final String START_COLUMN = "start";
    private static final String END_COLUMN = "stop";

    /**
     * read a list of targets with coverage from a tab-separated file with header line:
     * name contig  start   stop    <sample name>
     * where the final column is numerical data, generally coverage or normalized coverage
     */
    public static List<TargetCoverage> readTargetsWithCoverage(final File targetsFile) {
        try (final TableReader<TargetCoverage> reader = TableUtils.reader(targetsFile,
                (columns, formatExceptionFactory) -> {
                    if (!columns.containsAll(TARGET_NAME_COLUMN, CONTIG_COLUMN, START_COLUMN, END_COLUMN))
                        throw formatExceptionFactory.apply("Bad header");
                    //return the lambda to translate dataLines into targets
                    return (dataLine) -> new TargetCoverage(dataLine.get(TARGET_NAME_COLUMN),
                            new SimpleInterval(dataLine.get(CONTIG_COLUMN), dataLine.getInt(START_COLUMN), dataLine.getInt(END_COLUMN)), dataLine.getDouble(4));
                })) {
            return reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(targetsFile, e);
        }
    }

    /**
     * read a list of targets with coverage from a file into a TargetCollection
     */
    public static TargetCollection<TargetCoverage> readTargetsWithCoverageIntoTargetCollection(final File file) throws IOException {
        final List<TargetCoverage> targetList = readTargetsWithCoverage(file);
        return new HashedListTargetCollection<>(targetList);
    }

    public static TargetCollection<TargetCoverage> readModeledTargetFileIntoTargetCollection(final File file) {
        try (final TableReader<TargetCoverage> reader = TableUtils.reader(file,
                (columns, formatExceptionFactory) -> {
                    if (!columns.containsAll(TARGET_NAME_COLUMN, CONTIG_COLUMN, START_COLUMN, END_COLUMN) || (columns.columnCount() < 5))

                        throw formatExceptionFactory.apply("Bad header");
                    //return the lambda to translate dataLines into targets
                    //coverage is fifth column w/ header = <sample name>, so we use the column index.
                    return (dataLine) -> new TargetCoverage(dataLine.get(TARGET_NAME_COLUMN),
                            new SimpleInterval(dataLine.get(CONTIG_COLUMN), dataLine.getInt(START_COLUMN), dataLine.getInt(END_COLUMN)),
                            dataLine.getDouble(4));
                })) {
            return new HashedListTargetCollection<>(reader.stream().collect(Collectors.toList()));
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(file, e);
        }
    }

    public static TargetCollection<TargetCoverage> readReadCountCollectionIntoTargetCollection(final File file) {
        try (final TableReader<TargetCoverage> reader = TableUtils.reader(file,
                (columns, formatExceptionFactory) -> {
                    if (!columns.containsAll(ExomeReadCounts.EXON_NAME_COLUMN_NAME, ExomeReadCounts.EXON_CONTIG_COLUMN_NAME, ExomeReadCounts.EXON_START_COLUMN_NAME,
                            ExomeReadCounts.EXON_END_COLUMN_NAME) || (columns.columnCount() < 5))

                        throw formatExceptionFactory.apply("Bad header");
                    //return the lambda to translate dataLines into targets
                    //coverage is fifth column w/ header = <sample name>, so we use the column index.
                    return (dataLine) -> new TargetCoverage(dataLine.get(ExomeReadCounts.EXON_NAME_COLUMN_NAME),
                            new SimpleInterval(dataLine.get(ExomeReadCounts.EXON_CONTIG_COLUMN_NAME),
                                    dataLine.getInt(ExomeReadCounts.EXON_START_COLUMN_NAME), dataLine.getInt(ExomeReadCounts.EXON_END_COLUMN_NAME)),
                            dataLine.getDouble(4));
                })) {
            return new HashedListTargetCollection<>(reader.stream().collect(Collectors.toList()));
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(file, e);
        }
    }


    /**
     * write a list of targets with coverage to file
     */
    public static void writeTargetsWithCoverage(final File outFile, final String sampleName,
                                                final List<TargetCoverage> targets) {
        try (final TableWriter<TargetCoverage> writer = TableUtils.writer(outFile,
                new TableColumnCollection(TARGET_NAME_COLUMN, CONTIG_COLUMN, START_COLUMN, END_COLUMN, sampleName),
                //lambda for filling an initially empty DataLine
                (target, dataLine) -> {
                    final SimpleInterval interval = target.getInterval();
                    final String contig = interval.getContig();
                    final int start = interval.getStart();
                    final int end = interval.getEnd();
                    final String name = target.getName();
                    final double coverage = target.getCoverage();
                    dataLine.append(name, contig).append(start, end).append(coverage);
                })) {
            for (final TargetCoverage target : targets) {
                writer.writeRecord(target);
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e);
        }
    }
}
