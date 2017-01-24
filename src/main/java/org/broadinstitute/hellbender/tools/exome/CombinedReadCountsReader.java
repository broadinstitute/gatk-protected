package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.apache.logging.log4j.Logger;
/**
 * Created by asmirnov on 12/5/16.
 */
public class CombinedReadCountsReader {

    private List<File> coverageFiles;

    private final int maxMergeFiles;

    private TargetCollection<Target> targets;
    //what table column to read from (PCOV or RAW)

    private TargetTableColumn transform;


    //stores combined read counts from all samples
    private double[][] readCounts;

    //stores sample names
    private final List<String> columnNames;

    private int lastColumnIndex;


    /**
     * TODO finish java doc
     * @param coverageFiles
     * @param maxMergeSize
     * @param targets
     * @param transform
     */
    public CombinedReadCountsReader(final List<File> coverageFiles, int maxMergeSize, TargetCollection<Target> targets, TargetTableColumn transform) {
        this.readCounts = new double[targets.targetCount()][coverageFiles.size()];
        this.coverageFiles = Utils.nonNull(coverageFiles, "input file list cannot be null");
        this.targets = Utils.nonNull(targets, "target collection cannot be null");
        this.maxMergeFiles = maxMergeSize;
        this.transform = transform;

        this.columnNames = new ArrayList<>();
        this.lastColumnIndex = 0;
        //TODO check if the parameter has one of the appropriate values(create function in TargetTableColumn)

    }

    public ReadCountCollection read(final Logger logger) {
        final Queue<File> remainingFilesToMerge = new ArrayDeque<>(coverageFiles);
        while(!remainingFilesToMerge.isEmpty()) {
            logger.debug(String.format("Merging %d of %d",
                    Math.max(maxMergeFiles, remainingFilesToMerge.size()), remainingFilesToMerge.size()));
            final List<File> filesToMerge = removeFilesToMergeNext(maxMergeFiles, remainingFilesToMerge);
            doMerge(filesToMerge);
            lastColumnIndex += filesToMerge.size();
        }
        return new ReadCountCollection(Collections.unmodifiableList(targets.targets()), columnNames, new Array2DRowRealMatrix(readCounts));
    }

    /**
     * Extracts from the remaining files-to-merge queue the ones to be merged next.
     * <p>
     * This method modifies the input queue {@code filesToMerge} by removing the files to be merged next.
     * </p>
     * @param maximumMergingFileCount the maximum merge file group size.
     * @param remainingFilesToMerge the file-to-merge queue.
     * @return never {@code null}. The return list won't ever have more than {@code maximumMergingFileCount} members.
     */
    private List<File> removeFilesToMergeNext(final int maximumMergingFileCount, final Queue<File> remainingFilesToMerge) {
        final List<File> result = new ArrayList<>(maximumMergingFileCount);
        for (int i = 0; i < maximumMergingFileCount && !remainingFilesToMerge.isEmpty(); i++) {
            result.add(remainingFilesToMerge.remove());
        }
        return result;
    }

    /**
     * Returns the list of count column names (no target info related columns)
     * in the order they appear in the table column collection.
     * @param columns the input table column collection.
     * @return never {@code null} but perhaps empty.
     */
    private List<String> readCountColumnNames(final TableColumnCollection columns) {
        return columns.names().stream()
                .filter(n -> !TargetTableColumn.isMandatoryTargetColumnName(n))
                .collect(Collectors.toList());
    }

    /**
     * The actual merge operation.
     * @param filesToMerge input files to be merged.
     */
    private void doMerge(final List<File> filesToMerge) {
        final ReadCountReaderCollection readers = new ReadCountReaderCollection(filesToMerge, targets);
        columnNames.addAll(readers.countColumnNames);
        for(final ReadCountRecord record: readers) {
            int idx;
            for(int i = 0; i < readers.countColumnSourceIndexMap.length; i++) {
                idx = readers.countColumnSourceIndexMap[i];
                readCounts[targets.index(record.getTarget().getName())][lastColumnIndex + idx] = record.getDouble(idx);
            }
        }
    }

    /**
     * TODO finish the java doc
     */
    private final class ReadCountsTableReader extends TableReader<ReadCountRecord.SingleSampleRecord> {
        private boolean hasName;
        private boolean hasCoordinates;
        private TargetTableColumn columnToRead;
        //stores the index of the read count column we are reading
        private int columnToReadIndex;
        private String sampleName;
        private TargetCollection<Target> targets;

        public ReadCountsTableReader(File file, final TargetCollection<Target> targets, TargetTableColumn columnToRead) throws IOException {
            super(file);
            this.columnToRead = columnToRead;
            this.targets = targets;
        }

        @Override
        protected boolean isCommentLine(final String[] line) {
            boolean isCommentLine;
            if(isCommentLine = (line.length > 0 && line[0].startsWith(TableUtils.COMMENT_PREFIX))) {
                parseCommentLine(line);
            }
            return isCommentLine;
        }

        private void parseCommentLine(final String[] line){
            //TODO extract sample name from a key value pair(i.e. sample = "samplename")
            this.sampleName = null;
        }

        @Override
        public void processColumns(final TableColumnCollection columns) {
            hasCoordinates = columns.containsAll(TargetTableColumn.CONTIG.toString(), TargetTableColumn.START.toString(),
                    TargetTableColumn.END.toString());
            hasName = columns.contains(TargetTableColumn.NAME.toString());
            if (!hasCoordinates && !hasName) {
                throw formatException("header contains neither coordinates nor target name columns");
            }
            if(!columns.contains(columnToRead.toString())) {
                throw formatException("header must contain both RAW and PCOV columns");
            }
            columnToReadIndex = columns.indexOf(columnToRead.toString());
        }

        @Override
        protected ReadCountRecord.SingleSampleRecord createRecord(final DataLine dataLine) {
            final Target target = createTarget(dataLine);
            final double counts = dataLine.getDouble(columnToReadIndex);
            return new ReadCountRecord.SingleSampleRecord(target, counts);
        }

        /**
         * Extracts the target object out of a data input line.
         * @param dataLine the input data line.
         * @return never {@code null}.
         */
        private Target createTarget(final DataLine dataLine) {
            if (hasName) {
                final String name = dataLine.get(TargetTableColumn.NAME);
                final Target target = targets.target(name);
                final SimpleInterval interval = createInterval(dataLine);
                if (target == null) {
                    return new Target(name, createInterval(dataLine));
                } else if (interval != null && !interval.equals(target.getInterval())) {
                    throw new UserException.BadInput(String.format("invalid target '%s' coordinates: expected %s but found %s",
                            name, target.getInterval(), createInterval(dataLine)));
                } else {
                    return target;
                }
            } else { // hasCoordinates must be true.
                final SimpleInterval interval = createInterval(dataLine);
                final Optional<Target> target = targets.targets(interval).stream().findAny();
                if (!target.isPresent() || !target.get().getInterval().equals(interval)) {
                    throw formatException("target not found with coordinates " + interval);
                }
                return target.get();
            }
        }

        /**
         * Extract the interval out of a data line.
         * @param dataLine the input data line.
         * @return {@code null} if the interval cannot be determined from the input file alone.
         */
        private SimpleInterval createInterval(final DataLine dataLine) {
            if (hasCoordinates) {
                return new SimpleInterval(dataLine.get(TargetTableColumn.CONTIG),
                        dataLine.getInt(TargetTableColumn.START),
                        dataLine.getInt(TargetTableColumn.END));
            } else {
                return null;
            }
        }
    }

    /**
     * Collection of read-count readers that simultaneously parse several input read-count files.
     * <p>
     * All the readers progress through their corresponding input files at the same pace; thus sharing the
     * current target.
     * </p>
     */
    private final class ReadCountReaderCollection implements AutoCloseable, Iterator<ReadCountRecord>, Iterable<ReadCountRecord> {
        private final List<ReadCountsTableReader> readers;

        //this string supposedly contains all the sample names
        private List<String> countColumnNames;

        private int[] countColumnSourceIndexMap;
        private final TargetCollection<Target> targets;
        private int targetsProcessedCount = 0;

        /**
         * Counts buffer used to accumulate counts coming from different readers.
         */
        private final double[] countsBuffer;

        @Override
        public Iterator<ReadCountRecord> iterator() { return this; }

        @Override
        public boolean hasNext() {
            return targetsProcessedCount < targets.targetCount();
        }

        @Override
        public ReadCountRecord next() {
            while (true) {
                final List<ReadCountRecord> nextReadCounts = readers.stream().map(r -> getNextRecord(r)).collect(Collectors.toList());
                final Target targetInFirstReader = nextReadCounts.get(0).getTarget();
                if (nextReadCounts.stream().map(ReadCountRecord::getTarget).anyMatch(t -> !t.equals(targetInFirstReader))) {
                    throw new UserException.BadInput(String.format("Target in file %s is %s but at least one input file" +
                            "has a different target at this position", readers.get(0).getSource(), targetInFirstReader));
                }
                if (targets.index(targetInFirstReader) != -1) {
                    targetsProcessedCount++;
                    return composeMergedReadCountRecord(nextReadCounts);
                }
            }
        }

        public ReadCountReaderCollection(final List<File> mergeGroup, final TargetCollection<Target> targets, TargetTableColumn columnToRead) {
            this.targets = targets;
            readers = mergeGroup.stream().map(f -> { try { return new ReadCountsTableReader(f, targets, columnToRead); } catch (final IOException ex) {
                                                        throw new UserException.CouldNotReadInputFile(f, ex);}}).collect(Collectors.toList());
            composeCountColumnNamesAndSourceIndexMapping();
            // pre-allocate count array used to accumulate the counts from all readers.
            countsBuffer = new double[countColumnNames.size()];
        }

        /**
         * Initializes count column name data-structures.
         * <p>
         * Initializes {@link #countColumnNames} and {@link #countColumnSourceIndexMap}
         * based on the input readers headers.
         * </p>
         * <p>
         * This operation must be performed after we have found the first target
         * in all sources; headers may not be yet defined before then.
         * </p>
         */
        private void composeCountColumnNamesAndSourceIndexMapping() {
            final List<String> unsortedCountColumnNames = new ArrayList<>();
            readers.stream().forEach(r -> unsortedCountColumnNames.add(r.sampleName));
            if (unsortedCountColumnNames.isEmpty()) {
                throw new IllegalStateException("there must be at least one count column");
            }
            countColumnSourceIndexMap = IntStream.range(0, unsortedCountColumnNames.size()).boxed()
                    .sorted(Comparator.comparing(unsortedCountColumnNames::get))
                    .mapToInt(Integer::intValue).toArray();
            countColumnNames = IntStream.of(countColumnSourceIndexMap)
                    .mapToObj(unsortedCountColumnNames::get)
                    .collect(Collectors.toList());
            checkForRepeatedSampleNames();
        }

        /**
         * Makes sure that there are no repeated sample names in the input.
         *
         * <p>
         *     At this point {@link #countColumnNames} is assumed to have all sample names (more than one) sorted.
         * </p>
         */
        private void checkForRepeatedSampleNames() {
            String previous = countColumnNames.get(0);
            for (int i = 1; i < countColumnNames.size(); i++) {
                final String next = countColumnNames.get(i);
                if (next.equals(previous)) {
                    throw new UserException.BadInput("the input contains the sample repeated, e.g.:" + next);
                }
                previous = next;
            }
        }

        @Override
        public void close() {
            for (final ReadCountsTableReader reader : readers) {
                try {
                    reader.close();
                } catch (final IOException ex) {
                    throw new GATKException(String.format("problems closing a read-count reader for %s", reader.getSource()), ex);
                }
            }
        }

        /**
         * Compose a merged ReadCountRecord from a List of records, assuming these records share the same target
         *
         * @return never {@code null}.
         */
        private ReadCountRecord composeMergedReadCountRecord(final List<ReadCountRecord> records) {
            int nextIndex = 0;
            for (final ReadCountRecord record : records) {
                final int size = record.size();
                record.copyCountsTo(countsBuffer, nextIndex);
                nextIndex += size;
            }
            final double[] counts = Arrays.stream(countColumnSourceIndexMap).mapToDouble(idx -> countsBuffer[idx]).toArray();
            return new ReadCountRecord(records.get(0).getTarget(), counts);
        }
    }

    private static ReadCountRecord getNextRecord(final ReadCountsTableReader reader) {
        try {
            final ReadCountRecord record = reader.readRecord();
            if (record == null) {
                throw new UserException.BadInput(String.format("End of file %s reached without finding all requested targets.", reader.getSource()));
            }
            return record;
        } catch (final IOException e) {
            throw new UserException.BadInput(String.format("End of file %s reached without finding all requested targets.", reader.getSource()));
        }
    }

    private final class ReadCountRecordWithSample extends ReadCountRecord.SingleSampleRecord {
        private final String sampleName;

        public ReadCountRecordWithSample(final Target target, final double counts, String sampleName) {
            super(target, counts);
            this.sampleName = sampleName;
        }

        public String getSampleName() { return this.sampleName; }
    }



}
