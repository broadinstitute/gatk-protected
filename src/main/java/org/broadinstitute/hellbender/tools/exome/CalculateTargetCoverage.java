package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Calculates read-counts across targets for the exome copy number variant (CNV) calling workflow.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 *
 * <h3>Examples</h3>
 *
 * <p>
 *     The command encompasses empirically determined parameters for TCGA project data.
 *     You may obtain better results with different parameters.
 * </p>
 *
 * <p>For whole exome sequencing (WES) data: </p>
 *
 * <pre>
 * java -Xmx4g -jar $gatk_jar CalculateTargetCoverage \
 *   --input sample.bam \
 *   --targets padded_targets.tsv \
 *   --groupBy SAMPLE \
 *   --transform PCOV \
 *   --targetInformationColumns FULL \
 *   --disableReadFilter NotDuplicateReadFilter \
 *   --output sample.coverage.tsv
 * </pre>
 *
 * <p>
 *     The interval targets are exome target intervals padded, e.g. with 250 bases on either side.
 *     Target intervals do NOT overlap. Use the {@link PadTargets} tool to generate non-overlapping padded intervals from exome targets.
 *     Do NOT use BED format. See {@link ConvertBedToTargetFile}.
 * </p>
 *
 * <p>For whole genome sequencing (WGS) data, use {@link SparkGenomeReadCounts} instead.</p>
 *
 */
@CommandLineProgramProperties(
        summary = "Count overlapping reads target by target",
        oneLineSummary = "Count overlapping reads target by target",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class CalculateTargetCoverage extends ReadWalker {

    public static final String LINE_SEPARATOR = "\n";

    protected static final String TARGET_FILE_SHORT_NAME = ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME;
    protected static final String TARGET_FILE_FULL_NAME = ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME;


    @Argument(
            doc = "output tabular file with the counts",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    protected File output = null;

    @Argument(
            doc = "TSV file listing 1-based genomic intervals with specific column headers. Do NOT use BED format.",
            shortName = TARGET_FILE_SHORT_NAME,
            fullName = TARGET_FILE_FULL_NAME,
            optional = true
    )
    protected File targetsFile = null;

    /**
     * Default read count type
     */
    private ReadCountDataFactory.ReadCountType readCountType = ReadCountDataFactory.ReadCountType.RAW;

    /**
     * Writer to the main output file indicated by {@link #output}.
     */
    private PrintWriter outputWriter;

    /**
     * Reference to the logger.
     */
    private static final Logger logger = LogManager.getLogger(CalculateTargetCoverage.class);

    /**
     * Target database reference.
     */
    private TargetCollection<Target> targetCollection;

    /**
     * Map from targets to their corresponding read count data
     */
    private Map<Target, ReadCountData> readCountDataMap;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = new ArrayList<>(super.getDefaultReadFilters());
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);

        return filters;
    }

    @Override
    public void onTraversalStart() {
        if ( readArguments.getReadFilesNames().size() != 1 ) {
            throw new UserException.BadInput("This tool only accepts a single bam/sam/cram as input");
        }
        SampleCollection sampleCollection = new SampleCollection(getHeaderForReads());
        if (sampleCollection.sampleCount() > 1) {
            throw new UserException.BadInput("We do not support BAM files with multiple samples");
        }
        final String sampleName = sampleCollection.sampleIds().get(0);
        logger.log(Level.INFO, "Reading targets locations from intervals...");

        targetCollection = resolveTargetCollection();

        readCountDataMap = buildReadCountDataMap();

        // Open output files and write headers:
        outputWriter = openOutputWriter(output, composeMatrixOutputHeader(getCommandLine(), sampleName));

        // Next we start the traversal:
        logger.log(Level.INFO, "Collecting read counts ...");
    }

    /**
     * Builds the target collection given the values of user arguments.
     *
     * @throws UserException if there is some inconsistency in user arguments and inputs.
     * @throws GATKException if there was any problem parsing the content of the targets file.
     * @return never {@code null}.
     */
    private TargetCollection<Target> resolveTargetCollection() {
        final TargetCollection<Target> result;
        if (targetsFile != null) {
            result = resolveTargetsFromFile();
        } else if (hasIntervals()) {
            final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
            final List<SimpleInterval> intervals = intervalArgumentCollection.getIntervals(sequenceDictionary);

            //this constructor automatically generates target names from the intervals
            final List<Target> targets = intervals.stream().map(Target::new).collect(Collectors.toList());
            result = new HashedListTargetCollection<>(targets);
        } else {
            throw new UserException(String.format("You must indicate the set of target as input intervals (e.g. -L target-intervals.list) or a target feature file (e.g. -%s my-targets.tsv) ", TARGET_FILE_SHORT_NAME));
        }
        checkAllTargetsHaveName(result);
        return result;
    }

    /**
     * TODO
     * @return
     */
    private Map<Target,ReadCountData> buildReadCountDataMap() {
        final Map<Target, ReadCountData> result = new HashMap<>();
        targetCollection.targets().stream().forEach(target -> result.put(target, ReadCountDataFactory.getReadCountDataObject(readCountType, target)));
        return result;
    }

    /**
     * Checks whether all targets in the input target collection have a designated name.
     *
     * @param result the query target collection.
     * @throws UserException if there are some targets with no name in {@code result}.
     */
    private void checkAllTargetsHaveName(final TargetCollection<Target> result) {
        if (result.targets().stream().anyMatch(t -> t.getName() == null || t.getName().equals(""))) {
            throw new UserException(String.format("This tool requires that each target has a designated unique name/id but there are some with no names: %s",
                    result.targets().stream().filter(t -> t.getName() == null || t.getName().equals("")).limit(10).map(e -> result.location(e).toString()).collect(Collectors.joining(", "))));
        }
    }

    /**
     * Constructs the target collection from an target-file passed by the user.
     *
     * @return never {@code null}.
     */
    private TargetCollection<Target> resolveTargetsFromFile() {
        IOUtils.canReadFile(targetsFile);
        logger.log(Level.INFO,String.format("Reading target intervals from targets file '%s' ...", targetsFile.getAbsolutePath()));
        final List<Target> targets = TargetTableReader.readTargetFile(targetsFile);
        return new HashedListTargetCollection<>(targets);
    }

    /**
     * Opens the output file for writing with a print-writer.
     *
     * @param output     the output file.
     * @param headerText to be printed immediately after opening the writer.
     * @return never {@code null}.
     * @throws UserException.CouldNotCreateOutputFile if there was some problem creating or overwriting {@code output}.
     */
    private PrintWriter openOutputWriter(final File output, final String headerText) {
        try {
            final PrintWriter result = new PrintWriter(output);
            result.println(headerText);
            result.flush();
            return result;
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(output, e);
        }
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final SimpleInterval readLocation = referenceContext.getInterval();
        targetCollection.indexRange(readLocation).forEach(targetIndex -> readCountDataMap.get(targetCollection.target(targetIndex)).updateReadCount(read));
    }

    @Override
    public Object onTraversalSuccess() {
        logger.log(Level.INFO, "Collecting read counts done.");
        logger.log(Level.INFO, "Writing counts ...");

        List<String> readCountTableColumns = new ArrayList<>();
        readCountTableColumns.addAll(Arrays.asList(
                TargetTableColumn.CONTIG.toString(),
                TargetTableColumn.START.toString(),
                TargetTableColumn.END.toString(),
                TargetTableColumn.NAME.toString()));
        readCountTableColumns.addAll(ReadCountDataFactory.getColumnsOfReadCountType(readCountType));
        final TableColumnCollection columns = new TableColumnCollection(readCountTableColumns);

        try {
            TableWriter<ReadCountData> writer = getReadCountTableWriter(outputWriter, columns);
            for(Target target: targetCollection.targets()) {
                writer.writeRecord(readCountDataMap.get(target));
            }
        } catch(IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(output, "Could not create output file");
        }

        logger.log(Level.INFO, "Writing counts done.");

        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        if (outputWriter != null) {
            outputWriter.close();
        }
    }

    private TableWriter<ReadCountData> getReadCountTableWriter(final Writer writer, final TableColumnCollection tableColumns) throws IOException {
        return new TableWriter<ReadCountData>(writer, tableColumns) {
            @Override
            protected void composeLine(final ReadCountData record, final DataLine dataLine) {
                final SimpleInterval interval = record.getTarget().getInterval();
                if (interval == null) {
                    throw new IllegalStateException("invalid combination of targets with and without intervals defined");
                }
                dataLine.append(interval.getContig())
                        .append(interval.getStart())
                        .append(interval.getEnd())
                        .append(record.getTarget().getName());
                record.appendCountsTo(dataLine);
            }
        };
    }

    /**
     * Composes the main output header.
     *
     * @param commandLine      the tool command line.
     * @param sampleName        the name of the sample
     * @return never {@code null}.
     */
    private String composeMatrixOutputHeader(final String commandLine, final String sampleName) {
        final String formatString = String.join(LINE_SEPARATOR,
                "##fileFormat    = tsv",
                "##commandLine   = %s",
                "##title         = Read counts per target",
                ("##" + ReadCountFileHeaderKey.READ_COUNT_TYPE.getHeaderKeyName() + " = %s"),
                ("##" + ReadCountFileHeaderKey.SAMPLE_NAME.getHeaderKeyName() + "    = %s"));
        return String.format(formatString, commandLine, readCountType.toString(), sampleName);
    }

}
