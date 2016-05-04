package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.broadinstitute.hellbender.tools.exome.AllelicCountTableColumns.AllelicCountTableVerbosity;

import java.io.File;
import java.io.IOException;

/**
 * Writes {@link AllelicCount} instances to a tab-separated table file.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class AllelicCountWriter extends TableWriter<AllelicCount> {

    private final AllelicCountTableVerbosity verbosity;

    public AllelicCountWriter(final File file, final AllelicCountTableVerbosity verbosity)
            throws IOException {
        super(file, new TableColumnCollection(AllelicCountTableColumns.getColumns(verbosity)));
        this.verbosity = verbosity;
    }

    @Override
    protected void composeLine(final AllelicCount record, final DataLine dataLine) {
        switch (verbosity) {
            case BASIC:
                composeLineBasic(record, dataLine);
                break;
            case INTERMEDIATE:
                composeLineIntermediate(record, dataLine);
                break;
            case FULL:
                composeLineFull(record, dataLine);
                break;
        }
    }

    /**
     * Compose the basic-verbosity record: (contig, position, ref read count, alt read count)
     *
     * @param record the {@link AllelicCount} record
     * @param dataLine the {@link DataLine} to the composed
     */
    private static void composeLineBasic(final AllelicCount record, final DataLine dataLine) {
        dataLine.append(record.getInterval().getContig())
                .append(record.getInterval().getEnd())
                .append(record.getRefReadCount())
                .append(record.getAltReadCount());
    }

    /**
     * Compose the intermediate-verbosity record: (contig, position, ref read count, alt read count,
     * ref nucleotide, alt nucleotide, read depth)
     *
     * @param record the {@link AllelicCount} record
     * @param dataLine the {@link DataLine} to the composed
     */
    private static void composeLineIntermediate(final AllelicCount record, final DataLine dataLine) {
        dataLine.append(record.getInterval().getContig())
                .append(record.getInterval().getEnd())
                .append(record.getRefReadCount())
                .append(record.getAltReadCount())
                .append(record.getRefNucleotide().name())
                .append(record.getAltNucleotide().name())
                .append(record.getReadDepth());
    }

    /**
     * Compose the full-verbosity record: (contig, position, ref read count, alt read count,
     * ref nucleotide, alt nucleotide, read depth, heterozygosity log odds)
     *
     * @param record the {@link AllelicCount} record
     * @param dataLine the {@link DataLine} to the composed
     */
    private static void composeLineFull(final AllelicCount record, final DataLine dataLine) {
        dataLine.append(record.getInterval().getContig())
                .append(record.getInterval().getEnd())
                .append(record.getRefReadCount())
                .append(record.getAltReadCount())
                .append(record.getRefNucleotide().name())
                .append(record.getAltNucleotide().name())
                .append(record.getReadDepth())
                .append(String.format("%.4f", record.getHetLogOdds()));
    }
}
