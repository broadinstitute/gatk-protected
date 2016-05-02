package org.broadinstitute.hellbender.tools.exome;

import com.google.common.collect.Sets;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

/**
 * Reads {@link AllelicCount} instances from a tab separated table file.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class AllelicCountReader extends TableReader<AllelicCount> {

    private final AllelicCountTableColumns.AllelicCountTableVerbosity verbosity;

    /**
     * Opens a reader on an a preexisting file.
     *
     * @param file the source file where to read from.
     *
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws IOException if there is an issue trying to read the contents of the file.
     * @throws RuntimeException if there is a formatting issue within the file.
     * @throws UserException.BadInput if not all mandatory columns can be found.
     */
    public AllelicCountReader(final File file) throws IOException {
        super(file);

        /* detect verbosity level */
        if (columns().containsAll(AllelicCountTableColumns.FULL_COLUMN_NAME_ARRAY)) {
            verbosity = AllelicCountTableColumns.AllelicCountTableVerbosity.FULL;
        } else if (columns().containsAll(AllelicCountTableColumns.INTERMEDIATE_COLUMN_NAME_ARRAY)) {
            verbosity = AllelicCountTableColumns.AllelicCountTableVerbosity.INTERMEDIATE;
        } else if (columns().containsAll(AllelicCountTableColumns.BASIC_COLUMN_NAME_ARRAY)) {
            verbosity = AllelicCountTableColumns.AllelicCountTableVerbosity.BASIC;
        } else {
            final Set<String> missingColumns = Sets.difference(
                    new HashSet<>(AllelicCountTableColumns.BASIC_COLUMN_NAME_ARRAY), new HashSet<>(columns().names()));
            throw new UserException.BadInput("Bad header in AllelicCount file. Not all mandatory columns are present." +
                    " Missing: " + StringUtils.join(missingColumns, ", "));
        }

    }

    @Override
    protected AllelicCount createRecord(DataLine dataLine) {

        /* mandatory (basic) fields */
        final int position = dataLine.getInt(AllelicCountTableColumns.POSITION.name());
        final SimpleInterval interval = new SimpleInterval(
                dataLine.get(AllelicCountTableColumns.CONTIG.name()), position, position);
        final int refReadCount = dataLine.getInt(AllelicCountTableColumns.REF_COUNT.name());
        final int altReadCount = dataLine.getInt(AllelicCountTableColumns.ALT_COUNT.name());

        /* extra fields */
        final Nucleotide refNucleotide, altNucleotide;
        final int readDepth;
        final double hetLogOdds;

        switch (verbosity) {

            case BASIC:
                return new AllelicCount(interval, refReadCount, altReadCount);

            case INTERMEDIATE:
                refNucleotide = Nucleotide.valueOf(dataLine.get(AllelicCountTableColumns.REF_NUCLEOTIDE.name()).getBytes()[0]);
                altNucleotide = Nucleotide.valueOf(dataLine.get(AllelicCountTableColumns.ALT_NUCLEOTIDE.name()).getBytes()[0]);
                readDepth = dataLine.getInt(AllelicCountTableColumns.READ_DEPTH.name());
                return new AllelicCount(interval, refReadCount, altReadCount, refNucleotide, altNucleotide, readDepth);

            case FULL:
                refNucleotide = Nucleotide.valueOf(dataLine.get(AllelicCountTableColumns.REF_NUCLEOTIDE.name()).getBytes()[0]);
                altNucleotide = Nucleotide.valueOf(dataLine.get(AllelicCountTableColumns.ALT_NUCLEOTIDE.name()).getBytes()[0]);
                readDepth = dataLine.getInt(AllelicCountTableColumns.READ_DEPTH.name());
                hetLogOdds = dataLine.getDouble(AllelicCountTableColumns.HET_LOG_ODDS.name());
                return new AllelicCount(interval, refReadCount, altReadCount, refNucleotide, altNucleotide, readDepth,
                        hetLogOdds);

            default:
                throw new RuntimeException("The AllelicCount table verbosity level is not properly set.");
        }
    }

}
