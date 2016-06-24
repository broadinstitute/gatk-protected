package org.broadinstitute.hellbender.tools.exome.allelefraction;

import com.google.common.collect.Sets;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountTableColumn;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;
import java.util.AbstractMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Reads an {@link AllelicPanelOfNormals} from a tab-separated table file.
 *
 * @author Samuel Lee &lt;mehrtash@broadinstitute.org&gt;
 */
public class AllelicPanelOfNormalsReader extends TableReader<Map.Entry<SimpleInterval, AllelicPanelOfNormals.HyperparameterValues>> {
    /**
     * Opens a reader on an a pre-existing allelic panel of normals file.
     *
     * @param file the source file where to read from.
     *
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws IOException if there is an issue trying to read the contents of the file.
     * @throws RuntimeException if there is a formatting issue within the file.
     */
    public AllelicPanelOfNormalsReader(final File file) throws IOException {
        super(file); /* the constructor of TableReader parses the header */
    }

    @Override
    protected Map.Entry<SimpleInterval, AllelicPanelOfNormals.HyperparameterValues> createRecord(final DataLine dataLine) {
        final int position = dataLine.getInt(AllelicPanelOfNormalsTableColumn.POSITION);
        final SimpleInterval interval = new SimpleInterval(dataLine.get(AllelicPanelOfNormalsTableColumn.CONTIG), position, position);
        final double alpha = dataLine.getDouble(AllelicPanelOfNormalsTableColumn.ALPHA);
        final double beta = dataLine.getDouble(AllelicPanelOfNormalsTableColumn.BETA);
        final AllelicPanelOfNormals.HyperparameterValues hyperparameterValues = new AllelicPanelOfNormals.HyperparameterValues(alpha, beta);

        return new AbstractMap.SimpleEntry<>(interval, hyperparameterValues);
    }
}
