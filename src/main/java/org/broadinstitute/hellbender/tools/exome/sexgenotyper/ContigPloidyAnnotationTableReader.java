package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import com.google.common.collect.Sets;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import javax.annotation.Nonnull;
import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Reads contig ploidy annotations from tab-separated files and readers.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class ContigPloidyAnnotationTableReader extends TableReader<ContigPloidyAnnotation> {

    public static final Logger logger = LogManager.getLogger(ContigPloidyAnnotationTableReader.class);

    /**
     * The set of ploidy tags (= sex genotype identifer strings)
     */
    private final Set<String> ploidyTagsSet;

    /**
     * Public constructor.
     *
     * @param sourceName name of the source
     * @param sourceReader an instance of {@link Reader}
     * @throws IOException if a reading error occurs
     */
    public ContigPloidyAnnotationTableReader(final String sourceName, @Nonnull final Reader sourceReader)
            throws IOException {
        super(sourceName, sourceReader);
        TableUtils.checkMandatoryColumns(columns(), ContigPloidyAnnotationTableColumn.MANDATORY_CONTIG_ANNOTATION_COLUMNS,
                UserException.BadInput::new);
        ploidyTagsSet = Sets.difference(new HashSet<>(columns().names()),
                ContigPloidyAnnotationTableColumn.MANDATORY_CONTIG_ANNOTATION_COLUMNS_SET);
        if (ploidyTagsSet.isEmpty()) {
            throw new UserException.BadInput("At least one ploidy column is required!");
        }
        logger.info("Ploidy tags: " + ploidyTagsSet.stream().collect(Collectors.joining(", ")));
    }

    /**
     * Public constructor.
     *
     * @param sourceReader an instance of {@link Reader}
     * @throws IOException if a reading error occurs
     */
    public ContigPloidyAnnotationTableReader(@Nonnull final Reader sourceReader) throws IOException {
        this(null, sourceReader);
    }

    /**
     * Creates a {@link ContigPloidyAnnotation} instance from a {@link DataLine}
     *
     * @param dataLine a data line
     * @return an instance of {@link ContigPloidyAnnotation}
     */
    @Override
    protected ContigPloidyAnnotation createRecord(@Nonnull final DataLine dataLine) {
        final String contigName = dataLine.get(ContigPloidyAnnotationTableColumn.CONTIG_NAME);

        final ContigClass contigClass;
        final String contigClassString = dataLine.get(ContigPloidyAnnotationTableColumn.CONTIG_CLASS);
        if (!ContigClass.CONTIG_CLASS_NAMES_SET.contains(contigClassString)) {
            throw new UserException.BadInput("Bad contig class: provided value: " + contigClassString + ", acceptable values: " +
                    ContigClass.CONTIG_CLASS_NAMES_SET.stream().collect(Collectors.joining(", ", "[", "]")));
        } else {
            contigClass = ContigClass.valueOf(contigClassString);
        }

        final Map<String, Integer> ploidyMap = new HashMap<>();
        try {
            /* all lines must have all ploidy annotations defined in the header */
            ploidyTagsSet.forEach(tag -> ploidyMap.put(tag, dataLine.getInt(tag)));
        } catch (final IllegalArgumentException ex) {
            throw new UserException.BadInput("All lines in the contig ploidy annotation table must have values for all" +
                    " ploidy classes; " + ex.getMessage());
        }

        return new ContigPloidyAnnotation(contigName, contigClass, ploidyMap);
    }

    /**
     * Reads contig ploidy annotations from a file.
     *
     * @param contigPloidyAnnotationsFile a file containing contig ploidy annotations
     * @return a list of contig ploidy annotations
     * @throws IOException if a read error occurs
     */
    public static List<ContigPloidyAnnotation> readContigPloidyAnnotationsFromFile(@Nonnull final File contigPloidyAnnotationsFile)
            throws IOException {
        Utils.regularReadableUserFile(contigPloidyAnnotationsFile);
        try {
            return readContigPloidyAnnotationsFromReader(contigPloidyAnnotationsFile.getAbsolutePath(),
                    new FileReader(contigPloidyAnnotationsFile));
        } catch (final FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read contig ploidy annotations file " +
                    contigPloidyAnnotationsFile.getAbsolutePath());
        }
    }

    /**
     * Reads contig ploidy annotations from a {@link Reader}.
     *
     * @param contigAnnotationSourceName a string identifier for the reader.
     * @param contigAnnotationReader an instance of {@link Reader}.
     * @return list of contig ploidy annotations
     * @throws IOException if a read error occurs
     */
    public static List<ContigPloidyAnnotation> readContigPloidyAnnotationsFromReader(@Nonnull final String contigAnnotationSourceName,
                                                                                     @Nonnull final Reader contigAnnotationReader)
            throws IOException {
        /* read contig annotations */
        try (final ContigPloidyAnnotationTableReader reader =
                     new ContigPloidyAnnotationTableReader(contigAnnotationSourceName, contigAnnotationReader)) {
            return reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(contigAnnotationSourceName, e);
        }
    }
}
