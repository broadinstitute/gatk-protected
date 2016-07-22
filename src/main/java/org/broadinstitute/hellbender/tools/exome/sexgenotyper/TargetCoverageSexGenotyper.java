package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 * A command line tool for inferring sex genotypes from a tab-separated raw target coverage file.
 *
 * <p>
 *     In addition to the raw target coverage file, the user must provide a tab-separated contig
 *     ploidy annotation file. For homo sapiens, the annotation file may be as follows:
 *
 *     <pre>
 *         CONTIG_NAME    CONTIG_CLASS    SEX_XX    SEX_XY
 *         1              AUTOSOMAL       2          2
 *         2              AUTOSOMAL       2          2
 *                                   ...
 *                                   ...
 *                                   ...
 *         X              ALLOSOMAL       2          0
 *         Y              ALLOSOMAL       1          1
 *     </pre>
 *
 *     CONTIG_NAME column values must be same as contig names that appear in the read counts table.
 *     CONTIG_CLASS is either AUTOSOMAL or ALLOSOMAL. "SEX_XX" and "SEX_YY" are arbitrary sex genotype tags,
 *     along with their ploidies (= number of homologs). AUTOSOMAL contigs must have the same ploidy for all
 *     sexes. Every contig that appears in the read counts table must be annotated. One may include additional
 *     sex genotypes (along with contig ploidies) for specifies having more than two sexes by adding
 *     additional columns to the contig annotation file.
 * </p>
 *
 * <p>
 *     The provided target coverage file must contain both AUTOSOMAL and ALLOSOMAL targets. The tool
 *     infers the read depth density from AUTOSOMAL targets and uses this information to calculate the
 *     likelihood of sex genotypes.
 * </p>
 *
 * <p>
 *     Note: due to the uncertainty in the alignment of short reads, a small number of reads may be
 *     erroneously mapped to contigs with 0 actual ploidy (e.g. female homo sapiens samples may
 *     have a small number of reads aligned to the Y contig). The user must specifiy the typical mapping
 *     error probability for the tool to properly account for these errors.
 * </p>
 *
 * <p>
 *     Note: setting {@link TargetCoverageSexGenotyper#baselineMappingErrorProbability} to 0 (or to
 *     an unreasonably small number) will bias genotyping toward classes that cover more
 *     targets (SEX_XX < SEX_XY).
 * </p>
 *
 * <h2>Output File Format</h2>
 * <p>
 *     The inferred genotypes will be written to a tab-separated file such as:
 *
 *     <pre>
 *         SAMPLE_NAME                SEX_GENOTYPE      SEX_XX             SEX_XY
 *         arbitrary_XX_sample_name   SEX_XX            [log likelihood]   [log likelihood]
 *         arbitrary_XY_sample_name   SEX_XY            [log likelihood]   [log likelihood]
 *                                             ...
 *     </pre>
 * </p>
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Infers the sex genotypes from raw read counts.",
        oneLineSummary = "Infers the sex genotypes from raw read counts.",
        programGroup = CopyNumberProgramGroup.class
)

public class TargetCoverageSexGenotyper extends CommandLineProgram {

    private final Logger logger = LogManager.getLogger(TargetCoverageSexGenotyper.class);

    private static final String INPUT_READ_COUNT_COLLECTION_LONG_NAME = StandardArgumentDefinitions.INPUT_LONG_NAME;
    private static final String INPUT_READ_COUNT_COLLECTION_SHORT_NAME = StandardArgumentDefinitions.INPUT_SHORT_NAME;

    private static final String OUTPUT_SEX_GENOTYPE_LONG_NAME = StandardArgumentDefinitions.OUTPUT_LONG_NAME;
    private static final String OUTPUT_SEX_GENOTYPE_SHORT_NAME = StandardArgumentDefinitions.OUTPUT_SHORT_NAME;

    private static final String INPUT_CONTIG_ANNOTS_LONG_NAME = "contigAnnotations";
    private static final String INPUT_CONTIG_ANNOTS_SHORT_NAME = "annots";

    private static final String BASELINE_MAPPING_ERROR_PROBABILITY_LONG_NAME = "baselineMappingError";
    private static final String BASELINE_MAPPING_ERROR_PROBABILITY_SHORT_NAME = "mapErr";

    @Argument(
            doc = "Input raw read count collection tab-separated file.",
            fullName = INPUT_READ_COUNT_COLLECTION_LONG_NAME,
            shortName = INPUT_READ_COUNT_COLLECTION_SHORT_NAME,
            optional = false
    )
    protected File inputRawReadCountsFile;

    @Argument(
            doc = "Output sample sex genotype tab-separated file.",
            fullName = OUTPUT_SEX_GENOTYPE_LONG_NAME,
            shortName = OUTPUT_SEX_GENOTYPE_SHORT_NAME,
            optional = false
    )
    protected File outputSampleGenotypesFile;

    @Argument(
            doc = "Input contig annotations file.",
            fullName = INPUT_CONTIG_ANNOTS_LONG_NAME,
            shortName = INPUT_CONTIG_ANNOTS_SHORT_NAME,
            optional = false
    )
    protected File inputContigAnnotsFile;

    @Argument(
            doc = "Baseline mapping error probability.",
            fullName = BASELINE_MAPPING_ERROR_PROBABILITY_LONG_NAME,
            shortName = BASELINE_MAPPING_ERROR_PROBABILITY_SHORT_NAME,
            optional = true
    )
    protected double baselineMappingErrorProbability = 1e-4;

    @Override
    protected Object doWork() {
        /* check args */
        Utils.regularReadableUserFile(inputRawReadCountsFile);
        Utils.regularReadableUserFile(inputContigAnnotsFile);
        ParamUtils.inRange(baselineMappingErrorProbability, 0, 1, "Baseline mapping error probability must be a positive " +
                "real number between 0 and 1");

        /* read input data */
        ReadCountCollection rawReadCounts;
        List<ContigPloidyAnnotation> contigPloidyAnnotationList;
        try {
            logger.info("Parsing raw read count collection file...");
            rawReadCounts = ReadCountCollectionUtils.parse(inputRawReadCountsFile);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read raw read count collection file");
        }
        try {
            logger.info("Parsing contig genotype ploidy annotations file...");
            contigPloidyAnnotationList = ContigPloidyAnnotationTableReader.readContigPloidyAnnotationsFromFile(inputContigAnnotsFile);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read contig genotype ploidy annotations file");
        }

        /* perform genotyping */
        final TargetCoverageSexGenotypeCalculator genotyper = new TargetCoverageSexGenotypeCalculator(rawReadCounts,
                contigPloidyAnnotationList, baselineMappingErrorProbability);
        final SexGenotypeDataCollection sampleSexGenotypeCollection = genotyper.inferSexGenotypes();

        /* save results */
        try {
            sampleSexGenotypeCollection.write(new FileWriter(outputSampleGenotypesFile));
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile("Could not write inferred sample genotypes to file", ex);
        }

        return "SUCCESS";
    }
}
