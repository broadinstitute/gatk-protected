package org.broadinstitute.hellbender.tools.exome;

import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SparkToggleCommandLineProgram;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hdf5.HDF5File;
import org.broadinstitute.hellbender.utils.hdf5.HDF5PoN;
import org.broadinstitute.hellbender.utils.hdf5.HDF5PoNCreator;
import org.broadinstitute.hellbender.utils.hdf5.PoN;
import org.broadinstitute.hellbender.utils.text.XReadLines;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * TODO SL
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Creates an Allelic Panel of Normals (APoN), given het pulldowns for samples that are part of the panel, " +
                "and appends it to the coverage PoN file created by CreatePanelOfNormals.  Samples can differ from those " +
                "used to create the coverage PoN but should be representative of the same sequencing process.  " +
                "Supports Apache Spark for some operations.",
        oneLineSummary = "Creates an Allelic Panel of Normals.",
        programGroup = CopyNumberProgramGroup.class
)
public final class CreateAllelicPanelOfNormals extends SparkToggleCommandLineProgram {

    private static final long serialVersionUID = 42123132L;

    @Argument(
            doc = "Input files for all samples in the panel of normals.  " +
                    "Each file should contain either a pulldown or a list of pulldown files; the latter should have a \".list\" extension.  " +
                    "Duplicate samples are not removed.",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional = false
    )
    protected List<File> inputFiles = new ArrayList<>();

    @Argument(
            doc = "HDF5 file created by CreatePanelOfNormals.  Allelic panel of normals will be appended to this file.",
            shortName = ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.PON_FILE_LONG_NAME,
            optional = false
    )
    protected File ponFile;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        validateArguments();

        final String originalLogLevel =
                (ctx.getLocalProperty("logLevel") != null) ? ctx.getLocalProperty("logLevel") : "INFO";
        ctx.setLogLevel("WARN");

        logger.info("Parsing input files...");
        final List<File> pulldownFiles = parseInputFiles(inputFiles);

        logger.info("Starting allelic panel of normals creation...");
        final AllelicPanelOfNormals allelicPoN = new AllelicPanelOfNormalsCreator(ctx, pulldownFiles).create();

        logger.info("Appending allelic panel of normals to HDF5 file...");
        final PoN coveragePoN = loadCoveragePoN(ponFile, logger);

        ctx.setLogLevel(originalLogLevel);
        logger.info("SUCCESS: Allelic panel of normals created and appended to " + ponFile + ".");
    }

    private void validateArguments() {
        Utils.regularReadableUserFile(ponFile);
        inputFiles.stream().forEach(Utils::regularReadableUserFile);
    }

    /**
     * Replaces any .list files in rawFileList with the files named in said .list file.
     * @param inputFiles the original file list, possibly including .list files
     * @return           a new List, with .list files replaced
     */
    private static List<File> parseInputFiles(final List<File> inputFiles) {
        final List<File> result = new ArrayList<>(inputFiles.size());
        for (final File rawFile : inputFiles) {
            if (rawFile.getName().endsWith(".list")) {
                try {
                    for (final String line : new XReadLines(rawFile, true))
                        result.add(new File(line));
                } catch (final IOException e) {
                    throw new UserException.CouldNotReadInputFile(rawFile, e);
                }
            } else {
                result.add(rawFile);
            }
        }
        return result;
    }

    private static PoN loadCoveragePoN(final File coveragePoNFile, final Logger logger) {
        Utils.regularReadableUserFile(coveragePoNFile);
        try (final HDF5File ponReader = new HDF5File(coveragePoNFile)) {
            final PoN pon = new HDF5PoN(ponReader);
            // Test the version of the PoN
            if (pon.getVersion() < HDF5PoNCreator.CURRENT_PON_VERSION) {
                logger.warn("The version of the specified PoN (" + pon.getVersion() + ") is older than the latest version " +
                        "(" + HDF5PoNCreator.CURRENT_PON_VERSION + ").");
            }
            return pon;
        } catch (final IllegalArgumentException e) {
            throw new UserException.BadInput("Could not load coverage panel of normals from " + coveragePoNFile + ".");
        }
    }
}
