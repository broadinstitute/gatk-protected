package org.broadinstitute.hellbender.tools.genome;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.utils.gc.GcCorrector;

import java.io.File;
import java.io.IOException;

public final class CorrectGc extends CommandLineProgram{
    public static final String READ_COUNTS_FILE_FULL_NAME = StandardArgumentDefinitions.INPUT_LONG_NAME;
    public static final String READ_COUNTS_FILE_SHORT_NAME = StandardArgumentDefinitions.INPUT_SHORT_NAME;

    public static final String REFERENCE_FILE_FULL_NAME = StandardArgumentDefinitions.REFERENCE_LONG_NAME;
    public static final String REFERENCE_FILE_SHORT_NAME = StandardArgumentDefinitions.REFERENCE_SHORT_NAME;

    @Argument(
            doc = "read counts input file.  This can only contain one sample at a time.",
            shortName = READ_COUNTS_FILE_SHORT_NAME,
            fullName = READ_COUNTS_FILE_FULL_NAME,
            optional = false
    )
    protected File readCountsFile;

    @Argument(
            doc = "reference file in the fasta format.",
            shortName = REFERENCE_FILE_SHORT_NAME,
            fullName = REFERENCE_FILE_FULL_NAME,
            optional = false
    )
    protected File referenceFile;

    @Override
    protected Object doWork() {
        try {
            final ReadCountCollection readCountCollection = ReadCountCollectionUtils.parse(readCountsFile);
            final ReadCountCollection gcCorrectedreadCountCollection = GcCorrector.correctGc(readCountCollection, referenceFile);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(readCountsFile, ex.getMessage(), ex);
        }
        return null;
    }
}
