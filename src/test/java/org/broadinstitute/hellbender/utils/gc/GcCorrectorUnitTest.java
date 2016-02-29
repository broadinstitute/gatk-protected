package org.broadinstitute.hellbender.utils.gc;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class GcCorrectorUnitTest extends BaseTest{
    @Test()
    public void testReadCountCollectionGcCorrection() {
        final File INPUT_FILE = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/full-read-counts.txt");
        final File FASTA = new File("src/test/resources/org/broadinstitute/hg19mini.fasta");
        try {
            final ReadCountCollection readCountCollection = ReadCountCollectionUtils.parse(INPUT_FILE);
            final ReadCountCollection gcCorrectedreadCountCollection = GcCorrector.correctGc(readCountCollection, FASTA);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(INPUT_FILE, ex.getMessage(), ex);
        }
    }
}
