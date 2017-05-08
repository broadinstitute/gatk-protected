package org.broadinstitute.hellbender.tools.pon;

import au.com.bytecode.opencsv.CSVWriter;
import com.opencsv.CSVReader;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.tools.exome.CreatePanelOfNormals;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetArgumentCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.pon.coverage.pca.HDF5PCACoveragePoN;
import org.broadinstitute.hellbender.tools.pon.coverage.pca.HDF5PCACoveragePoNCreationUtils;
import org.broadinstitute.hellbender.tools.pon.coverage.pca.PCACoveragePoN;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.testng.Assert;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.OptionalInt;

/**
 * Static class to help test PoN functionality
 */
public final class PoNTestUtils {
    private static final double DOUBLE_ARRAY_TOLERANCE = 1E-8;
    private static final double DOUBLE_MATRIX_TOLERANCE = 1E-8;

    private PoNTestUtils() {}

    /** Creates a HDF5 PoN (using {@link HDF5PCACoveragePoNCreationUtils} ).  Parameters use the
     * current {@link CreatePanelOfNormals} defaults.
     * @param inputPCovFile regular readable file that could be used to create a PoN.  Must be same format as output of
     *  {@link org.broadinstitute.hellbender.tools.exome.CombineReadCounts}
     * @param numEigensamples number of desired eigensamples in the PoN reduction
     * @return HDF5 File.  Never {@code null}
     */
    public static File createDummyHDF5FilePoN(final File inputPCovFile, final int numEigensamples) {
        IOUtils.canReadFile(inputPCovFile);
        ParamUtils.isPositive(numEigensamples, "Number of eigensamples must be greater than zero.");
        final File outputFile = IOUtils.createTempFile("dummy-pon-", ".pon");
        final TargetCollection<Target> targets = TargetArgumentCollection.readTargetCollection(inputPCovFile);
        HDF5PCACoveragePoNCreationUtils.create(null, outputFile, HDF5File.OpenMode.CREATE, inputPCovFile, targets, new ArrayList<>(), CreatePanelOfNormals.DEFAULT_TARGET_FACTOR_THRESHOLD_PERCENTILE, CreatePanelOfNormals.DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_COLUMN, CreatePanelOfNormals.DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_TARGET, CreatePanelOfNormals.DEFAULT_COLUMN_OUTLIER_DROP_THRESHOLD_PERCENTILE, CreatePanelOfNormals.DEFAULT_OUTLIER_TRUNCATE_PERCENTILE_THRESHOLD, OptionalInt.of(numEigensamples), false
        );
        return outputFile;
    }

    /**
     * Test whether two matrices are equal (within 1e-4)
     * @param left never {@code null}
     * @param right never {@code null}
     * @param isAllowNegatedValues whether values that are just negated are still considered equal.  True is useful for
     *                       the outputs of some matrix operations, such as SVD.
     */
    public static void assertEqualsMatrix(final RealMatrix left, final RealMatrix right, final boolean isAllowNegatedValues) {
        Assert.assertEquals(left.getRowDimension(), right.getRowDimension());
        Assert.assertEquals(left.getColumnDimension(), right.getColumnDimension());
        for (int i = 0; i < left.getRowDimension(); i++) {
            final double[] leftRow = left.getRow(i);
            final double[] rightRow = right.getRow(i);
            for (int j = 0; j < leftRow.length; j++) {
                if (isAllowNegatedValues) {
                    Assert.assertEquals(Math.abs(leftRow[j]), Math.abs(rightRow[j]), DOUBLE_MATRIX_TOLERANCE);
                } else {
                    Assert.assertEquals(leftRow[j], rightRow[j], DOUBLE_MATRIX_TOLERANCE);
                }
            }
        }
    }

    /**
     * Test whether two double arrays are equal.  For this method NaN is considered to equal NaN
     *
     * @param actual never {@code null}
     * @param gt never {@code null}
     */
    public static void assertEqualsDoubleArrays(final double[] actual, final double[] gt, final double tolerance) {
        Assert.assertEquals(actual.length, gt.length);
        for (int i = 0; i < actual.length; i++) {
            if (Double.isNaN(gt[i])) {
                Assert.assertTrue(Double.isNaN(actual[i]));
            } else {
                Assert.assertEquals(actual[i], gt[i], tolerance, String.format("Arrays were not equal (within tolerance %s) at index %d.", Double.toString(tolerance), i));
            }
        }
    }

    public static void assertEqualsDoubleArrays(final double[] actual, final double[] gt) {
        assertEqualsDoubleArrays(actual, gt, DOUBLE_ARRAY_TOLERANCE);
    }

    /**
     * Make sure that two PoNs are effectively the same.
     *
     * @param left never {@code null}
     * @param right never {@code null}
     */
    public static void assertEquivalentPoN(final File left, final File right) {
        IOUtils.canReadFile(left);
        IOUtils.canReadFile(right);
        try (final HDF5File leftFile = new HDF5File(left);
             final HDF5File rightFile = new HDF5File(right)) {
            final HDF5PCACoveragePoN leftPoN = new HDF5PCACoveragePoN(leftFile);
            final HDF5PCACoveragePoN rightPoN = new HDF5PCACoveragePoN(rightFile);
            assertEquivalentPoN(leftPoN, rightPoN);
        }
    }

    /**
     * Make sure that two PoNs are effectively the same.
     *
     * @param leftPoN never {@code null}
     * @param rightPoN never {@code null}
     */
    public static void assertEquivalentPoN(final PCACoveragePoN leftPoN, final PCACoveragePoN rightPoN) {
        Utils.nonNull(leftPoN, "Left PoN is null.");
        Utils.nonNull(rightPoN, "Right PoN is null.");

        Assert.assertEquals(leftPoN.getTargets(), rightPoN.getTargets());
        Assert.assertEquals(leftPoN.getRawTargets(), rightPoN.getRawTargets());
        Assert.assertEquals(leftPoN.getPanelTargets(), rightPoN.getPanelTargets());

        PoNTestUtils.assertEqualsDoubleArrays(leftPoN.getTargetFactors(), rightPoN.getTargetFactors());
        PoNTestUtils.assertEqualsDoubleArrays(leftPoN.getTargetVariances(), rightPoN.getTargetVariances());

        PoNTestUtils.assertEqualsMatrix(leftPoN.getNormalizedCounts(), rightPoN.getNormalizedCounts(), false);
        PoNTestUtils.assertEqualsMatrix(leftPoN.getLogNormalizedCounts(), rightPoN.getLogNormalizedCounts(), false);
        PoNTestUtils.assertEqualsMatrix(leftPoN.getLogNormalizedPInverseCounts(), rightPoN.getLogNormalizedPInverseCounts(), false);
        PoNTestUtils.assertEqualsMatrix(leftPoN.getReducedPanelCounts(), rightPoN.getReducedPanelCounts(), true);
        PoNTestUtils.assertEqualsMatrix(leftPoN.getReducedPanelPInverseCounts(), rightPoN.getReducedPanelPInverseCounts(), true);

        Assert.assertEquals(leftPoN.getTargetNames(), rightPoN.getTargetNames());
        Assert.assertEquals(leftPoN.getRawTargetNames(), rightPoN.getRawTargetNames());
        Assert.assertEquals(leftPoN.getPanelTargetNames(), rightPoN.getPanelTargetNames());

        Assert.assertEquals(leftPoN.getSampleNames(), rightPoN.getSampleNames());
        Assert.assertEquals(leftPoN.getPanelSampleNames(), rightPoN.getPanelSampleNames());
    }

    /**
     * Reads a very basic tsv (numbers separated by tabs) into a RealMatrix.
     * <p>Very little error checking happens in this method</p>
     *
     * @param inputFile readable file.  Not {@code null}
     * @return never {@code null}
     */
    public static RealMatrix readTsvIntoMatrix(final File inputFile) {
        IOUtils.canReadFile(inputFile);
        final List<double []> allData = new ArrayList<>();
        int ctr = 0;
        try {

            final CSVReader reader = new CSVReader(new FileReader(inputFile), '\t', CSVWriter.NO_QUOTE_CHARACTER);
            String[] nextLine;
            while ((nextLine = reader.readNext()) != null) {
                ctr++;
                allData.add(Arrays.stream(nextLine).filter(s -> StringUtils.trim(s).length() > 0).map(s -> Double.parseDouble(StringUtils.trim(s))).mapToDouble(d -> d).toArray());
            }
        } catch (final IOException ioe) {
            Assert.fail("Could not open test file: " + inputFile, ioe);
        }
        final RealMatrix result = new Array2DRowRealMatrix(allData.size(), allData.get(0).length);
        for (int i = 0; i < result.getRowDimension(); i++) {
            result.setRow(i, allData.get(i));
        }
        return result;
    }
}
