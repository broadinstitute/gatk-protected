package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ReadCountFileHeaderKey;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.*;

/**
 * Integration tests for {@link CalculateTargetCoverage}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class CalculateTargetCoverageIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/calculatetargetcoverage/temp/");
    private static File testFile(final String fileName) {
        return new File(TEST_DIR, fileName);
    }

    //////////////////////
    // Test input files //
    //////////////////////

    private final static File INTERVALS_LIST = testFile("exome-read-counts-intervals.list");

    private final static File INTERVALS_BED = testFile("exome-read-counts-intervals.tsv");

    private final static File NA12878_BAM = testFile("exome-read-counts-NA12878.bam");

    private final static File NA12778_BAM = testFile("exome-read-counts-NA12778.bam");

    private final static File NA12872_BAM = testFile("exome-read-counts-NA12872.bam");

    private final static File NA12878_NA12778_BAM_MERGED = testFile("exome-read-counts-NA12878-NA12778-merged.bam");

    private final static File INTERVALS_BED_MISSING_NAMES = testFile("exome-read-counts-intervals-missing-names.tsv");

    private final static File INTERVALS_LIST_DUPS = testFile("exome-read-counts-intervals_dups.list");

    private final static File INTERVALS_BED_DUPS = testFile("exome-read-counts-intervals_dups.tsv");

    ////////////////////////
    //  Test output files //
    ////////////////////////

    private final static File NA12878_RAW_COUNT_EXPECTED_OUTPUT_INTERVALS = testFile("NA12878_raw_counts_intervals.tsv");

    @Override
    public String getTestedClassName() {
        return CalculateTargetCoverage.class.getSimpleName();
    }

    /**
     * Creates a temp file making sure it is going to be removed after testing VM finishes.
     * @param id temporal file name identifiable part.
     * @return never {@code null}.
     */
    private File createTempFile(final String id) {
        return BaseTest.createTempFile("GATK4-ERCTest-" + id,".tmp");
    }

    @DataProvider(name="correctRunData")
    public Object[][] correctRunData() {
        return new Object[][] {
                {
                        new File[]{ NA12878_BAM },
                        INTERVALS_LIST,
                        NA12878_RAW_COUNT_EXPECTED_OUTPUT_INTERVALS,
                        new String[0]
                }

        };
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testMultipleBamsInput() {
        testCorrectRun(
                new File[]{ NA12878_BAM, NA12778_BAM, NA12872_BAM },
                INTERVALS_LIST,
                NA12878_RAW_COUNT_EXPECTED_OUTPUT_INTERVALS,
                new String[0]
        );
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testMultipleSamplesBamInput() {
        testCorrectRun(
                new File[]{ NA12878_NA12778_BAM_MERGED },
                INTERVALS_LIST,
                NA12878_RAW_COUNT_EXPECTED_OUTPUT_INTERVALS,
                new String[0]
        );
    }



    @Test(dataProvider = "correctRunData")
    public void testCorrectRun(final File[] bamFiles, final File intervalFile, final File expectedOutputFile, final String[] additionalArguments) {
        final File outputFile = createTempFile("cohort-output");
        final List<String> headerKeysToCompare = Arrays.asList(ReadCountFileHeaderKey.READ_COUNT_TYPE.getHeaderKeyName(),
                ReadCountFileHeaderKey.SAMPLE_NAME.getHeaderKeyName());
        final List<String> arguments = new ArrayList<>(Arrays.asList(
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        ));
        if (intervalFile != null) {
            arguments.add("-L");
            arguments.add(intervalFile.getAbsolutePath());
        }
        Arrays.asList(bamFiles).forEach(bam -> {
            arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
            arguments.add(bam.getAbsolutePath());
        });
        arguments.addAll(Arrays.asList(additionalArguments));
        System.err.println("COMMAND LINE: " + Arrays.toString(arguments.toArray()));
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists(), "output file does not exist: " + outputFile);
        Assert.assertTrue(outputFile.isFile(), "output file is not a regular file: " + outputFile);
        Assert.assertTrue(outputFile.canRead(), "output file cannot be read: " + outputFile);
        compareTableFiles(outputFile, expectedOutputFile, headerKeysToCompare, "matrix output: " + outputFile + " " + expectedOutputFile);
    }

    private void compareTableFiles(File outputFile, File expectedOutputFile, List<String> headerKeysToCompare, final String role) {
        try {
            final Table actualTable = Table.fromFile(outputFile);
            final Table expectedTable = Table.fromFile(expectedOutputFile);
            Table.assertEquals(actualTable, expectedTable, outputFile, expectedOutputFile, headerKeysToCompare);
        } catch (final IOException ex) {
            Assert.fail("Failed to read actual or expected " +  role + " table files ", ex);
        }
    }

    public static class Table {

        private String[] values;
        private final int columnCount;
        private final int rowCount;
        private String[] columnNames;
        private Map<String, String> headerValuesMap;

        public static Table fromFile(final File file) throws IOException {
            try (final FileReader fr = new FileReader(file)) {
                return new Table(fr);
            }
        }

        public Table(final Reader reader) throws IOException {
            Utils.nonNull(reader, "the reader cannot be null");
            final BufferedReader lineReader = new BufferedReader(reader);
            headerValuesMap = new HashMap<>(10);
            String lastLine;
            String headerKeyValuePattern = "^[# ]*(\\S*)[\\s\t]*=[\\s\t]*(.*\\S)[\\s\t]*$";
            while ((lastLine = lineReader.readLine()) != null) {
                if (lastLine.matches("^#.*$")) {
                    headerValuesMap.put(lastLine.replaceAll(headerKeyValuePattern, "$1"), lastLine.replaceAll(headerKeyValuePattern, "$2"));
                } else {
                    break;
                }
            }
            final String[] header = lastLine == null ? null : lastLine.split("\\t");
            Utils.nonNull(header, "the table does not have a header");
            columnCount = header.length;

            final List<String> values = new ArrayList<>(100);
            int lineNumber = 1 + 1 + headerValuesMap.size();
            while ((lastLine = lineReader.readLine()) != null) {
                final String[] lineValues = lastLine.split("\\t");
                if (lineValues.length != columnCount) {
                    throw new IllegalArgumentException("line " + lineNumber + " does not has the expected number of columns " + columnCount + ": " + lastLine);
                }
                values.addAll(Arrays.asList(lineValues));
            }

            this.values = values.toArray(new String[values.size()]);
            rowCount = values.size() / columnCount;
            columnNames = header;
        }

        public static void assertEquals(final Table tb1, final Table tb2, final File f1, final File f2, List<String> headerKeysToCompare) {

            //check that mandatory header key value pairs are the same
            headerKeysToCompare.stream().forEach(key -> {
                Assert.assertTrue(tb1.headerValuesMap.containsKey(key));
                Assert.assertTrue(tb2.headerValuesMap.containsKey(key));
                Assert.assertEquals(tb1.headerValuesMap.get(key), tb2.headerValuesMap.get(key));
                });

            Assert.assertEquals(tb1.columnCount,tb2.columnCount,"Different number of columns");
            Assert.assertEquals(tb1.rowCount,tb2.rowCount,"Different number of rows");
            Assert.assertEquals(tb1.columnNames,tb2.columnNames,"Differences in header");
            assertValuesAreEqual(tb1.values, tb2.values, 0.01, "Difference in values between " + f1 + " and " + f2);

        }

        private static void assertValuesAreEqual(final String[] v1, final String[] v2, double epsilon, String message) {
            Assert.assertEquals(v1.length, v2.length, message);
            for (int i = 0; i < v1.length; i++) {
                final boolean v1IsDouble = isADouble(v1[i]);
                final boolean v2IsDouble = isADouble(v2[i]);
                Assert.assertEquals(v1IsDouble,v2IsDouble,message);
                if (v1IsDouble) {
                    if (Double.isNaN(Double.valueOf(v1[i]))) {
                        Assert.assertTrue(Double.isNaN(Double.valueOf(v2[i])));
                    } else {
                        Assert.assertEquals(Double.valueOf(v1[i]), Double.valueOf(v2[i]), epsilon, message);
                    }
                } else {
                    Assert.assertEquals(v1[i],v2[i],message);
                }
            }
        }

        private static boolean isADouble(final String s) {
            if (s == null) {
                return false;
            } else {
                try {
                    Double.parseDouble(s);
                } catch (final NumberFormatException n) {
                    return false;
                }
                return true;
            }
        }

        public String getString(int r, int c) {
            return values[r * columnCount + c];
        }

        public long getLong(int r, int c) {
            return Long.parseLong(values[r * columnCount + c]);
        }

        public double getDouble(int r, int c) {
            return Double.parseDouble(values[r * columnCount + c]);
        }
    }
}
