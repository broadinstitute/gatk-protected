package org.broadinstitute.hellbender.tools.exome.allelefraction;

import htsjdk.samtools.util.Log;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Tests for {@link AllelicPanelOfNormalsCreator}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicPanelOfNormalsCreatorUnitTest extends BaseTest {
    private static final String TEST_DIR = "/home/slee/working/ipython/alleliccapseg/allelic-pon/wgs-pon-liberal-MQ60-BEAGLEmaf10/pulldown";
    private static final List<File> PULLDOWN_FILES = Arrays.asList(
            new File(TEST_DIR, "TCGA-02-2483-10A-01D-1494-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-02-2485-10A-01D-1494-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-05-4396-10A-01D-1855-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-05-5429-10A-01D-1625-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-06-0157-10A-01D-1491-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-06-0214-10A-01D-1491-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-06-0686-10A-01D-1492-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-06-0744-10A-01D-1492-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-06-0745-10A-01D-1492-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-06-2557-10A-01D-1494-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-06-2570-10A-01D-1495-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-06-5415-10A-01D-1486-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-14-1823-10A-01D-1494-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-14-2554-10A-01D-1494-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-19-2620-10A-01D-1495-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-19-2624-10A-01D-1495-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-26-5132-10A-01D-1486-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-26-5135-10A-01D-1486-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-27-2523-10A-01D-1494-08.hets.tsv"),
            new File(TEST_DIR, "TCGA-27-2528-10A-01D-1494-08.hets.tsv"));

    @Test
    public void testCreate() {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final AllelicPanelOfNormalsCreator allelicPoNCreator = new AllelicPanelOfNormalsCreator(ctx, PULLDOWN_FILES);
        final double siteFrequency = 0.2;
        final AllelicPanelOfNormals allelicPoN = allelicPoNCreator.create(siteFrequency);
    }
}