package org.broadinstitute.hellbender.tools.exome.allelefraction;

import htsjdk.samtools.util.Log;
import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.exome.allelefraction.pon.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.tools.exome.allelefraction.pon.AllelicPanelOfNormalsCreator;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileFilter;
import java.util.Arrays;
import java.util.List;

/**
 * Tests for {@link AllelicPanelOfNormalsCreator}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicPanelOfNormalsCreatorUnitTest extends BaseTest {
    private static final String TEST_DIR = "/home/slee/working/ipython/alleliccapseg/allelic-pon/wgs-pon-chip/pulldown";
    private static final FileFilter tsvFileFilter = new WildcardFileFilter("*tsv");
    private static final List<File> PULLDOWN_FILES = Arrays.asList(new File(TEST_DIR).listFiles(tsvFileFilter));

    @Test
    public void testCreate() {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final AllelicPanelOfNormalsCreator allelicPoNCreator = new AllelicPanelOfNormalsCreator(ctx, PULLDOWN_FILES);
        final double siteFrequency = 0.2;
        final AllelicPanelOfNormals allelicPoN = allelicPoNCreator.create(siteFrequency);
    }
}