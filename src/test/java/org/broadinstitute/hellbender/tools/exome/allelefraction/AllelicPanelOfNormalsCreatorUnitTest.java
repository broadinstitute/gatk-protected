package org.broadinstitute.hellbender.tools.exome.allelefraction;

import htsjdk.samtools.util.Log;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.exome.CopyNumberTriStateAllele;
import org.broadinstitute.hellbender.tools.exome.Genome;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Tests for {@link AllelicPanelOfNormalsCreator}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicPanelOfNormalsCreatorUnitTest extends BaseTest {
//    private static final String TEST_DIR = "/home/slee/working/ipython/alleliccapseg/allelic-pon/wgs-pon-liberal-MQ60-BEAGLEmaf10/pulldown";
//    private static final String TEST_DIR = "/home/slee/working/ipython/alleliccapseg/allelic-pon/wgs-chip/pulldown";
//    private static final String TEST_DIR = "/home/slee/working/ipython/alleliccapseg/allelic-pon/thca-pon-liberal-MQ60-BEAGLEmaf10/pulldown";
//    private static final FileFilter tsvFileFilter = new WildcardFileFilter("*tsv");
//    private static final List<File> PULLDOWN_FILES = Arrays.asList(new File(TEST_DIR).listFiles(tsvFileFilter)).subList(0, 100);
//    private static final List<File> PULLDOWN_FILES = Arrays.asList(new File(TEST_DIR).listFiles(tsvFileFilter));


//    private static final FileFilter normalFileFilter = new WildcardFileFilter("*10A*tsv");
//    private static final List<File> PULLDOWN_FILES = Arrays.asList(new File(TEST_DIR).listFiles(normalFileFilter));
//    private static final FileFilter tumorFileFilter = new WildcardFileFilter("*01A*tsv");
//    private static final List<File> TUMOR_PULLDOWN_FILES = Arrays.asList(new File(TEST_DIR).listFiles(tumorFileFilter));

//    private static final String TEST_DIR = "/home/slee/working/ipython/alleliccapseg/allelic-pon/luad-alpha-1.0-BEAGLE/pulldown";
//    private static final FileFilter tsvFileFilter = new WildcardFileFilter("*tsv");
//    private static final List<File> PULLDOWN_FILES = Arrays.asList(new File(TEST_DIR).listFiles(tsvFileFilter));
//    private static final FileFilter tumorFileFilter = new WildcardFileFilter("*01A*tsv");
//    private static final List<File> TUMOR_PULLDOWN_FILES = Arrays.asList(new File(TEST_DIR).listFiles(tumorFileFilter));

    private static final String TEST_DIR = "/home/slee/working/ipython/alleliccapseg/allelic-pon/prad-bayes-s30/pulldown";
    private static final FileFilter tsvFileFilter = new WildcardFileFilter("*tsv");
    private static final List<File> PULLDOWN_FILES = Arrays.asList(new File(TEST_DIR).listFiles(tsvFileFilter));
    private static final FileFilter tumorFileFilter = new WildcardFileFilter("*01A*tsv");
    private static final List<File> TUMOR_PULLDOWN_FILES = Arrays.asList(new File(TEST_DIR).listFiles(tumorFileFilter));

    @Test
    public void testCreate() {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final AllelicPanelOfNormalsCreator allelicPoNCreator = new AllelicPanelOfNormalsCreator(PULLDOWN_FILES);
        final double siteFrequency = 0.1;
        final AllelicPanelOfNormals allelicPoN = allelicPoNCreator.create(siteFrequency);
        allelicPoN.write(new File(TEST_DIR, "test.pon"));
    }

    @Test
    public void testMAFPoNNormalize() throws IOException {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final AllelicPanelOfNormals allelicPoN = new AllelicPanelOfNormals(new File(TEST_DIR, "test.pon"));

        for (final File file : ListUtils.union(PULLDOWN_FILES, TUMOR_PULLDOWN_FILES)) {
            final String sampleName = file.getName().split("[.]")[0];
            final AllelicCountCollection sampleCounts = new AllelicCountCollection(file);
            System.out.println(sampleName);

            final Genome genome = new Genome(AlleleFractionSimulatedData.TRIVIAL_TARGETS, sampleCounts.getCounts(), "sample");
            final List<SimpleInterval> sites = sampleCounts.getCounts().stream().map(AllelicCount::getInterval).collect(Collectors.toList());
            final SegmentedGenome segmentedGenome = new SegmentedGenome(sites, genome);
            final AlleleFractionData data = new AlleleFractionData(segmentedGenome);
            final AlleleFractionData dataWithPoN = new AlleleFractionData(segmentedGenome, allelicPoN);

            final AlleleFractionInitializer initializerWithPon = new AlleleFractionInitializer(dataWithPoN);
            final AlleleFractionState stateWithPoN = initializerWithPon.getInitializedState();
            writeMAFs(new File(TEST_DIR, sampleName + ".with.pon.maf"), sites, stateWithPoN);

            final AlleleFractionInitializer initializer = new AlleleFractionInitializer(data);
            final AlleleFractionState state = initializer.getInitializedState();
            writeMAFs(new File(TEST_DIR, sampleName + ".maf"), sites, state);
        }
        System.out.println("woo");
    }

    @Test
    public void test() {
        CopyNumberTriStateAllele.ALL_ALLELES.add(CopyNumberTriStateAllele.DEL);
    }

    private void writeMAFs(final File outputFile, final List<SimpleInterval> sites, final AlleleFractionState state) throws IOException {
        final FileWriter writer = new FileWriter(outputFile);
        final List<Double> mafs = IntStream.range(0, sites.size()).boxed().map(state::segmentMinorFraction).collect(Collectors.toList());
        for (final double maf : mafs) {
            writer.write(maf + "\n");
        }
        writer.close();
    }
}