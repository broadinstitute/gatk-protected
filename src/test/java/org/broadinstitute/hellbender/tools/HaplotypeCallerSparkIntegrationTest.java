package org.broadinstitute.hellbender.tools;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.apache.spark.SparkException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.SparkTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

import static org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerIntegrationTest.TEST_FILES_DIR;
import static org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerIntegrationTest.calculateConcordance;

public class HaplotypeCallerSparkIntegrationTest extends CommandLineProgramTest {

    /*
    * Test that in VCF mode we're >= 99% concordant with GATK3.5 results
    */
    @Test
    public void testVCFModeIsConcordantWithGATK3_5Results() throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConcordantWithGATK3Results", ".vcf");
        //Created by running:
        //java -jar ~/bin/GenomeAnalysisTK-3.5.0/GenomeAnalysisTK.jar -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // --out a.gatk3.5.noDownsample.vcf -G StandardHC -G Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testVCFMode.gatk3.5.vcf");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_2bit_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.5 in VCF mode is < 99%");
    }

    /**
     * Test that in VCF mode we're >= 99% concordant with GATK3.5 results
     * THIS TEST explodes with an exeption because Allele-Specific annotations are not supported in vcf mode yet.
     * It's included to parallel the matching (also exploding) test for the non-spark HaplotypeCaller
     * {@link org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerIntegrationTest#testVCFModeIsConcordantWithGATK3_5ResultsAlleleSpecificAnnotations()}
     */
    @Test(expectedExceptions = SparkException.class) //this should be a UserException, but spark exceptions are not unwrapped yet
    public void testVCFModeIsConcordantWithGATK3_5ResultsAlleleSpecificAnnotations() throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConcordantWithGATK3_5ResultsAlleleSpecificAnnotations", ".vcf");

        //Created by running
        //java -jar ~/bin/GenomeAnalysisTK-3.5.0/GenomeAnalysisTK.jar -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // --out as.gatk3.5.noDownsample.vcf -G StandardHC -G Standard -G AS_Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testVCFMode.gatk3.5.alleleSpecific.vcf");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_2bit_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-G", "StandardAnnotation",
                "-G", "AS_StandardAnnotation",
                "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.5 in AS VCF mode is < 99%");
    }

    /*
   * Test that in GVCF mode we're >= 99% concordant with GATK3 results
   */
    @Test
    public void testGVCFModeIsConcordantWithGATK3_5Results() throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGVCFModeIsConcordantWithGATK3Results", ".g.vcf");
        //Created by running:
        //java -jar ~/bin/GenomeAnalysisTK-3.5.0/GenomeAnalysisTK.jar -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // -ERC GVCF --out a.gatk3.5.noDownsample.g.vcf -G StandardHC -G Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testGVCFMode.gatk3.5.g.vcf");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_2bit_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-ERC", "GVCF",
                "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.5 in GVCF mode is < 99%");
    }


    @Test
    public void testGVCFModeIsConcordantWithGATK3_5AlelleSpecificResults() throws Exception {
        Utils.resetRandomGenerator();
        final File output = createTempFile("testGVCFModeIsConcordantWithGATK3_5AlelleSpecificResults", ".g.vcf");

        //Created by running:
        // java -jar ~/bin/GenomeAnalysisTK-3.5.0/GenomeAnalysisTK.jar -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // -ERC GVCF --out as.gatk3.5.noDownsample.g.vcf -G StandardHC -G Standard -G AS_Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testGVCFMode.gatk3.5.alleleSpecific.g.vcf");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_2bit_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-G", "StandardAnnotation",
                "-G", "AS_StandardAnnotation",
                "-ERC", "GVCF",
                "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.5 in AS GVCF mode is < 99%");
    }

    @Test
    public void testReferenceAdapterIsSerializable() throws IOException {
        final AuthHolder auth = new AuthHolder("name", "somestring");
        final ReferenceMultiSource referenceMultiSource = new ReferenceMultiSource(auth, b37_2bit_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION);
        SparkTestUtils.roundTripInKryo(referenceMultiSource, ReferenceMultiSource.class, SparkContextFactory.getTestSparkContext().getConf());
        final HaplotypeCallerSpark.ReferenceMultiSourceAdaptor adaptor = new HaplotypeCallerSpark.ReferenceMultiSourceAdaptor(referenceMultiSource, auth);
        SparkTestUtils.roundTripInKryo(adaptor, HaplotypeCallerSpark.ReferenceMultiSourceAdaptor.class, SparkContextFactory.getTestSparkContext().getConf());

    }

    @Test
    public void testGenotypeCalculationArgumentCollectionIsSerializable() {
        final GenotypeCalculationArgumentCollection args = new GenotypeCalculationArgumentCollection();
        SparkTestUtils.roundTripInKryo(args, GenotypeCalculationArgumentCollection.class, SparkContextFactory.getTestSparkContext().getConf());

    }

    @Test
    public void testHaplotypeCallerArgsAreSerializable() {
        final HaplotypeCallerArgumentCollection args = new HaplotypeCallerArgumentCollection();
        SparkTestUtils.roundTripInKryo(args, HaplotypeCallerArgumentCollection.class, SparkContextFactory.getTestSparkContext().getConf());
    }


    @Test
    public void testReferenceMultiSourceIsSerializable() {
        final ReferenceMultiSource args = new ReferenceMultiSource((PipelineOptions) null, BaseTest.b37_2bit_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION);
        SparkTestUtils.roundTripInKryo(args, ReferenceMultiSource.class, SparkContextFactory.getTestSparkContext().getConf());
    }


    @Test
    public void testBroadcastHcArgs() {
        SparkContextFactory.getTestSparkContext().broadcast(new HaplotypeCallerArgumentCollection());
    }

    @Test
    public void testFastGenotypeIsSerializable() {
        SparkTestUtils.roundTripInKryo(Collections.nCopies(2, "value"), Collections.nCopies(2, "value").getClass(), SparkContextFactory.getTestSparkContext().getConf());
    }
}
