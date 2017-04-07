package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.commons.codec.digest.DigestUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.*;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

public class GenotypeGVCFsIntegrationTest extends CommandLineProgramTest {

    private static final List<String> NO_EXTRA_ARGS = Collections.emptyList();
    private static final String b38_reference_20_21 =largeFileTestDir + "Homo_sapiens_assembly38.20.21.fasta";

    private static <T> void assertForEachElementInLists(final List<T> actual, final List<T> expected, final BiConsumer<T, T> assertion) {
        Assert.assertEquals(actual.size(), expected.size(), "different number of elements in lists:\n"
                + actual.stream().map(Object::toString).collect(Collectors.joining("\n","actual:\n","\n"))
        +  expected.stream().map(Object::toString).collect(Collectors.joining("\n","expected:\n","\n")));
        for (int i = 0; i < actual.size(); i++) {
            assertion.accept(actual.get(i), expected.get(i));
        }
    }

    @DataProvider(name = "gvcfsToGenotype")
    public Object[][] gvcfsToGenotype() {
        String basePairGVCF = "gvcf.basepairResolution.gvcf";
        return new Object[][]{
//                {"combine.single.sample.pipeline.1.vcf", null, Arrays.asList("-V", getTestFile("combine.single.sample.pipeline.2.vcf").toString() , "-V", getTestFile("combine.single.sample.pipeline.3.vcf").toString()), b37_reference_20_21}, //combine not supported yet
                {getTestFile(basePairGVCF), getTestFile( "gvcf.basepairResolution.output.vcf"), NO_EXTRA_ARGS, b37_reference_20_21}, //base pair level gvcf
                {getTestFile("testUpdatePGT.gvcf"), getTestFile( "testUpdatePGT.output.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},   //testUpdatePGT
                {getTestFile("gvcfExample1.vcf"), getTestFile( "gvcfExample1.vcf.expected.vcf"), NO_EXTRA_ARGS, b37_reference_20_21}, //single sample vcf
                {getTestFile("gvcfExample1.vcf"), getTestFile( "gvcfExample1.vcf.expected.vcf"), Arrays.asList("-L", "20"), b37_reference_20_21}, //single sample vcf with -L
                {getTestFile("combined_genotype_gvcf_exception.vcf"), getTestFile( "combined_genotype_gvcf_exception.output.vcf"), NO_EXTRA_ARGS, b37_reference_20_21}, //test that an input vcf with 0/0 already in GT field is overwritten
                {getTestFile("combined_genotype_gvcf_exception.original.vcf"), getTestFile( "combined_genotype_gvcf_exception.output.vcf"), NO_EXTRA_ARGS, b37_reference_20_21}, //test that an input vcf with 0/0 already in GT field is overwritten
                {getTestFile("combined_genotype_gvcf_exception.nocall.vcf"), getTestFile( "combined_genotype_gvcf_exception.output.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},  //same test as above but with ./.
                {getTestFile(basePairGVCF), getTestFile( "ndaTest.expected.vcf"), Collections.singletonList("-nda"), b37_reference_20_21},  //annotating with the number of alleles discovered option
                {getTestFile(basePairGVCF), getTestFile( "maxAltAllelesTest.expected.vcf"), Arrays.asList("-maxAltAlleles", "1"), b37_reference_20_21 }, //restricting the max number of alt alleles
                {getTestFile(basePairGVCF), getTestFile( "standardConfTest.expected.vcf"), Arrays.asList("-stand_call_conf","300"),b37_reference_20_21}, //set minimum calling threshold
                {getTestFile("spanningDel.combined.g.vcf"), getTestFile( "spanningDel.combined.g.vcf.expected.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},
                {getTestFile("spanningDel.delOnly.g.vcf"), getTestFile( "spanningDel.delOnly.g.vcf.expected.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},
                {getTestFile("spanningDel.depr.delOnly.g.vcf"), getTestFile( "spanningDel.depr.delOnly.g.vcf.expected.vcf" ), NO_EXTRA_ARGS, b37_reference_20_21},
                {getTestFile("ad-bug-input.vcf"), getTestFile( "ad-bug-output.vcf"), NO_EXTRA_ARGS, b37_reference_20_21}, //Bad AD Propagation Haploid Bug
                {getTestFile("CEUTrio.20.21.gatk3.4.g.vcf"), getTestFile( "CEUTrio.20.21.expected.vcf"), Arrays.asList("--dbsnp", "src/test/resources/large/dbsnp_138.b37.20.21.vcf"), b37_reference_20_21},
                {getTestFile("CEUTrio.20.21.missingIndel.g.vcf"), getTestFile( "CEUTrio.20.21.missingIndel.expected.vcf"), Arrays.asList("--dbsnp", "src/test/resources/large/dbsnp_138.b37.20.21.vcf"), b37_reference_20_21},
                {new File(largeFileTestDir + "gvcfs/cutDown.24_sample.21.g.vcf"), new File( largeFileTestDir + "gvcfs/cutDown.24_sample.21.expected.vcf"), NO_EXTRA_ARGS, b38_reference_20_21},
                {getTestFile("chr21.bad.pl.g.vcf"), getTestFile( "chr21.bad.pl.expected.vcf"), Arrays.asList("-L", "chr21:28341770-28341790"), b38_reference_20_21},
                {new File(largeFileTestDir + "gvcfs/combined.gatk3.7.g.vcf.gz"),  new File(largeFileTestDir + "gvcfs/combined.gatk3.7.expected.vcf.gz"), NO_EXTRA_ARGS, b38_reference_20_21}
                //{getTestFile(basePairGVCF), getTestFile( "gvcf.basepairResolution.includeNonVariantSites.expected.vcf"), Collections.singletonList("--includeNonVariantSites") //allsites not supported yet
        };
    }


    @Test(dataProvider = "gvcfsToGenotype", enabled = false) // this is useful for development
    public void compareToGATK3(File input, File outputFile, List<String> extraArgs, String reference) throws IOException, NoSuchAlgorithmException {
        final String GATK3_PATH = "gatk3.7-30-ga4f720357.jar";
        final String params = GATK3_PATH + input.getAbsolutePath() + extraArgs.stream().collect(Collectors.joining()) + reference;
        final String md5 = DigestUtils.md5Hex(params);
        final File gatk3ResultsDir = new File("gatk3results");
        if(! gatk3ResultsDir.exists()){
            Assert.assertTrue(gatk3ResultsDir.mkdir());
        }
        final File gatk3Result = new File(gatk3ResultsDir, md5 + ".vcf");
        if (!gatk3Result.exists()) {
            List<String> gatk3Command = new ArrayList<>(
                    Arrays.asList("java", "-jar", GATK3_PATH, "-T", "GenotypeGVCFs"));
            gatk3Command.add("-V");
            gatk3Command.add(input.getAbsolutePath());
            gatk3Command.add("-o");
            gatk3Command.add(gatk3Result.getAbsolutePath());
            gatk3Command.add("-R");
            gatk3Command.add(reference);
            gatk3Command.addAll(extraArgs);

            runProcess(new ProcessController(), gatk3Command.toArray(new String[gatk3Command.size()]));
        } else {
            System.out.println("Found precomputed gatk3Result");
        }
        final Path expectedResultsDir = Paths.get("expectedResults");
        if ( !Files.exists(expectedResultsDir)) {
            Files.createDirectory(expectedResultsDir);
        }
        Files.copy(gatk3Result.toPath(), expectedResultsDir.resolve(outputFile.getName()), StandardCopyOption.REPLACE_EXISTING);

        assertGenotypesMatch(input, gatk3Result, extraArgs, reference);
        assertVariantContextsMatch(input, gatk3Result, extraArgs, reference);
    }


    private static void runProcess(ProcessController processController, String[] command) {
            final ProcessSettings prs = new ProcessSettings(command);
            prs.getStderrSettings().printStandard(true);
            prs.getStdoutSettings().printStandard(true);
            final ProcessOutput output = processController.exec(prs);
            Assert.assertEquals(output.getExitValue(), 0, "Process exited with non-zero value. Command: "+ Arrays.toString(command) + "\n");
    }

    @Test(dataProvider = "gvcfsToGenotype")
    public void testGenotypesOnly(File input, File expected, List<String> extraArgs, String reference) throws IOException {
        assertGenotypesMatch(input, expected, extraArgs, reference);
    }

    @Test(dataProvider = "gvcfsToGenotype")
    public void testEntireVariantContext(File input, File expected, List<String> extraArgs, String reference) throws IOException {
        assertVariantContextsMatch(input, expected, extraArgs, reference);
    }

    private void assertVariantContextsMatch(File input, File expected, List<String> extraArgs, String reference) throws IOException {
        runGenotypeGVCFSAndAssertSomething(input, expected, extraArgs, (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e,
                Arrays.asList("QD",//TODO QD has a cap value and anything that reaches that is randomized.  It's difficult to reproduce the same random numbers accross gatk3 -> 4
                              "FS")),//TODO There's some bug in either gatk3 or gatk4 fisherstrand that's making them not agree still, I'm not sure which is correct
                 reference);
    }

    private void assertGenotypesMatch(File input, File expected, List<String> additionalArguments, String reference) throws IOException {
        runGenotypeGVCFSAndAssertSomething(input, expected, additionalArguments, VariantContextTestUtils::assertVariantContextsHaveSameGenotypes,
                                           reference);
    }

    private void runGenotypeGVCFSAndAssertSomething(File input, File expected, List<String> additionalArguments, BiConsumer<VariantContext, VariantContext> assertion, String reference) throws IOException {
        final File output = createTempFile("genotypegvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(reference))
                .addVCF(input)
                .addOutput(output);

        additionalArguments.forEach(args::add);

        Utils.resetRandomGenerator();
        runCommandLine(args);

        final List<VariantContext> expectedVC = getVariantContexts(expected);
        final List<VariantContext> actualVC = getVariantContexts(output);
        assertForEachElementInLists(actualVC, expectedVC, assertion); //TODO Inbreeding calculation changed between 3.4 and now
    }

    /**
     * Returns a list of VariantContext records from a VCF file
     *
     * @param vcfFile VCF file
     * @return list of VariantContext records
     * @throws IOException if the file does not exist or can not be opened
     */
    private static List<VariantContext> getVariantContexts(final File vcfFile) throws IOException {
        final VCFCodec codec = new VCFCodec();
        final FileInputStream s = new FileInputStream(vcfFile);
        final LineIterator lineIteratorVCF = codec.makeSourceFromStream(new PositionalBufferedStream(s));
        codec.readHeader(lineIteratorVCF);

        final List<VariantContext> VCs = new ArrayList<>();
        while (lineIteratorVCF.hasNext()) {
            final String line = lineIteratorVCF.next();
            Assert.assertFalse(line == null);
            VCs.add(codec.decode(line));
        }

        return VCs;
    }
}
