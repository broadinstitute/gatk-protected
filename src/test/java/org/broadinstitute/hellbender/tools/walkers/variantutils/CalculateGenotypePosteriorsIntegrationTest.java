package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Collections;

public final class CalculateGenotypePosteriorsIntegrationTest extends CommandLineProgramTest {

    private final String largeDir = largeFileTestDir;
    private final String dir= getToolTestDataDir();

    private String CEUtrioFamilyFile = dir + "CEUtrio.ped";
    private String threeMemberNonTrioFamilyFile = dir + "threeMemberNonTrio.ped";

    //Note: these files were made by taking GATK3 files and changing the chromosome from 1 to 20
    // (because our reference is on 20 and it does not matter for this tool anyway)
    private String CEUtrioTest = dir + "CEUtrioTest_chr20.vcf";
    private String CEUtrioPopPriorsTest = dir + "CEUtrioPopPriorsTest_chr20.vcf";
    private String getThreeMemberNonTrioTest = dir + "threeMemberNonTrioTest_chr20.vcf";

    @Test
    //use the first 20 variants to save time; they have a nice range of AC from 4 to over 4000
    public void testUsingDiscoveredAF() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -O %s" +
                        " -R " + b37_reference_20_21 +    //NOTE: we need a reference for -L
                        " -L 20:10,000,000-10,001,432" +
                        " -V " + largeDir + "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf",
                Collections.singletonList(dir + "expectedCGP_testUsingDiscoveredAF.vcf")
        );
        spec.executeTest("testUsingDiscoveredAF", this);
    }

    @Test
    //this test ignores discovered AC and has no external priors, so it should only apply the PP tag with values the same as PLs
    //only test the first 20 variants to save time
    public void testMissingPriors() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        "-useACoff" +
                        " -O %s" +
                        " -R " + b37_reference_20_21 +    //NOTE: we need a reference for -L
                        " -L 20:10,000,000-10,001,432" +
                        " -V " + largeDir + "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf",
                Collections.singletonList(dir + "expectedCGP_testMissingPriors.vcf")
        );
        spec.executeTest("testMissingPriors", this);
    }

    @Test
    public void testInputINDELs() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        "-useACoff" +
                        " -O %s" +
                        " -R " + b37_reference_20_21 +    //NOTE: we need a reference for -L
                        " -L 20:10,000,000-10,100,000" +
                        " -V " + dir + "NA12878.Jan2013.haplotypeCaller.subset.indels.vcf" +
                        " -supporting " + largeDir + "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf",
                Collections.singletonList(dir + "expectedCGP_testInputINDELs.vcf")
        );
        spec.executeTest("testInputINDELs", this);
    }

    @Test
    public void testFamilyPriors() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        "-useACoff" +
                        " -O %s" +
                        " -ped " + CEUtrioFamilyFile +
                        " -V " + CEUtrioTest +
                        " -supporting " + CEUtrioPopPriorsTest,
                Collections.singletonList(dir + "expectedCGP_testFamilyPriors.vcf")
        );
        spec.executeTest("testFamilyPriors", this);
    }

    @Test
    public void testSingleParentFamily() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -O %s" +
                        " -ped " + threeMemberNonTrioFamilyFile +
                        " -V " + getThreeMemberNonTrioTest +
                        " -skipPop",
                Collections.singletonList(dir + "expectedCGP_testSingleParentFamily.vcf")
        );
        spec.executeTest("testFamilyPriors", this);
    }
}
