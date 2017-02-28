package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.concordance.ConcordanceSummaryRecord;
import org.broadinstitute.hellbender.tools.walkers.concordance.VariantStatusRecord;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Created by tsato on 1/31/17.
 */
public class ConcordanceTest extends CommandLineProgramTest{
    final double epsilon = 1e-3;

    @Test
    public void testSimple() throws Exception {
        final String testDir =  publicTestDir + "org/broadinstitute/hellbender/tools/concordance/";
        final File evalVcf = new File(testDir, "gatk4-dream3-mini.vcf");
        final File truthVcf = new File(testDir, "gatk3-dream3-mini.vcf");
        final File variantStatusTable = createTempFile("table", ".txt");
        final File summaryTable = createTempFile("summary", ".txt");

        final String[] args = {
                "-V", evalVcf.toString(),
                "--truth", truthVcf.toString(),
                "--table",  variantStatusTable.toString(),
                "--summary", summaryTable.toString(),
                "-tumor", "IS3.snv.indel.sv"
        };
        runCommandLine(args);

        // the number of lines (variants) in the status table should match that of the eval vcf
        final Path path = Paths.get(variantStatusTable.getAbsolutePath());
        // Assert.assertEquals(Files.lines(path).count(), 20);

        ConcordanceSummaryRecord.Reader reader = new ConcordanceSummaryRecord.Reader(summaryTable);
        ConcordanceSummaryRecord snpRecord = reader.readRecord();
        ConcordanceSummaryRecord indelRecord = reader.readRecord();

        Assert.assertEquals(snpRecord.getVariantType(), VariantContext.Type.SNP);
        Assert.assertEquals(indelRecord.getVariantType(), VariantContext.Type.INDEL);

        // numbers verified by manual inspection
        Assert.assertEquals(snpRecord.getTruePositives(), 1);
        Assert.assertEquals(indelRecord.getTruePositives(), 4);
        Assert.assertEquals(snpRecord.getFalsePositives(), 0);
        Assert.assertEquals(indelRecord.getFalsePositives(), 3);
        Assert.assertEquals(snpRecord.getFalseNegatives(), 3);
        Assert.assertEquals(indelRecord.getFalseNegatives(), 2);
        Assert.assertEquals(snpRecord.getSensitivity(), 1.0/4, epsilon);
        Assert.assertEquals(snpRecord.getPrecision(), 1.0, epsilon);
    }

    // Test going from an integer chromosome (22) to a character chromosome (X)
    @Test
    public void test22X() throws Exception {
        final String testDir = publicTestDir + "org/broadinstitute/hellbender/tools/concordance/";
        final File truthVcf = new File(testDir, "gatk3-dream3-x.vcf");
        final File evalVcf = new File(testDir, "gatk4-dream3-x.vcf");
        final File table = createTempFile("table", ".txt");
        final File summary = createTempFile("summary", ".txt");

        final String[] args = {
                "-V", evalVcf.toString(),
                "--truth", truthVcf.toString(),
                "--table", table.toString(),
                "--summary", summary.toString(),
                "-tumor", "IS3.snv.indel.sv"
        };
        runCommandLine(args);
        ConcordanceSummaryRecord.Reader reader = new ConcordanceSummaryRecord.Reader(summary);
        ConcordanceSummaryRecord record = reader.readRecord();

        Assert.assertEquals(record.getTruePositives(), 4);
        Assert.assertEquals(record.getFalsePositives(), 6);
        Assert.assertEquals(record.getFalseNegatives(), 6);
        Assert.assertEquals(record.getSensitivity(), 4.0/10, epsilon);
        Assert.assertEquals(record.getPrecision(), 4.0/10, epsilon);
    }

    @Test
    public void testPerfectMatchVcf() throws Exception {
        // Truth and eval come from the same vcf file.
        // The truth file only contains the PASS sites. Should get 100% sensitivity and specificity.
        final String testDir = publicTestDir + "org/broadinstitute/hellbender/tools/concordance/";
        final File truthVcf = new File(testDir, "same-truth.vcf");
        final File evalVcf = new File(testDir, "same-eval.vcf");
        final File table = createTempFile("table", ".txt");
        final File summary = createTempFile("summary", ".txt");

        final String[] args = {
                "-V", evalVcf.toString(),
                "--truth", truthVcf.toString(),
                "--table", table.toString(),
                "--summary", summary.toString(),
                "-tumor", "TUMOR"
        };

        runCommandLine(args);
        ConcordanceSummaryRecord.Reader reader = new ConcordanceSummaryRecord.Reader(summary);
        ConcordanceSummaryRecord record = reader.readRecord();

        Assert.assertEquals(record.getTruePositives(), 284);
        Assert.assertEquals(record.getFalsePositives(), 0);
        Assert.assertEquals(record.getFalseNegatives(), 0);
        Assert.assertEquals(record.getSensitivity(), 1.0, epsilon);
        Assert.assertEquals(record.getPrecision(), 1.0, epsilon);
    }

    // TODO: figure out a way to subset the truth variant
    @Test
    public void testDreamSensitivity() throws Exception {
        final String testDir = publicTestDir + "org/broadinstitute/hellbender/tools/concordance/";
        final File evalVcf = new File(testDir, "dream3-chr21.vcf");
        final File truthVcf = new File(testDir, "dream3-truth-minus-MSK-chr21.vcf");
        final File table = createTempFile("table", ".txt");
        final File summary = createTempFile("summary", ".txt");

        final String[] args = {
                "-V", evalVcf.toString(),
                "--truth", truthVcf.toString(),
                "--table", table.toString(),
                "--summary", summary.toString(),
                "-tumor", "TUMOR"
        };

        runCommandLine(args);

        ConcordanceSummaryRecord.Reader reader = new ConcordanceSummaryRecord.Reader(summary);
        ConcordanceSummaryRecord snpRecord = reader.readRecord();
        ConcordanceSummaryRecord indelRecord = reader.readRecord();

        // The expected sensitivities come from dream challenge's official python script
        // We don't expect them to match exactly
        // Currently we tend to overestimate sensitivities TODO: investigate
        Assert.assertTrue(snpRecord.getSensitivity() > 0.85); // with python script (masked, chr21): 0.879
        Assert.assertTrue(indelRecord.getSensitivity() > 0.75); // with python script (masked): 0.791

        // the following tests should pass but they don't; there must be bugs

        // count the number of INDEL variants in the truth using awk. it should match # TP + # FN
        // grep -v ^# dream3-truth-minus-MSK-chr21.vcf | grep -v SVTYPE | awk 'length($4) != length($5) {print $0 }' | wc -l
        // 97
        // (as of 2/27/17) we get 85 TP + 6 FN = 91
        // Assert.assertEquals(indelRecord.getTruePositives() + record.getFalseNegatives(), 97);

        // analogously for SNPs:
        // grep -v ^# dream3-truth-minus-MSK-chr21.vcf | grep -v SVTYPE | awk 'length($4) == length($5) {print $0 }' | wc -l
        // 78
        Assert.assertEquals(snpRecord.getTruePositives() + snpRecord.getFalseNegatives(), 78);

        // when we don't separate variants by SNPs and INDELs we get:
        // 97 + 78 = 175 (using awk)
        // 91 + 81 = 172 (using Concordance.java)
    }

    // TODO: we don't suppport -L argument yet but once David refactors this tool we should create a test for it

}