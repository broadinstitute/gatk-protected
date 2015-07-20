package org.broadinstitute.hellbender.tools.exome;


import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Created by davidben on 7/1/15.
 */
public class TargetCoverageTest extends BaseTest {

    @Test
    public void testTargetCoverageConstructor() {
        SimpleInterval interval = new SimpleInterval("chr",1,2);
        String name = "arbitrary";
        double coverage = 0.1;

        TargetCoverage target = new TargetCoverage(name, interval, coverage);
        Assert.assertEquals(target.getName(), name);
        Assert.assertEquals(target.getInterval(), interval);
        Assert.assertEquals(target.getCoverage(), coverage);
    }

    @Test
    public void testReadTargetsWithCoverage() {
        final File TEST_DIR = new File("src/test/resources/org/broadinstitute/tools/exome/caller");
        final File TEST_TARGETS = new File(TEST_DIR,"targets.tsv");

        try {
            List<TargetCoverage> targets = TargetCoverageUtils.readTargetsWithCoverage(TEST_TARGETS);
            Assert.assertEquals(targets.size(), 52);
            Assert.assertEquals(targets.get(0).getName(), "arbitrary_name");
            Assert.assertEquals(targets.get(1).getInterval().getContig(), "chr");
            Assert.assertEquals(targets.get(2).getCoverage(), 1.0);

        } catch (IOException e) {
            Assert.fail(e.getMessage());
        }
    }
}
