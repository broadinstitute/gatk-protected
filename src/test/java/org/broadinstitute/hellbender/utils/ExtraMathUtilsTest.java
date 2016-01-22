package org.broadinstitute.hellbender.utils;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

import static org.testng.Assert.*;

/**
 * Created by davidben on 1/22/16.
 */
public class ExtraMathUtilsTest {
    @Test
    public void naturalLogSumNaturalLogTest() {
        final double[] realSpaceValues = {0.1, 0.2, 1.4, 5.9};
        final double[] logSpaceValues= Arrays.stream(realSpaceValues).map(Math::log).toArray();

        final double expected = Math.log(Arrays.stream(realSpaceValues).sum());
        final double actual = ExtraMathUtils.naturalLogSumNaturalLog(logSpaceValues);
        Assert.assertEquals(actual, expected, 1e-10);
    }

}