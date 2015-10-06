package org.broadinstitute.hellbender.utils.param;

import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Created by lichtens on 8/25/15.
 */
public class ParamUtilsUnitTest {

    @Test
    public void testInRangeSuccess(){
        Assert.assertTrue(4 == ParamUtils.inRange(4, 3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(4.1 == ParamUtils.inRange(4.1, 3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(4.1 == ParamUtils.inRange(4.1, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(0 == ParamUtils.inRange(0, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(0.0 == ParamUtils.inRange(0.0, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(0 == ParamUtils.inRange(0, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(0.0 == ParamUtils.inRange(0.0, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(-1 == ParamUtils.inRange(-1, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(-1.5 == ParamUtils.inRange(-1.5, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(-1 == ParamUtils.inRange(-1, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureAllPositiveInt(){
        ParamUtils.inRange(4, 7, 10, "Range calculation did not work properly");
    }
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureAllPositiveDouble(){
        ParamUtils.inRange(4.1, 7.2, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureMinNegativeDouble(){
        ParamUtils.inRange(45.2, -7.2, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureValNegativeDouble(){
        ParamUtils.inRange(-10.3, -7.2, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureValNegativeInt(){
        ParamUtils.inRange(-10, -7, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureValNegativeInt2(){
        ParamUtils.inRange(-10, -7.2, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureValMinMaxUnreasonable(){
        ParamUtils.inRange(5, 10, -7, "Range calculation did not work properly.  Min was greater than max, so will always throw exception.");
    }

    @Test
    public void testInRangeLongDoubleDouble(){
        ParamUtils.inRange(5, -1.2, 7.1, "Range calculation did not work properly!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is infinite!")
    public void testInfinite(){
        ParamUtils.isFinite(Double.POSITIVE_INFINITY, "That is infinite!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is infinite!")
    public void testNegInfinite(){
        ParamUtils.isFinite(Double.NEGATIVE_INFINITY, "That is infinite!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is negative!")
    public void testNegativeIsNotPositive(){
        ParamUtils.isPositiveOrZero(Double.NEGATIVE_INFINITY, "That is negative!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is negative!")
    public void testNegativeIsNotPositive2(){
        ParamUtils.isPositiveOrZero(-4.2, "That is negative!");
    }
}
