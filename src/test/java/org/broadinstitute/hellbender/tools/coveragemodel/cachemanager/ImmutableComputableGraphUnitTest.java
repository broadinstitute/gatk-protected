package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import junit.framework.AssertionFailedError;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link ImmutableComputableGraph}
 *
 * TODO github/gatk-protected issue #803
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class ImmutableComputableGraphUnitTest extends BaseTest {

//    public static ComputableNodeFunction f_computation_function = new ComputableNodeFunction() {
//        @Override
//        public Duplicable apply(Map<String, Duplicable> parents) throws ParentValueNotFoundException {
//            final INDArray x = fetchINDArray("x", parents);
//            final double y = fetch()
//            return new DuplicableNDArray(x.add(y));
//        }
//    };

//    public static Function<Map<String, ? extends Duplicable>, ? extends Duplicable> g_computation_function = parents -> {
//        final INDArray y = DuplicableNDArray.of(col.get("X_prod_Y"));
//        final double y = DuplicableNumber.of(col.get("Y"));
//        return new DuplicableNDArray(xProdY.add(y));
//    };

    @Test(enabled = false)
    public void testMissingParents() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test(enabled = false)
    public void testDuplicateNodeKeys() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test(enabled = false)
    public void testCyclicGraphException() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test(enabled = false)
    public void testOutdatedCaches() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test(enabled = false)
    public void testUpToDateCaches() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test(enabled = false)
    public void testComputeOnDemandNodes() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test(enabled = false)
    public void testPrimitiveUpdating() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test(enabled = false)
    public void testExternallyComputableUpdating() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test(enabled = false)
    public void testCacheByTag() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test(enabled = false)
    public void testCacheByNode() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test(enabled = false)
    public void testCacheAutoUpdate() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test(enabled = false)
    public void testUnchangedNodesSameReferenceAfterUpdate() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test(enabled = false)
    public void testNoRedundantComputation() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    public void testUninitializedPrimitiveNode() {

    }

    public void testUninitializedExternallyComputedNode() {

    }
}
