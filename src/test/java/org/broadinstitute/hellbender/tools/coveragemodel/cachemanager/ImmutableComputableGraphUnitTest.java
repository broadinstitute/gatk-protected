package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import avro.shaded.com.google.common.collect.Sets;
import junit.framework.AssertionFailedError;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Set;

/**
 * Unit tests for {@link ImmutableComputableGraph}
 *
 *
 * Example:
 *
 *     x    y    z
 *      \  / \  /
 *       \/   \/
 *       f    g
 *        \  /
 *         \/
 *         h
 *
 *  x stores DuplicableNDArray
 *  y stores DuplicableNumber
 *  z stores DuplicableNDArray
 *
 *  f = x+y
 *  g = y.z
 *  h = log(f.g) = (x+y).(y.z)
 *
 *  All operations are point-wise
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class ImmutableComputableGraphUnitTest extends BaseTest {

    /**
     * A simple helper class for keeping track of function evaluations
     */
    private static class Counter {
        final Map<String, Integer> counts;

        Counter(String ... keys) {
            counts = new HashMap<>();
            for (final String key : keys) {
                counts.put(key, 0);
            }
        }

        private Counter(final Map<String, Integer> otherCounts) {
            counts = new HashMap<>(otherCounts.size());
            counts.putAll(otherCounts);
        }

        public void increment(final String key) {
            counts.put(key, getCount(key) + 1);
        }

        public int getCount(final String key) {
            return counts.get(key);
        }

        public Set<String> getKeys() {
            return counts.keySet();
        }

        public Counter copy() {
            return new Counter(counts);
        }

        public Counter diff(final Counter oldCounter) {
            Utils.validateArg(Sets.symmetricDifference(oldCounter.getKeys(), getKeys()).isEmpty(),
                    "the counters must have the same keys");
            final Map<String, Integer> diffMap = new HashMap<>(getKeys().size());
            getKeys().forEach(key -> diffMap.put(key, getCount(key) - oldCounter.getCount(key)));
            return new Counter(diffMap);
        }
    }

    private static Counter counter = new Counter("f", "g", "h");
    private static int[] TEST_INDARRAY_SHAPE = new int[] {2, 3};
    private static Random rng = new Random(1984);

    private static INDArray getRandomINDArray() {
        return Nd4j.rand(TEST_INDARRAY_SHAPE);
    }

    private static double getRandomDouble() {
        return rng.nextDouble();
    }

    private static Counter getCounterInstance() {
        return counter.copy();
    }

    public static ComputableNodeFunction f_computation_function = new ComputableNodeFunction() {
        @Override
        public Duplicable apply(Map<String, Duplicable> parents) throws ParentValueNotFoundException {
            final INDArray x = fetchINDArray("x", parents);
            final double y = fetchDouble("y", parents);
            counter.increment("f");
            return new DuplicableNDArray(x.add(y));
        }
    };

    public static ComputableNodeFunction g_computation_function = new ComputableNodeFunction() {
        @Override
        public Duplicable apply(Map<String, Duplicable> parents) throws ParentValueNotFoundException {
            final double y = fetchDouble("y", parents);
            final INDArray z = fetchINDArray("z", parents);
            counter.increment("g");
            return new DuplicableNDArray(z.mul(y));
        }
    };

    public static ComputableNodeFunction h_computation_function = new ComputableNodeFunction() {
        @Override
        public Duplicable apply(Map<String, Duplicable> parents) throws ParentValueNotFoundException {
            final INDArray f = fetchINDArray("f", parents);
            final INDArray g = fetchINDArray("g", parents);
            counter.increment("h");
            return new DuplicableNDArray(f.mul(g));
        }
    };

    public static INDArray g_external_computer(final double y, final INDArray z) {
        return z.mul(y);
    }

    /**
     * Calculates f, g, and h directly
     */
    public static ImmutableTriple<INDArray, INDArray, INDArray> CALCULATE_DIRECT(final INDArray x, final double y,
                                                                                 final INDArray z) {
        final INDArray f = x.add(y);
        final INDArray g = z.mul(y);
        final INDArray h = f.mul(g);
        return ImmutableTriple.of(f, g, h);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMissingParents() {
        ImmutableComputableGraph.builder()
                .addPrimitiveNode("x", new String[] {}, new DuplicableNDArray())
                .addPrimitiveNode("y", new String[] {}, new DuplicableNumber<Double>())
                .addPrimitiveNode("z", new String[] {}, new DuplicableNDArray())
                /* q is undefined */
                .addComputableNode("f", new String[] {}, new String[] {"x", "y", "q"}, f_computation_function, true)
                .addComputableNode("g", new String[] {}, new String[] {"y", "z"}, g_computation_function, true)
                .addComputableNode("h", new String[] {}, new String[] {"f", "g"}, h_computation_function, true)
                .build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDuplicatePrimitiveNode() {
        ImmutableComputableGraph.builder()
                .addPrimitiveNode("x", new String[] {}, new DuplicableNDArray())
                .addPrimitiveNode("x", new String[] {}, new DuplicableNDArray())
                .build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDuplicateComputableNode() {
        ImmutableComputableGraph.builder()
                .addComputableNode("x", new String[] {}, new String[] {}, null, true)
                .addComputableNode("x", new String[] {}, new String[] {}, null, true)
                .build();
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
