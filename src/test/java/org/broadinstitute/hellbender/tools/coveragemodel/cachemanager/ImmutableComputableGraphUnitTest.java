package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import avro.shaded.com.google.common.collect.ImmutableMap;
import avro.shaded.com.google.common.collect.Sets;
import junit.framework.AssertionFailedError;
import org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jApacheAdapterUtils;
import org.broadinstitute.hellbender.utils.MathObjectAsserts;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link ImmutableComputableGraph}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class ImmutableComputableGraphUnitTest extends BaseTest {

    private static final Random rng = new Random(1984);

    private static final double EPSILON = 1e-16;

    private static final int NUM_TRIALS = 10;

    /*
     * A working example:
     *
     *     x    y    z
     *     /\  / \  /
     *     \ \/   \/
     *      \f    g
     *       \\  /
     *        \\/
     *         h
     *
     *  x stores DuplicableNDArray
     *  y stores DuplicableNumber
     *  z stores DuplicableNDArray
     *
     *  f = x+y
     *  g = y.z
     *  h = f.g - x
     *
     *  All operations are point-wise
     */

    private static final Counter counter = new Counter("f", "g", "h");

    /**
     * Shape of "x" and "z"
     */
    private static final int[] TEST_NDARRAY_SHAPE = new int[] {2, 3};

    private static INDArray getRandomINDArray() {
        return Nd4j.rand(TEST_NDARRAY_SHAPE);
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
            final INDArray x = fetchINDArray("x", parents);
            counter.increment("h");
            return new DuplicableNDArray(f.mul(g).sub(x));
        }
    };

    private static ImmutableComputableGraphUtils.ImmutableComputableGraphBuilder getTestICGBuilder(
            final boolean f_caching, final boolean f_external,
            final boolean g_caching, final boolean g_external,
            final boolean h_caching, final boolean h_external,
            final String[] x_tags, final String[] y_tags, final String[] z_tags,
            final String[] f_tags, final String[] g_tags, final String[] h_tags) {
        return ImmutableComputableGraph.builder()
                .addPrimitiveNode("x", x_tags, new DuplicableNDArray())
                .addPrimitiveNode("y", y_tags, new DuplicableNumber<Double>())
                .addPrimitiveNode("z", z_tags, new DuplicableNDArray())
                .addComputableNode("f", f_tags, new String[]{"x", "y"},
                        f_external ? null : f_computation_function, f_caching)
                .addComputableNode("g", g_tags, new String[]{"y", "z"},
                        g_external ? null : g_computation_function, g_caching)
                .addComputableNode("h", h_tags, new String[]{"f", "g", "x"},
                        h_external ? null : h_computation_function, h_caching);
    }

    private static ImmutableComputableGraphUtils.ImmutableComputableGraphBuilder getTestICGBuilder(
            final boolean f_caching, final boolean f_external,
            final boolean g_caching, final boolean g_external,
            final boolean h_caching, final boolean h_external) {
        return getTestICGBuilder(f_caching, f_external, g_caching, g_external, h_caching, h_external,
                new String[] {}, new String[] {}, new String[] {},
                new String[] {}, new String[] {}, new String[] {});
    }

    /**
     * Calculates f, g, and h directly
     */
    private static void assertCorrectness(final INDArray x, final double y, final INDArray z,
                                          final INDArray xICG, final double yICG, final INDArray zICG,
                                          final INDArray fICG, final INDArray gICG, final INDArray hICG) {
        final INDArray fExpected = x.add(y);
        final INDArray gExpected = z.mul(y);
        final INDArray hExpected = fExpected.mul(gExpected).sub(x);

        MathObjectAsserts.assertRealMatrixEquals(
                Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(xICG),
                Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(x));
        MathObjectAsserts.assertRealMatrixEquals(
                Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(zICG),
                Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(z));
        Assert.assertEquals(yICG, y, EPSILON);
        MathObjectAsserts.assertRealMatrixEquals(
                Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(fICG),
                Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(fExpected));
        MathObjectAsserts.assertRealMatrixEquals(
                Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(gICG),
                Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(gExpected));
        MathObjectAsserts.assertRealMatrixEquals(
                Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(hICG),
                Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(hExpected));
    }

    /**
     * Tests a fully automated auto-updating ICG
     */
    @Test
    public void testAutoUpdateCache() {
        final ImmutableComputableGraph icg = getTestICGBuilder(true, false, true, false, true, false)
                .enableCacheAutoUpdate().build();
        for (int i = 0; i < NUM_TRIALS; i++) {
            final INDArray x = getRandomINDArray();
            final double y = getRandomDouble();
            final INDArray z = getRandomINDArray();

            final Counter startCounts = getCounterInstance();
            final ImmutableComputableGraph initializedICG = icg
                    .setValue("x", new DuplicableNDArray(x))
                    .setValue("y", new DuplicableNumber<>(y))
                    .setValue("z", new DuplicableNDArray(z));
            final INDArray xICG = DuplicableNDArray.strip(initializedICG.getValueDirect("x"));
            final double yICG = DuplicableNumber.strip(initializedICG.getValueDirect("y"));
            final INDArray zICG = DuplicableNDArray.strip(initializedICG.getValueDirect("z"));
            final INDArray fICG = DuplicableNDArray.strip(initializedICG.getValueDirect("f"));
            final INDArray gICG = DuplicableNDArray.strip(initializedICG.getValueDirect("g"));
            final INDArray hICG = DuplicableNDArray.strip(initializedICG.getValueDirect("h"));
            final Counter diffCounts = getCounterInstance().diff(startCounts);

            assertCorrectness(x, y, z, xICG, yICG, zICG, fICG, gICG, hICG);

            /* each function must be calculated only once; otherwise, ICG is doing redundant computations */
            Assert.assertEquals(diffCounts.getCount("f"), 1);
            Assert.assertEquals(diffCounts.getCount("g"), 1);
            Assert.assertEquals(diffCounts.getCount("h"), 1);
        }
    }

    /**
     * Tests bookkeeping of outdated nodes
     */
    @Test
    public void testBookkeeping() {
        for (final boolean f_caching : new boolean[] {true, false})
            for (final boolean f_external : f_caching ? new boolean[] {true, false} : new boolean[] {false})
                for (final boolean g_caching : new boolean[] {true, false})
                    for (final boolean g_external : g_caching ? new boolean[] {true, false} : new boolean[] {false})
                        for (final boolean h_caching : new boolean[] {true, false})
                            for (final boolean h_external : h_caching ? new boolean[] {true, false} : new boolean[] {false})
                                performBookkeepingTest(f_caching, f_external, g_caching, g_external, h_caching, h_external);
    }

    private void performBookkeepingTest(final boolean f_caching, final boolean f_external,
                                        final boolean g_caching, final boolean g_external,
                                        final boolean h_caching, final boolean h_external) {
        for (int i = 0; i < NUM_TRIALS; i++) {
            final ImmutableComputableGraph icg_0 = getTestICGBuilder(f_caching, f_external, g_caching, g_external,
                    h_caching, h_external).build();

            Assert.assertTrue(!icg_0.isValueDirectlyAvailable("x"));
            Assert.assertTrue(!icg_0.isValueDirectlyAvailable("y"));
            Assert.assertTrue(!icg_0.isValueDirectlyAvailable("z"));
            Assert.assertTrue(!icg_0.isValueDirectlyAvailable("f"));
            Assert.assertTrue(!icg_0.isValueDirectlyAvailable("g"));
            Assert.assertTrue(!icg_0.isValueDirectlyAvailable("h"));

            ImmutableComputableGraph icg_tmp = icg_0;
            icg_tmp = icg_tmp.setValue("x", new DuplicableNDArray(getRandomINDArray()));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("x"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("y"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("z"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("f"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("g"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("h"));

            icg_tmp = icg_tmp.setValue("y", new DuplicableNumber<>(getRandomDouble()));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("x"));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("y"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("z"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("f"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("g"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("h"));

            icg_tmp = icg_tmp.setValue("z", new DuplicableNDArray(getRandomINDArray()));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("x"));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("y"));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("z"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("f"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("g"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("h"));

            try {
                icg_tmp = icg_tmp.updateAllCaches();
            } catch (final Exception ex) {
                if (!f_external && !g_external && !h_external) {
                    throw new AssertionError("Could not update all caches but it should have been possible");
                } else {
                    icg_tmp = icg_tmp.updateAllCachesIfPossible(); /* this will not throw exception by design */
                }
            }

            Assert.assertTrue(!f_caching ||
                    (f_external && !icg_tmp.isValueDirectlyAvailable("f")) ||
                    (!f_external && icg_tmp.isValueDirectlyAvailable("f")));

            Assert.assertTrue(!g_caching ||
                    (g_external && !icg_tmp.isValueDirectlyAvailable("g")) ||
                    (!g_external && icg_tmp.isValueDirectlyAvailable("g")));

            if (!f_external && !g_external) {
                Assert.assertTrue(!h_caching ||
                        (h_external && !icg_tmp.isValueDirectlyAvailable("h")) ||
                        (!h_external && icg_tmp.isValueDirectlyAvailable("h")));
            } else {
                Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("h"));
            }

            /* fill in the external values */
            if (f_external) {
                icg_tmp = icg_tmp.setValue("f", f_computation_function.apply(
                        ImmutableMap.of("x", icg_tmp.getValueDirect("x"), "y", icg_tmp.getValueDirect("y"))));
                Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("f"));
            }

            if (g_external) {
                icg_tmp = icg_tmp.setValue("g", g_computation_function.apply(
                        ImmutableMap.of("y", icg_tmp.getValueDirect("y"), "z", icg_tmp.getValueDirect("z"))));
                Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("g"));
            }

            if (h_external) {
                icg_tmp = icg_tmp.setValue("h", h_computation_function.apply(ImmutableMap.of(
                        "f", icg_tmp.getValueWithRequiredEvaluations("f"),
                        "g", icg_tmp.getValueWithRequiredEvaluations("g"),
                        "x", icg_tmp.getValueDirect("x"))));
                Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("h"));
            }

            /* since all externally computed nodes are initialized, a call to updateAllCaches() must succeed */
            icg_tmp = icg_tmp.updateAllCaches();

            /* at this point, every caching node must be up-to-date */
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("x"));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("y"));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("z"));
            Assert.assertTrue(!f_caching || icg_tmp.isValueDirectlyAvailable("f"));
            Assert.assertTrue(!g_caching || icg_tmp.isValueDirectlyAvailable("g"));
            Assert.assertTrue(!h_caching || icg_tmp.isValueDirectlyAvailable("h"));

            /* update x -- f and h must go out of date */
            ImmutableComputableGraph icg_tmp_x = icg_tmp.setValue("x", new DuplicableNDArray(getRandomINDArray()));
            Assert.assertTrue(!icg_tmp_x.isValueDirectlyAvailable("f"));
            Assert.assertTrue(!g_caching || icg_tmp_x.isValueDirectlyAvailable("g"));
            Assert.assertTrue(!icg_tmp_x.isValueDirectlyAvailable("h"));

            /* update y -- f, g and h must go out of date */
            ImmutableComputableGraph icg_tmp_y = icg_tmp.setValue("y", new DuplicableNumber<>(getRandomDouble()));
            Assert.assertTrue(!icg_tmp_y.isValueDirectlyAvailable("f"));
            Assert.assertTrue(!icg_tmp_y.isValueDirectlyAvailable("g"));
            Assert.assertTrue(!icg_tmp_y.isValueDirectlyAvailable("h"));

            /* update z -- g and h must go out of date */
            ImmutableComputableGraph icg_tmp_z = icg_tmp.setValue("z", new DuplicableNDArray(getRandomINDArray()));
            Assert.assertTrue(!f_caching || icg_tmp_z.isValueDirectlyAvailable("f"));
            Assert.assertTrue(!icg_tmp_z.isValueDirectlyAvailable("g"));
            Assert.assertTrue(!icg_tmp_z.isValueDirectlyAvailable("h"));

            /* update x and y -- f, g and h must go out of date */
            ImmutableComputableGraph icg_tmp_xy = icg_tmp
                    .setValue("x", new DuplicableNDArray(getRandomINDArray()))
                    .setValue("y", new DuplicableNumber<>(getRandomDouble()));
            Assert.assertTrue(!icg_tmp_xy.isValueDirectlyAvailable("f"));
            Assert.assertTrue(!icg_tmp_xy.isValueDirectlyAvailable("g"));
            Assert.assertTrue(!icg_tmp_xy.isValueDirectlyAvailable("h"));

            /* update x and z -- f, g and h must go out of date */
            ImmutableComputableGraph icg_tmp_xz = icg_tmp
                    .setValue("x", new DuplicableNDArray(getRandomINDArray()))
                    .setValue("z", new DuplicableNDArray(getRandomINDArray()));
            Assert.assertTrue(!icg_tmp_xz.isValueDirectlyAvailable("f"));
            Assert.assertTrue(!icg_tmp_xz.isValueDirectlyAvailable("g"));
            Assert.assertTrue(!icg_tmp_xz.isValueDirectlyAvailable("h"));

            /* update x and z -- f, g and h must go out of date */
            ImmutableComputableGraph icg_tmp_xyz = icg_tmp
                    .setValue("x", new DuplicableNDArray(getRandomINDArray()))
                    .setValue("y", new DuplicableNumber<>(getRandomDouble()))
                    .setValue("z", new DuplicableNDArray(getRandomINDArray()));
            Assert.assertTrue(!icg_tmp_xyz.isValueDirectlyAvailable("f"));
            Assert.assertTrue(!icg_tmp_xyz.isValueDirectlyAvailable("g"));
            Assert.assertTrue(!icg_tmp_xyz.isValueDirectlyAvailable("h"));

            if (f_external) {
                /* update f -- h must go out of date */
                ImmutableComputableGraph icg_tmp_f = icg_tmp
                        .setValue("f", new DuplicableNDArray(getRandomINDArray()));
                Assert.assertTrue(!g_caching || icg_tmp_f.isValueDirectlyAvailable("g"));
                Assert.assertTrue(!icg_tmp_f.isValueDirectlyAvailable("h"));
            }

            if (g_external) {
                /* update g -- h must go out of date */
                ImmutableComputableGraph icg_tmp_g = icg_tmp
                        .setValue("g", new DuplicableNDArray(getRandomINDArray()));
                Assert.assertTrue(!f_caching || icg_tmp_g.isValueDirectlyAvailable("f"));
                Assert.assertTrue(!icg_tmp_g.isValueDirectlyAvailable("h"));
            }

            if (f_external && g_external) {
                /* update f and g -- h must go out of date */
                ImmutableComputableGraph icg_tmp_fg = icg_tmp
                        .setValue("f", new DuplicableNDArray(getRandomINDArray()))
                        .setValue("g", new DuplicableNDArray(getRandomINDArray()));
                Assert.assertTrue(!icg_tmp_fg.isValueDirectlyAvailable("h"));
            }
        }
    }

    /**
     * Tests propagation of tags from descendents to parents
     */
    @Test
    public void testTagPropagation() {
        final int MAX_TAGS_PER_NODE = 5;
        for (int i = 0; i < NUM_TRIALS; i++) {
            final Set<String> x_tags = getRandomSetOfStrings(rng.nextInt(MAX_TAGS_PER_NODE));
            final Set<String> y_tags = getRandomSetOfStrings(rng.nextInt(MAX_TAGS_PER_NODE));
            final Set<String> z_tags = getRandomSetOfStrings(rng.nextInt(MAX_TAGS_PER_NODE));
            final Set<String> f_tags = getRandomSetOfStrings(rng.nextInt(MAX_TAGS_PER_NODE));
            final Set<String> g_tags = getRandomSetOfStrings(rng.nextInt(MAX_TAGS_PER_NODE));
            final Set<String> h_tags = getRandomSetOfStrings(rng.nextInt(MAX_TAGS_PER_NODE));
            final ImmutableComputableGraph icg = getTestICGBuilder(true, true, true, true, true, true,
                    x_tags.toArray(new String[0]), y_tags.toArray(new String[0]),
                    z_tags.toArray(new String[0]), f_tags.toArray(new String[0]),
                    g_tags.toArray(new String[0]), h_tags.toArray(new String[0])).build();

            final Set<String> all_x_tags = Sets.union(Sets.union(x_tags, f_tags), h_tags);
            final Set<String> all_y_tags = Sets.union(Sets.union(Sets.union(y_tags, f_tags), g_tags), h_tags);
            final Set<String> all_z_tags = Sets.union(Sets.union(z_tags, g_tags), h_tags);
            final Set<String> all_f_tags = Sets.union(f_tags, h_tags);
            final Set<String> all_g_tags = Sets.union(g_tags, h_tags);
            final Set<String> all_h_tags = h_tags;

            final Set<String> all_x_tags_actual = icg.getComputableGraphStructure().getAllTagsForNode("x");
            final Set<String> all_y_tags_actual = icg.getComputableGraphStructure().getAllTagsForNode("y");
            final Set<String> all_z_tags_actual = icg.getComputableGraphStructure().getAllTagsForNode("z");
            final Set<String> all_f_tags_actual = icg.getComputableGraphStructure().getAllTagsForNode("f");
            final Set<String> all_g_tags_actual = icg.getComputableGraphStructure().getAllTagsForNode("g");
            final Set<String> all_h_tags_actual = icg.getComputableGraphStructure().getAllTagsForNode("h");

            Assert.assertTrue(all_x_tags.equals(all_x_tags_actual));
            Assert.assertTrue(all_y_tags.equals(all_y_tags_actual));
            Assert.assertTrue(all_z_tags.equals(all_z_tags_actual));
            Assert.assertTrue(all_f_tags.equals(all_f_tags_actual));
            Assert.assertTrue(all_g_tags.equals(all_g_tags_actual));
            Assert.assertTrue(all_h_tags.equals(all_h_tags_actual));
        }
    }

    private Set<String> getRandomSetOfStrings(final int num) {
        return IntStream.range(0, num)
                .mapToObj(n -> String.valueOf(rng.nextLong()))
                .collect(Collectors.toSet());
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
    public void testUnchangedNodesSameReferenceAfterUpdate() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test(enabled = false)
    public void testUninitializedPrimitiveNode() {
    }

    @Test(enabled = false)
    public void testUninitializedExternallyComputedNode() {
    }

    @Test(enabled = false)
    public void testNonRedundantComputation() {
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

    @Test(expectedExceptions = ComputableGraphStructure.CyclicGraphException.class)
    public void testCyclicGraphException_1() {
        ImmutableComputableGraph.builder()
                .addPrimitiveNode("x", new String[] {}, new DuplicableNDArray())
                .addComputableNode("y", new String[] {}, new String[] {"x", "w"}, null, true) /* cycle */
                .addComputableNode("z", new String[] {}, new String[] {"y"}, null, true)
                .addComputableNode("w", new String[] {}, new String[] {"z"}, null, true)
                .build();
    }

    @Test(expectedExceptions = ComputableGraphStructure.CyclicGraphException.class)
    public void testCyclicGraphException_2() {
        ImmutableComputableGraph.builder()
                .addPrimitiveNode("x", new String[] {}, new DuplicableNDArray())
                .addPrimitiveNode("y", new String[] {}, new DuplicableNDArray())
                .addPrimitiveNode("z", new String[] {}, new DuplicableNDArray())
                .addComputableNode("f", new String[] {}, new String[] {"x", "y", "h"}, null, true) /* cycle */
                .addComputableNode("g", new String[] {}, new String[] {"y", "z"}, null, true)
                .addComputableNode("h", new String[] {}, new String[] {"f", "g"}, null, true)
                .build();
    }


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

}
