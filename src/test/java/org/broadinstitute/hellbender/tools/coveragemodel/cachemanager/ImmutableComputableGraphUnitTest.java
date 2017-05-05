package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import avro.shaded.com.google.common.collect.ImmutableMap;
import avro.shaded.com.google.common.collect.Sets;
import org.apache.commons.lang.RandomStringUtils;
import org.broadinstitute.hellbender.utils.MathObjectAsserts;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import javax.annotation.Nonnull;
import java.util.*;
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
     *      \ f   g
     *       \ \  /
     *        \ \/
     *          h
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

    private static final Set<String> ALL_NODES = new HashSet<>(Arrays.asList("x", "y", "z", "f", "g", "h"));
    private static final Set<String> ALL_PRIMITIVE_NODES = new HashSet<>(Arrays.asList("x", "y", "z"));
    private static final Set<String> ALL_COMPUTABLE_NODES = new HashSet<>(Arrays.asList("f", "g", "h"));

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

        MathObjectAsserts.assertNDArrayEquals(xICG, x);
        Assert.assertEquals(yICG, y, EPSILON);
        MathObjectAsserts.assertNDArrayEquals(zICG, z);
        MathObjectAsserts.assertNDArrayEquals(fICG, fExpected);
        MathObjectAsserts.assertNDArrayEquals(gICG, gExpected);
        MathObjectAsserts.assertNDArrayEquals(hICG, hExpected);
    }

    private boolean assertIntactReferences(@Nonnull final ImmutableComputableGraph original,
                                           @Nonnull final ImmutableComputableGraph other,
                                           @Nonnull final Set<String> unaffectedNodeKeys) {
        final Set<String> affectedNodeKeys = unaffectedNodeKeys.stream()
                .filter(nodeKey -> original.getCacheNode(nodeKey) != other.getCacheNode(nodeKey))
                .collect(Collectors.toSet());
        if (!affectedNodeKeys.isEmpty()) {
            throw new AssertionError("Some of the node references have changed but they were supposed to remain" +
                    " intact: " + affectedNodeKeys.stream().collect(Collectors.joining(", ")));
        }
        return true;
    }

    private boolean assertChangedReferences(@Nonnull final ImmutableComputableGraph original,
                                            @Nonnull final ImmutableComputableGraph other,
                                            @Nonnull final Set<String> affectedNodeKeys) {
        final Set<String> unaffectedNodeKeys = affectedNodeKeys.stream()
                .filter(nodeKey -> original.getCacheNode(nodeKey) == other.getCacheNode(nodeKey))
                .collect(Collectors.toSet());
        if (!unaffectedNodeKeys.isEmpty()) {
            throw new AssertionError("Some of the node references have not changed but they were supposed to change: " +
                    unaffectedNodeKeys.stream().collect(Collectors.joining(", ")));
        }
        return true;
    }

    private boolean assertIntactReferences(@Nonnull final ImmutableComputableGraph original,
                                           @Nonnull final ImmutableComputableGraph other,
                                           @Nonnull final String... unaffectedNodeKeys) {
        return assertIntactReferences(original, other, Arrays.stream(unaffectedNodeKeys).collect(Collectors.toSet()));
    }

    private boolean assertChangedReferences(@Nonnull final ImmutableComputableGraph original,
                                            @Nonnull final ImmutableComputableGraph other,
                                            @Nonnull final String... affectedNodeKeys) {
        return assertChangedReferences(original, other, Arrays.stream(affectedNodeKeys).collect(Collectors.toSet()));
    }

    /**
     * Tests a fully automated auto-updating {@link ImmutableComputableGraph}
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
            ImmutableComputableGraph initializedICG = icg
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

    @DataProvider(name = "allPossibleNodeFlags")
    public Object[][] getAllPossibleNodeFlags() {
        final List<Object[]> data = new ArrayList<>();
        for (final boolean f_caching : new boolean[] {true, false})
            for (final boolean f_external : f_caching ? new boolean[] {true, false} : new boolean[] {false})
                for (final boolean g_caching : new boolean[] {true, false})
                    for (final boolean g_external : g_caching ? new boolean[] {true, false} : new boolean[] {false})
                        for (final boolean h_caching : new boolean[] {true, false})
                            for (final boolean h_external : h_caching ? new boolean[] {true, false} : new boolean[] {false})
                                data.add(new Object[] {f_caching, f_external, g_caching, g_external, h_caching, h_external});
        return data.toArray(new Object[data.size()][6]);
    }

    /**
     * Tests bookkeeping of outdated nodes
     */
    @Test(dataProvider = "allPossibleNodeFlags")
    public void testBookkeeping(final boolean f_caching, final boolean f_external,
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
            assertIntactReferences(icg_0, icg_tmp, "y", "z", "g");

            ImmutableComputableGraph icg_tmp_old = icg_tmp;
            icg_tmp = icg_tmp.setValue("y", new DuplicableNumber<>(getRandomDouble()));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("x"));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("y"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("z"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("f"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("g"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("h"));
            assertIntactReferences(icg_tmp_old, icg_tmp, "x", "z");

            icg_tmp_old = icg_tmp;
            icg_tmp = icg_tmp.setValue("z", new DuplicableNDArray(getRandomINDArray()));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("x"));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("y"));
            Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("z"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("f"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("g"));
            Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("h"));
            assertIntactReferences(icg_tmp_old, icg_tmp, "x", "y", "f");

            icg_tmp_old = icg_tmp;
            try {
                icg_tmp = icg_tmp.updateAllCaches();
            } catch (final Exception ex) {
                if (!f_external && !g_external && !h_external) {
                    throw new AssertionError("Could not update all caches but it should have been possible");
                } else {
                    icg_tmp = icg_tmp.updateAllCachesIfPossible(); /* this will not throw exception by design */
                }
            }
            assertIntactReferences(icg_tmp_old, icg_tmp, "x", "y", "z");

            Assert.assertTrue((!f_caching && assertIntactReferences(icg_tmp_old, icg_tmp, "f")) ||
                    (f_external && !icg_tmp.isValueDirectlyAvailable("f") && assertIntactReferences(icg_tmp_old, icg_tmp, "f")) ||
                    (!f_external && icg_tmp.isValueDirectlyAvailable("f") && assertChangedReferences(icg_tmp_old, icg_tmp, "f")));

            Assert.assertTrue((!g_caching && assertIntactReferences(icg_tmp_old, icg_tmp, "g")) ||
                    (g_external && !icg_tmp.isValueDirectlyAvailable("g") && assertIntactReferences(icg_tmp_old, icg_tmp, "g")) ||
                    (!g_external && icg_tmp.isValueDirectlyAvailable("g") && assertChangedReferences(icg_tmp_old, icg_tmp, "g")));

            if (!f_external && !g_external) {
                Assert.assertTrue((!h_caching && assertIntactReferences(icg_tmp_old, icg_tmp, "h")) ||
                        (h_external && !icg_tmp.isValueDirectlyAvailable("h") && assertIntactReferences(icg_tmp_old, icg_tmp, "h")) ||
                        (!h_external && icg_tmp.isValueDirectlyAvailable("h") && assertChangedReferences(icg_tmp_old, icg_tmp, "h")));
            } else {
                Assert.assertTrue(!icg_tmp.isValueDirectlyAvailable("h") && assertIntactReferences(icg_tmp_old, icg_tmp, "h"));
            }

            /* fill in the external values */
            if (f_external) {
                icg_tmp_old = icg_tmp;
                icg_tmp = icg_tmp.setValue("f", f_computation_function.apply(
                        ImmutableMap.of("x", icg_tmp.getValueDirect("x"), "y", icg_tmp.getValueDirect("y"))));
                Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("f"));
                assertIntactReferences(icg_tmp_old, icg_tmp, "x", "y", "z", "g");
                assertChangedReferences(icg_tmp_old, icg_tmp, "f", "h");
            }

            if (g_external) {
                icg_tmp_old = icg_tmp;
                icg_tmp = icg_tmp.setValue("g", g_computation_function.apply(
                        ImmutableMap.of("y", icg_tmp.getValueDirect("y"), "z", icg_tmp.getValueDirect("z"))));
                Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("g"));
                assertIntactReferences(icg_tmp_old, icg_tmp, "x", "y", "z", "f");
                assertChangedReferences(icg_tmp_old, icg_tmp, "g", "h");
            }

            if (h_external) {
                icg_tmp_old = icg_tmp;
                icg_tmp = icg_tmp.setValue("h", h_computation_function.apply(ImmutableMap.of(
                        "f", icg_tmp.getValueWithRequiredEvaluations("f"),
                        "g", icg_tmp.getValueWithRequiredEvaluations("g"),
                        "x", icg_tmp.getValueDirect("x"))));
                Assert.assertTrue(icg_tmp.isValueDirectlyAvailable("h"));
                assertIntactReferences(icg_tmp_old, icg_tmp, "x", "y", "z", "f", "g");
                assertChangedReferences(icg_tmp_old, icg_tmp, "h");
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
    @Test(dataProvider = "allPossibleNodeFlags")
    public void testTagPropagation(final boolean f_caching, final boolean f_external,
                                   final boolean g_caching, final boolean g_external,
                                   final boolean h_caching, final boolean h_external) {
        for (int i = 0; i < NUM_TRIALS; i++) {
            final Set<String> x_tags = getRandomSetOfTags();
            final Set<String> y_tags = getRandomSetOfTags();
            final Set<String> z_tags = getRandomSetOfTags();
            final Set<String> f_tags = getRandomSetOfTags();
            final Set<String> g_tags = getRandomSetOfTags();
            final Set<String> h_tags = getRandomSetOfTags();
            final ImmutableComputableGraph icg = getTestICGBuilder(
                    f_caching, f_external, g_caching, g_external, h_caching, h_external,
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

    private Set<String> getRandomSetOfTags() {
        final int MAX_NUM_TAGS = 5;
        final int TAG_LENGTH = 12;
        return IntStream.range(0, rng.nextInt(MAX_NUM_TAGS))
                .mapToObj(n -> RandomStringUtils.randomAlphanumeric(TAG_LENGTH))
                .collect(Collectors.toSet());
    }

    private Map<String, INDArray> getExpectedComputableNodeValues(final Duplicable x, final Duplicable y, final Duplicable z) {
        final INDArray xVal = DuplicableNDArray.strip(x);
        final Double yVal = DuplicableNumber.strip(y);
        final INDArray zVal = DuplicableNDArray.strip(z);
        final INDArray fExpected = xVal.add(yVal);
        final INDArray gExpected = zVal.mul(yVal);
        final INDArray hExpected = fExpected.mul(gExpected).sub(xVal);
        return ImmutableMap.of("f", fExpected, "g", gExpected, "h", hExpected);
    }

    /**
     * Tests {@link ImmutableComputableGraph#updateCachesForTag(String)}}
     */
    @Test(dataProvider = "allPossibleNodeFlags")
    public void testUpdateCachesByTag(final boolean f_caching, final boolean f_external,
                                      final boolean g_caching, final boolean g_external,
                                      final boolean h_caching, final boolean h_external) {
        for (int i = 0; i < NUM_TRIALS; i++) {
            final ImmutableComputableGraph icg_empty = getTestICGBuilder(
                    f_caching, f_external, g_caching, g_external, h_caching, h_external,
                    getRandomSetOfTags().toArray(new String[0]), getRandomSetOfTags().toArray(new String[0]),
                    getRandomSetOfTags().toArray(new String[0]), getRandomSetOfTags().toArray(new String[0]),
                    getRandomSetOfTags().toArray(new String[0]), getRandomSetOfTags().toArray(new String[0])).build();

            final Set<String> all_x_tags = icg_empty.getComputableGraphStructure().getAllTagsForNode("x");
            final Set<String> all_y_tags = icg_empty.getComputableGraphStructure().getAllTagsForNode("y");
            final Set<String> all_z_tags = icg_empty.getComputableGraphStructure().getAllTagsForNode("z");
            final Set<String> all_f_tags = icg_empty.getComputableGraphStructure().getAllTagsForNode("f");
            final Set<String> all_g_tags = icg_empty.getComputableGraphStructure().getAllTagsForNode("g");
            final Set<String> all_h_tags = icg_empty.getComputableGraphStructure().getAllTagsForNode("h");
            final Set<String> all_tags = new HashSet<>();
            all_tags.addAll(all_x_tags); all_tags.addAll(all_y_tags); all_tags.addAll(all_z_tags);
            all_tags.addAll(all_f_tags); all_tags.addAll(all_g_tags); all_tags.addAll(all_h_tags);

            final INDArray x = getRandomINDArray();
            final double y = getRandomDouble();
            final INDArray z = getRandomINDArray();
            final ImmutableComputableGraph icg_0 = icg_empty
                    .setValue("x", new DuplicableNDArray(x))
                    .setValue("y", new DuplicableNumber<>(y))
                    .setValue("z", new DuplicableNDArray(z));
            final Map<String, INDArray> expectedComputableNodeValues = getExpectedComputableNodeValues(
                    icg_0.getValueDirect("x"), icg_0.getValueDirect("y"), icg_0.getValueDirect("z"));

            for (final String tag : all_tags) {
                ImmutableComputableGraph icg_1;
                Counter startCounter;
                try {
                    startCounter = getCounterInstance();
                    icg_1 = icg_0.updateCachesForTag(tag);
                } catch (final Exception ex) { /* should fail only if some of the tagged nodes are external */
                    if (!f_external && !g_external && !h_external) {
                        throw new AssertionError("Could not update tagged nodes but it should have been possible");
                    }
                    startCounter = getCounterInstance();
                    icg_1 = icg_0.updateCachesForTagIfPossible(tag);
                }
                final Counter evalCounts = getCounterInstance().diff(startCounter);

                /* check updated caches */
                final Set<String> updatedNodesExpected = new HashSet<>();
                if (!f_external && f_caching && all_f_tags.contains(tag)) {
                    updatedNodesExpected.add("f");
                }
                if (!g_external && g_caching && all_g_tags.contains(tag)) {
                    updatedNodesExpected.add("g");
                }
                if (!h_external && !f_external && !g_external && h_caching && all_h_tags.contains(tag)) {
                    updatedNodesExpected.add("h");
                }
                assertChangedReferences(icg_0, icg_1, updatedNodesExpected);
                assertIntactReferences(icg_0, icg_1, Sets.difference(ALL_NODES, updatedNodesExpected));

                for (final String nodeKey : updatedNodesExpected) {
                    Assert.assertTrue(icg_1.isValueDirectlyAvailable(nodeKey));
                    MathObjectAsserts.assertNDArrayEquals(DuplicableNDArray.strip(icg_1.getValueDirect(nodeKey)),
                            expectedComputableNodeValues.get(nodeKey));
                }
                for (final String nodeKey : Sets.difference(ALL_COMPUTABLE_NODES, updatedNodesExpected)) {
                    Assert.assertTrue(!icg_1.isValueDirectlyAvailable(nodeKey));
                }

                /* check function evaluation counts */
                if ((!f_external && all_f_tags.contains(tag)) /* f is computable and caching */ ||
                        (all_h_tags.contains(tag) && !f_external && !g_external && !h_external) /* h, as a descendant, is computable */) {
                    Assert.assertEquals(evalCounts.getCount("f"), 1);
                } else {
                    Assert.assertEquals(evalCounts.getCount("f"), 0);
                }
                if ((!g_external && all_g_tags.contains(tag)) /* g is computable and caching */ ||
                        (all_h_tags.contains(tag) && !g_external && !f_external && !h_external) /* h, as a descendant, is computable */) {
                    Assert.assertEquals(evalCounts.getCount("g"), 1);
                } else {
                    Assert.assertEquals(evalCounts.getCount("g"), 0);
                }
                if (all_h_tags.contains(tag) && !f_external && !g_external && !h_external) {
                    Assert.assertEquals(evalCounts.getCount("h"), 1);
                } else {
                    Assert.assertEquals(evalCounts.getCount("h"), 0);
                }
            }
        }
    }

    /**
     * Tests {@link ImmutableComputableGraph#updateCachesForNode(String)}}
     */
    @Test(dataProvider = "allPossibleNodeFlags")
    public void testUpdateCacheByNode(final boolean f_caching, final boolean f_external,
                                      final boolean g_caching, final boolean g_external,
                                      final boolean h_caching, final boolean h_external) {
        for (int i = 0; i < NUM_TRIALS; i++) {
            final ImmutableComputableGraph icg_empty = getTestICGBuilder(
                    f_caching, f_external, g_caching, g_external, h_caching, h_external).build();

            final INDArray x = getRandomINDArray();
            final double y = getRandomDouble();
            final INDArray z = getRandomINDArray();
            final ImmutableComputableGraph icg_0 = icg_empty
                    .setValue("x", new DuplicableNDArray(x))
                    .setValue("y", new DuplicableNumber<>(y))
                    .setValue("z", new DuplicableNDArray(z));
            final Map<String, INDArray> expectedComputableNodeValues = getExpectedComputableNodeValues(
                    icg_0.getValueDirect("x"), icg_0.getValueDirect("y"), icg_0.getValueDirect("z"));

            for (final String nodeKey : ALL_PRIMITIVE_NODES) {
                Counter startCounter = getCounterInstance();
                ImmutableComputableGraph icg_1 = icg_0.updateCachesForNode(nodeKey);
                final Counter evalCounts = getCounterInstance().diff(startCounter);
                assertIntactReferences(icg_0, icg_1, ALL_NODES);
                evalCounts.assertZero();
            }

            for (final String nodeKey : ALL_COMPUTABLE_NODES) {
                ImmutableComputableGraph icg_1;
                Counter startCounter;
                try {
                    startCounter = getCounterInstance();
                    icg_1 = icg_0.updateCachesForNode(nodeKey);
                } catch (final Exception ex) { /* should fail only if some of the tagged nodes are external */
                    if (!f_external && !g_external && !h_external) {
                        throw new AssertionError("Could not update tagged nodes but it should have been possible");
                    }
                    startCounter = getCounterInstance();
                    icg_1 = icg_0.updateCachesForNodeIfPossible(nodeKey);
                }
                final Counter evalCounts = getCounterInstance().diff(startCounter);

                final boolean isExternal = icg_1.getCacheNode(nodeKey).isExternallyComputed();
                final boolean isCaching = ((ComputableCacheNode)icg_1.getCacheNode(nodeKey)).isCaching();

                switch (nodeKey) {
                    case "f":
                        assertIntactReferences(icg_0, icg_1, "x", "y", "z", "g");
                        if (isExternal) {
                            assertIntactReferences(icg_0, icg_1, ALL_NODES);
                            evalCounts.assertZero();
                        } else {
                            Assert.assertEquals(evalCounts.getCount("f"), 1);
                            Assert.assertEquals(evalCounts.getCount("g"), 0);
                            Assert.assertEquals(evalCounts.getCount("h"), 0);
                            if (isCaching) {
                                Assert.assertTrue(icg_1.isValueDirectlyAvailable("f"));
                                MathObjectAsserts.assertNDArrayEquals(DuplicableNDArray.strip(icg_1.getValueDirect("f")),
                                        expectedComputableNodeValues.get("f"));
                            } else {
                                final Counter before = getCounterInstance();
                                Assert.assertTrue(!icg_1.isValueDirectlyAvailable("f"));
                                MathObjectAsserts.assertNDArrayEquals(DuplicableNDArray.strip(icg_1.getValueWithRequiredEvaluations("f")),
                                        expectedComputableNodeValues.get("f"));
                                final Counter diff = getCounterInstance().diff(before);
                                Assert.assertEquals(diff.getCount("f"), 1);
                                Assert.assertEquals(diff.getCount("g"), 0);
                                Assert.assertEquals(diff.getCount("h"), 0);
                            }
                        }
                        break;

                    case "g":
                        assertIntactReferences(icg_0, icg_1, "x", "y", "z", "f");
                        if (isExternal) {
                            assertIntactReferences(icg_0, icg_1, ALL_NODES);
                            evalCounts.assertZero();
                        } else {
                            Assert.assertEquals(evalCounts.getCount("f"), 0);
                            Assert.assertEquals(evalCounts.getCount("g"), 1);
                            Assert.assertEquals(evalCounts.getCount("h"), 0);
                            if (isCaching) {
                                Assert.assertTrue(icg_1.isValueDirectlyAvailable("g"));
                                MathObjectAsserts.assertNDArrayEquals(DuplicableNDArray.strip(icg_1.getValueDirect("g")),
                                        expectedComputableNodeValues.get("g"));
                            } else {
                                Assert.assertTrue(!icg_1.isValueDirectlyAvailable("g"));
                                final Counter before = getCounterInstance();
                                MathObjectAsserts.assertNDArrayEquals(DuplicableNDArray.strip(icg_1.getValueWithRequiredEvaluations("g")),
                                        expectedComputableNodeValues.get("g"));
                                final Counter diff = getCounterInstance().diff(before);
                                Assert.assertEquals(diff.getCount("f"), 0);
                                Assert.assertEquals(diff.getCount("g"), 1);
                                Assert.assertEquals(diff.getCount("h"), 0);
                            }
                        }
                        break;

                    case "h": /* we just check important cases here */
                        assertIntactReferences(icg_0, icg_1, "x", "y", "z");
                        if (h_external && f_external && g_external) {
                            assertIntactReferences(icg_0, icg_1, ALL_NODES);
                            evalCounts.assertZero();
                        } else if (!h_external && !f_external && !g_external) {
                            Assert.assertEquals(evalCounts.getCount("f"), 1);
                            Assert.assertEquals(evalCounts.getCount("g"), 1);
                            Assert.assertEquals(evalCounts.getCount("h"), 1);
                            if (h_caching) {
                                Assert.assertTrue(icg_1.isValueDirectlyAvailable("h"));
                                MathObjectAsserts.assertNDArrayEquals(
                                        DuplicableNDArray.strip(icg_1.getValueDirect("h")),
                                        expectedComputableNodeValues.get("h"));
                            } else {
                                Assert.assertTrue(!icg_1.isValueDirectlyAvailable("h"));
                                final Counter before = getCounterInstance();
                                MathObjectAsserts.assertNDArrayEquals(
                                        DuplicableNDArray.strip(icg_1.getValueWithRequiredEvaluations("h")),
                                        expectedComputableNodeValues.get("h"));
                                final Counter diff = getCounterInstance().diff(before);
                                Assert.assertEquals(diff.getCount("f"), f_caching ? 0 : 1);
                                Assert.assertEquals(diff.getCount("g"), g_caching ? 0 : 1);
                                Assert.assertEquals(diff.getCount("h"), 1);
                            }
                        }
                        break;
                }
            }
        }
    }

    @Test
    public void testUninitializedPrimitiveNode() {
        final ImmutableComputableGraph icg = getTestICGBuilder(true, false, true, false, true, false).build()
                .setValue("x", new DuplicableNDArray(getRandomINDArray()))
                .setValue("y", new DuplicableNumber<>(getRandomDouble()));
        boolean failed = false;
        try {
            icg.updateAllCaches();
        } catch (final PrimitiveCacheNode.PrimitiveValueNotInitializedException ex) {
            failed = true;
        }
        if (!failed) {
            throw new AssertionError("Expected PrimitiveValueNotInitializedException but it was not thrown");
        }

        icg.updateCachesForNode("f"); /* should not fail */

        failed = false;
        try {
            icg.updateCachesForNode("g");
        } catch (final PrimitiveCacheNode.PrimitiveValueNotInitializedException ex) {
            failed = true;
        }
        if (!failed) {
            throw new AssertionError("Expected PrimitiveValueNotInitializedException but it was not thrown");
        }

        failed = false;
        try {
            icg.updateCachesForNode("h");
        } catch (final PrimitiveCacheNode.PrimitiveValueNotInitializedException ex) {
            failed = true;
        }
        if (!failed) {
            throw new AssertionError("Expected PrimitiveValueNotInitializedException but it was not thrown");
        }
    }

    @Test
    public void testUninitializedExternallyComputedNode() {
        final ImmutableComputableGraph icg = getTestICGBuilder(true, true, true, false, true, false).build()
                .setValue("x", new DuplicableNDArray(getRandomINDArray()))
                .setValue("y", new DuplicableNumber<>(getRandomDouble()))
                .setValue("z", new DuplicableNDArray(getRandomINDArray()));
        boolean failed = false;
        try {
            icg.updateAllCaches();
        } catch (final ComputableCacheNode.ExternallyComputableNodeValueUnavailableException ex) {
            failed = true;
        }
        if (!failed) {
            throw new AssertionError("Expected ExternallyComputableNodeValueUnavailableException but it was not thrown");
        }

        icg.updateCachesForNode("g"); /* should not fail */

        failed = false;
        try {
            icg.updateCachesForNode("f");
        } catch (final ComputableCacheNode.ExternallyComputableNodeValueUnavailableException ex) {
            failed = true;
        }
        if (!failed) {
            throw new AssertionError("Expected ExternallyComputableNodeValueUnavailableException but it was not thrown");
        }

        failed = false;
        try {
            icg.updateCachesForNode("h");
        } catch (final ComputableCacheNode.ExternallyComputableNodeValueUnavailableException ex) {
            failed = true;
        }
        if (!failed) {
            throw new AssertionError("Expected ExternallyComputableNodeValueUnavailableException but it was not thrown");
        }

        /* supply f */
        ImmutableComputableGraph icg_1 = icg.setValue("f", f_computation_function.apply(
                ImmutableMap.of("x", icg.getValueDirect("x"), "y", icg.getValueDirect("y"))));
        Assert.assertTrue(icg_1.isValueDirectlyAvailable("f"));

        /* cache g */
        Assert.assertTrue(!icg_1.isValueDirectlyAvailable("g"));
        Counter before = getCounterInstance();
        icg_1 = icg_1.updateCachesForNode("g");
        Assert.assertTrue(icg_1.isValueDirectlyAvailable("g"));
        Counter diff = getCounterInstance().diff(before);
        Assert.assertEquals(diff.getCount("f"), 0);
        Assert.assertEquals(diff.getCount("g"), 1);
        Assert.assertEquals(diff.getCount("h"), 0);

        /* cache h -- now, it is computable */
        Assert.assertTrue(!icg_1.isValueDirectlyAvailable("h"));
        before = getCounterInstance();
        icg_1 = icg_1.updateCachesForNode("h");
        Assert.assertTrue(icg_1.isValueDirectlyAvailable("h"));
        diff = getCounterInstance().diff(before);
        Assert.assertEquals(diff.getCount("f"), 0);
        Assert.assertEquals(diff.getCount("g"), 0);
        Assert.assertEquals(diff.getCount("h"), 1);

        /* updating all caches must have no effect */
        before = getCounterInstance();
        ImmutableComputableGraph icg_2 = icg_1.updateAllCaches();
        getCounterInstance().diff(before).assertZero();
        Assert.assertTrue(icg_2.isValueDirectlyAvailable("f"));
        Assert.assertTrue(icg_2.isValueDirectlyAvailable("g"));
        Assert.assertTrue(icg_2.isValueDirectlyAvailable("h"));
        assertIntactReferences(icg_1, icg_2, ALL_NODES);
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

        public void assertZero() {
            Assert.assertTrue(counts.values().stream().allMatch(val -> val == 0));
        }
    }

}
