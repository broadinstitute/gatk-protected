package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import avro.shaded.com.google.common.collect.ImmutableMap;
import avro.shaded.com.google.common.collect.Sets;
import org.apache.commons.lang.RandomStringUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link ComputableGraphStructure}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class ComputableGraphStructureUnitTest extends BaseTest {

    private static final Random rng = new Random(1984);
    private static final int MAX_DAG_DEPTH = 10;
    private static final int MAX_NODES_PER_LAYER = 10;
    private static final int MAX_TAGS_PER_NODE = 10;
    private static final int MAX_PARENTS_PER_NODE = 10;
    private static final int NUM_TRIALS = 10;

    @Test(expectedExceptions = ComputableGraphStructure.NonexistentParentNodeKey.class)
    public void testMissingParents() {
        ImmutableComputableGraph.builder()
                .primitiveNode("x", new String[] {}, new DuplicableNDArray())
                .primitiveNode("y", new String[] {}, new DuplicableNumber<Double>())
                .primitiveNode("z", new String[] {}, new DuplicableNDArray())
                .computableNode("f", new String[] {}, new String[] {"x", "y", "q"}, null, true) /* q is undefined */
                .computableNode("g", new String[] {}, new String[] {"y", "z"}, null, true)
                .computableNode("h", new String[] {}, new String[] {"f", "g"}, null, true)
                .build();
    }

    @Test(expectedExceptions = ComputableGraphStructure.CyclicGraphException.class)
    public void testCyclicGraphException_1() {
        ImmutableComputableGraph.builder()
                .primitiveNode("x", new String[] {}, new DuplicableNDArray())
                .computableNode("y", new String[] {}, new String[] {"x", "w"}, null, true) /* cycle */
                .computableNode("z", new String[] {}, new String[] {"y"}, null, true)
                .computableNode("w", new String[] {}, new String[] {"z"}, null, true)
                .build();
    }

    @Test(expectedExceptions = ComputableGraphStructure.CyclicGraphException.class)
    public void testCyclicGraphException_2() {
        ImmutableComputableGraph.builder()
                .primitiveNode("x", new String[] {}, new DuplicableNDArray())
                .primitiveNode("y", new String[] {}, new DuplicableNDArray())
                .primitiveNode("z", new String[] {}, new DuplicableNDArray())
                .computableNode("f", new String[] {}, new String[] {"x", "y", "h"}, null, true) /* cycle */
                .computableNode("g", new String[] {}, new String[] {"y", "z"}, null, true)
                .computableNode("h", new String[] {}, new String[] {"f", "g"}, null, true)
                .build();
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testNodeTagsAndKeys() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        final Set<String> cgsNodeTagsSet = cgs.getNodeTagsSet();
        final Set<String> dagNodeTagsSet = dag.tagsSet;
        final Set<String> cgsNodeKeysSet = cgs.getNodeKeysSet();
        final Set<String> dagNodeKeysSet = dag.nodeKeysSet;
        Assert.assertTrue(cgsNodeKeysSet.equals(dagNodeKeysSet));
        Assert.assertTrue(cgsNodeTagsSet.equals(dagNodeTagsSet));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testTopologicalOrder() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(dag.topologicalOrderMap.get(nodeKey) ==
                cgs.getTopologicalOrder(nodeKey)));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testAncestors() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(dag.getAncestors(nodeKey).equals(cgs.getAncestors(nodeKey))));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testDescendants() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(dag.getDescendents(nodeKey).equals(cgs.getDescendants(nodeKey))));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testParents() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(dag.getParents(nodeKey).equals(cgs.getParents(nodeKey))));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testChildren() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(dag.getChildren(nodeKey).equals(cgs.getChildren(nodeKey))));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testInducedTags() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(dag.getInducedTags(nodeKey).equals(cgs.getInducedTagsForNode(nodeKey))));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testTopologicalOrderForNodeEvaluation() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(isTopologicallyEquivalent(dag.getTopologicalOrderForNodeEvaluation(nodeKey),
                cgs.getTopologicalOrderForNodeEvaluation(nodeKey), dag.topologicalOrderMap)));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testTopologicalOrderForNodeMutation() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.nodeKeysSet.forEach(nodeKey -> Assert.assertTrue(isTopologicallyEquivalent(dag.getTopologicalOrderForNodeMutation(nodeKey),
                cgs.getTopologicalOrderForNodeMutation(nodeKey), dag.topologicalOrderMap)));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testTopologicalOrderForTagEvaluation() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        dag.tagsSet.forEach(tag -> Assert.assertTrue(isTopologicallyEquivalent(dag.getTopologicalOrderForTagEvaluation(tag),
                cgs.getTopologicalOrderForTagEvaluation(tag), dag.topologicalOrderMap)));
    }

    @Test(invocationCount = NUM_TRIALS)
    public void testTopologicalOrderForCompleteEvaluation() {
        final RandomDAG dag = RandomDAG.getRandomDAG();
        final ComputableGraphStructure cgs = new ComputableGraphStructure(dag.getEquivalentCacheNodeSet());
        final List<String> orderedNodes = new ArrayList<>(dag.nodeKeysSet);
        orderedNodes.sort(Comparator.comparingInt(dag.topologicalOrderMap::get));
        Assert.assertTrue(isTopologicallyEquivalent(cgs.getTopologicalOrderForCompleteEvaluation(), orderedNodes,
                dag.topologicalOrderMap));
    }

    @Test
    public void testAssertTopologicallyEquivalentLists() {
        final Map<String, Integer> topologicalOrderMap = ImmutableMap.<String, Integer>builder()
                .put("a0", 0).put("b0", 0).put("c0", 0)
                .put("a1", 1).put("b1", 1).put("c1", 1)
                .put("a2", 2).put("b2", 2).put("c2", 2).put("d2", 2).build();
        Assert.assertTrue(isTopologicallyEquivalent(
                Arrays.asList("a0", "c0", "a2", "d2"),
                Arrays.asList("c0", "a0", "a2", "d2"),
                topologicalOrderMap));
        Assert.assertTrue(isTopologicallyEquivalent(
                Arrays.asList("a0", "c0", "b1", "a1", "c1", "a2", "d2"),
                Arrays.asList("c0", "a0", "b1", "c1", "a1", "a2", "d2"),
                topologicalOrderMap));
        Assert.assertTrue(!isTopologicallyEquivalent(
                Arrays.asList("c0", "a2", "d2"),
                Arrays.asList("c0", "a0", "a2", "d2"),
                topologicalOrderMap));
        Assert.assertTrue(!isTopologicallyEquivalent(
                Arrays.asList("a0", "c0", "b1", "a1", "c1", "a2", "c2"),
                Arrays.asList("c0", "a0", "b1", "c1", "a1", "a2", "d2"),
                topologicalOrderMap));
    }

    /**
     * This test helper class creates a random DAG starting from a topological order. All helper methods are
     * implemented in a brute-force manner.
     */
    private static final class RandomDAG {
        private static final int TAG_LENGTH = 32;
        private static final int NODE_KEY_LENGTH = 32;

        final Map<Integer, Set<String>> nodesByTopologicalOrder;
        final Map<String, Integer> topologicalOrderMap;
        final Map<String, Set<String>> parentsMap;
        final Map<String, Set<String>> tagsMap;
        final Set<String> nodeKeysSet;
        final Set<String> tagsSet;

        private RandomDAG(final int depth, final int maxNodesPerLayer, final int maxParentsPerNode, final int maxTagsPerNode) {
            Utils.validateArg(depth >= 0, "DAG depth must be  >= 0");
            Utils.validateArg(maxNodesPerLayer > 0, "Max nodes per layer must be positive");
            Utils.validateArg(maxParentsPerNode > 0, "Max parents per node must be positive");
            Utils.validateArg(maxTagsPerNode > 0, "Max tags per node must be positive");

            nodesByTopologicalOrder = new HashMap<>();
            parentsMap = new HashMap<>();
            tagsMap = new HashMap<>();
            nodeKeysSet = new HashSet<>();

            for (int d = 0; d <= depth; d++) {
                final Set<String> randomNodes = getRandomNodeKeys(maxNodesPerLayer);
                nodeKeysSet.addAll(randomNodes);
                nodesByTopologicalOrder.put(d, randomNodes);
                randomNodes.forEach(nodeKey -> tagsMap.put(nodeKey, getRandomTags(maxTagsPerNode)));
                randomNodes.forEach(nodeKey -> parentsMap.put(nodeKey, new HashSet<>()));
                if (d > 0) {
                    final Set<String> possibleAncestors = IntStream.range(0, d)
                            .mapToObj(nodesByTopologicalOrder::get)
                            .flatMap(Set::stream)
                            .collect(Collectors.toSet());
                    for (final String nodeKey : randomNodes) {
                        final String randomParent = getRandomElement(nodesByTopologicalOrder.get(d - 1));
                        final int numParents = rng.nextInt(maxParentsPerNode);
                        final Set<String> randomAncestors = IntStream.range(0, numParents)
                                .mapToObj(i -> getRandomElement(possibleAncestors))
                                .collect(Collectors.toSet());
                        parentsMap.put(nodeKey, new HashSet<>());
                        parentsMap.get(nodeKey).add(randomParent);
                        parentsMap.get(nodeKey).addAll(randomAncestors);
                    }
                }
            }
            topologicalOrderMap = new HashMap<>();
            nodesByTopologicalOrder.entrySet().forEach(entry -> entry.getValue()
                    .forEach(nodeKey -> topologicalOrderMap.put(nodeKey, entry.getKey())));
            tagsSet = tagsMap.values().stream().flatMap(Set::stream).collect(Collectors.toSet());
        }

        static RandomDAG getRandomDAG() {
            return new RandomDAG(rng.nextInt(MAX_DAG_DEPTH),
                    1 + rng.nextInt(MAX_NODES_PER_LAYER),
                    1 + rng.nextInt(MAX_PARENTS_PER_NODE),
                    1 + rng.nextInt(MAX_TAGS_PER_NODE));
        }

        Set<String> getParents(final String nodeKey) {
            return parentsMap.get(nodeKey);
        }

        /**
         * Brute-force method
         */
        Set<String> getAncestors(final String nodeKey) {
            final Set<String> parents = getParents(nodeKey);
            return Sets.union(parents, parents.stream()
                    .map(this::getAncestors)
                    .flatMap(Set::stream)
                    .collect(Collectors.toSet()));
        }

        /**
         * Brute-force method
         */
        Set<String> getChildren(final String nodeKey) {
            return nodeKeysSet.stream()
                    .filter(key -> getParents(key).contains(nodeKey))
                    .collect(Collectors.toSet());
        }

        /**
         * Brute-force method
         */
        Set<String> getDescendents(final String nodeKey) {
            final Set<String> children = getChildren(nodeKey);
            return Sets.union(children, children.stream()
                    .map(this::getDescendents)
                    .flatMap(Set::stream)
                    .collect(Collectors.toSet()));
        }

        Set<String> getInducedTags(final String nodeKey) {
            final Set<String> allNodes = Sets.union(getDescendents(nodeKey), Collections.singleton(nodeKey));
            return allNodes.stream()
                    .map(tagsMap::get)
                    .flatMap(Set::stream)
                    .collect(Collectors.toSet());
        }

        List<String> getTopologicalOrderForNodeEvaluation(final String nodeKey) {
            final List<String> sortedNodes = new ArrayList<>(Sets.union(getAncestors(nodeKey),
                    Collections.singleton(nodeKey)));
            sortedNodes.sort(Comparator.comparingInt(topologicalOrderMap::get));
            return sortedNodes;
        }

        /**
         * Brute-force method
         */
        List<String> getTopologicalOrderForNodeMutation(final String nodeKey) {
            final Set<String> mutatedNodeAndDescendants = Sets.union(getDescendents(nodeKey),
                    Collections.singleton(nodeKey));
            final Set<String> ancestorsOfDescendants = getDescendents(nodeKey)
                    .stream()
                    .map(this::getAncestors)
                    .flatMap(Set::stream)
                    .collect(Collectors.toSet());
            final Set<String> involvedNodes = Sets.union(mutatedNodeAndDescendants, ancestorsOfDescendants);
            final List<String> topologicallySortedInvolvedNodes = new ArrayList<>(involvedNodes);
            topologicallySortedInvolvedNodes.sort(Comparator.comparingInt(topologicalOrderMap::get));
            return topologicallySortedInvolvedNodes;
        }

        /**
         * Brute-force method
         */
        List<String> getTopologicalOrderForTagEvaluation(final String tag) {
            final Set<String> taggedNodes = nodeKeysSet.stream()
                    .filter(nodeKey -> getInducedTags(nodeKey).contains(tag))
                    .collect(Collectors.toSet());
            final Set<String> taggedNodesAndTheirAncestors = Sets.union(taggedNodes,
                    taggedNodes.stream()
                            .map(this::getAncestors)
                            .flatMap(Set::stream)
                            .collect(Collectors.toSet()));
            final List<String> topologicallySortedtaggedNodesAndTheirAncestors = new ArrayList<>(taggedNodesAndTheirAncestors);
            topologicallySortedtaggedNodesAndTheirAncestors.sort(Comparator.comparingInt(topologicalOrderMap::get));
            return topologicallySortedtaggedNodesAndTheirAncestors;
        }

        Set<CacheNode> getEquivalentCacheNodeSet() {
            final List<String> shuffledNodeList = new ArrayList<>(nodeKeysSet);
            Collections.shuffle(shuffledNodeList, rng);
            return shuffledNodeList.stream()
                    .map(nodeKey -> topologicalOrderMap.get(nodeKey) == 0
                            ? new PrimitiveCacheNode(nodeKey, tagsMap.get(nodeKey), null)
                            : new ComputableCacheNode(nodeKey, tagsMap.get(nodeKey), getParents(nodeKey), null, true))
                    .collect(Collectors.toSet());
        }

        private static Set<String> getRandomNodeKeys(final int maxNodes) {
            return IntStream.range(0, rng.nextInt(maxNodes) + 1)
                    .mapToObj(i -> "NODE_KEY_" + RandomStringUtils.randomAlphanumeric(NODE_KEY_LENGTH))
                    .collect(Collectors.toSet());
        }

        private static Set<String> getRandomTags(final int maxTags) {
            return IntStream.range(0, rng.nextInt(maxTags))
                    .mapToObj(i -> "TAG_" + RandomStringUtils.randomAlphanumeric(TAG_LENGTH))
                    .collect(Collectors.toSet());
        }

        private static String getRandomElement(final Set<String> set) {
            final List<String> list = new ArrayList<>(set);
            return list.get(rng.nextInt(list.size()));
        }
    }

    private boolean isTopologicallyEquivalent(final List<String> actual, final List<String> expected,
                                              final Map<String, Integer> topologicalOrderMap) {
        if (actual == null || expected == null || !(new HashSet<>(actual).equals(new HashSet<>(expected)))) {
            return false;
        }
        Utils.validateArg(topologicalOrderMap.keySet().containsAll(actual), "Some strings have unknown topological order");
        return IntStream.range(0, expected.size()).allMatch(i ->
                topologicalOrderMap.get(actual.get(i)).equals(topologicalOrderMap.get(expected.get(i))));
    }
}
