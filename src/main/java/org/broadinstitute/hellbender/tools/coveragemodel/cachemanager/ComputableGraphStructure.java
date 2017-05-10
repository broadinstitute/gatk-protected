package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import avro.shaded.com.google.common.collect.Sets;
import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class performs consistency checks and computes several structural properties for the DAG specified by a
 * set of {@link CacheNode}s. These include:
 *
 * - assertion for existence of no cycles
 * - construction of maps of nodes to their descendants and ancestors
 * - propagation of tags from descendants to ancestors
 * - topological order for evaluating a computable node
 * - topological order for mutating a primitive/externally-computed node
 * - topological order for evaluating all nodes associated to a tag (see {@link ImmutableComputableGraph})
 * - topological order for complete computation of the graph
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class ComputableGraphStructure implements Serializable {

    private static final long serialVersionUID = -3124293279477371159L;

    private final Set<String> nodeKeysSet;
    private final Set<String> nodeTagsSet;
    private final Map<String, Set<String>> inducedTagsMap;
    private final Map<String, Set<String>> descendantsMap;
    private final Map<String, Set<String>> childrenMap;
    private final Map<String, Set<String>> parentsMap;
    private final Map<String, Set<String>> ancestorsMap;
    private final Map<String, Integer> topologicalOrderMap;
    private final Map<String, List<String>> topologicalOrderForNodeEvaluation;
    private final Map<String, List<String>> topologicalOrderForNodeMutation;
    private final Map<String, List<String>> topologicalOrderForTagEvaluation;
    private final List<String> topologicalOrderForCompleteEvaluation;

    /**
     * An arbitrary negative number to denote the to-be-determined topological order a node
     */
    private static final int UNDEFINED_TOPOLOGICAL_ORDER = -1;

    /**
     * Package-private constructor from a set of {@link CacheNode}s. The graph is specified by the immediate
     * parents and descendants of each node. A {@link CyclicGraphException} is thrown of the graph has a cycle.
     *
     * @param nodeSet a set of {@link CacheNode}s
     */
    ComputableGraphStructure(@Nonnull final Set<CacheNode> nodeSet) {

        /* create the set of keys */
        Utils.nonNull(nodeSet, "The given set of nodes must be non-null");
        nodeKeysSet = extractKeys(nodeSet);
        nodeTagsSet = extractTags(nodeSet);
        assertParentKeysExist(nodeSet, nodeKeysSet);

        /* nodeKey -> set of parents' keys */
        parentsMap = getParentsMap(nodeSet);

        /* nodeKey -> set of children's keys */
        childrenMap = getChildrenMap(nodeSet);

        /* nodeKey -> topological order */
        topologicalOrderMap = getTopologicalOrderMap(parentsMap, nodeKeysSet);

        /* topological order -> set of node keys */
        final Map<Integer, Set<String>> nodesByTopologicalOrderMap = getNodesByTopologicalOrderMap(topologicalOrderMap,
                nodeKeysSet);

        /* nodeKey -> set of descendants' keys */
        descendantsMap = getDescendantsMap(childrenMap, nodesByTopologicalOrderMap);

        /* nodeKey -> set of ancestors' keys */
        ancestorsMap = getAncestorsMap(parentsMap, nodesByTopologicalOrderMap);

        /* nodeKey -> set of upward-propagated tags */
        inducedTagsMap = getInducedTagsMap(nodeSet, nodesByTopologicalOrderMap, ancestorsMap);

        /* tag -> set of tagged nodes, including upward-propagation */
        final Map<String, Set<String>> nodesByInducedTagMap = getNodesByInducedTagMap(inducedTagsMap, nodeKeysSet,
                nodeTagsSet);

        /* topological order for evaluating a single node */
        topologicalOrderForNodeEvaluation = getTopologicalOrderForNodeEvaluation(nodeKeysSet, topologicalOrderMap,
                ancestorsMap);

        /* topological order for evaluating all nodes associated to a tag */
        topologicalOrderForTagEvaluation = getTopologicalOrderForTagEvaluation(nodeTagsSet, topologicalOrderMap,
                nodesByInducedTagMap, ancestorsMap);

        /* topological order for evaluating all nodes */
        topologicalOrderForCompleteEvaluation = getTopologicalOrderForCompleteEvaluation(nodeKeysSet,
                topologicalOrderMap);

        /* topological order for updating the descendants of a mutated node */
        topologicalOrderForNodeMutation = getTopologicalOrderForNodeMutation(nodeKeysSet, ancestorsMap,
                descendantsMap, topologicalOrderMap);
    }

    private static void assertParentKeysExist(@Nonnull final Set<CacheNode> nodeSet,
                                              @Nonnull final Set<String> nodeKeysSet) {
        for (final CacheNode node : nodeSet) {
            Utils.validateArg(nodeKeysSet.containsAll(node.getParents()), () -> {
                final Set<String> undefinedParents = Sets.difference(new HashSet<>(node.getParents()), nodeKeysSet);
                return "Node " + ImmutableComputableGraphUtils.quote(node.getKey()) + " depends on undefined parent(s): " +
                        undefinedParents.stream().map(ImmutableComputableGraphUtils::quote).collect(Collectors.joining(", "));
            });
        }
    }

    private static Set<String> extractTags(@Nonnull Set<CacheNode> nodeSet) {
        return nodeSet.stream().map(CacheNode::getTags).flatMap(Collection::stream).collect(Collectors.toSet());
    }

    private static Set<String> extractKeys(@Nonnull Set<CacheNode> nodeSet) {
        return nodeSet.stream().map(CacheNode::getKey).collect(Collectors.toSet());
    }

    private static Map<String, Integer> getTopologicalOrderMap(@Nonnull final Map<String,Set<String>> immediateParentsMap,
                                                               @Nonnull final Set<String> nodeKeysSet) {
        final Map<String, Integer> topologicalOrderMap = new HashMap<>();
        nodeKeysSet.forEach(key -> topologicalOrderMap.put(key, UNDEFINED_TOPOLOGICAL_ORDER));
        nodeKeysSet.forEach(nodeKey -> updateDepth(nodeKey, 0, nodeKeysSet, immediateParentsMap, topologicalOrderMap));
        return topologicalOrderMap;
    }

    private static Map<Integer, Set<String>> getNodesByTopologicalOrderMap(Map<String, Integer> topologicalOrderMap,
                                                                          Set<String> nodeKeysSet) {
        final int maxDepth = Collections.max(topologicalOrderMap.values());
        final Map<Integer, Set<String>> nodesByTopologicalOrderMap = new HashMap<>();
        IntStream.range(0, maxDepth + 1).forEach(depth ->
                nodesByTopologicalOrderMap.put(depth,
                        nodeKeysSet.stream().filter(node ->
                                topologicalOrderMap.get(node) == depth).collect(Collectors.toSet())));
        return nodesByTopologicalOrderMap;
    }

    private static Map<String, Set<String>> getParentsMap(@Nonnull final Set<CacheNode> nodeSet) {
        return nodeSet.stream()
                .collect(Collectors.toMap(CacheNode::getKey, node -> new HashSet<>(node.getParents())));
    }

    private static Map<String, Set<String>> getChildrenMap(@Nonnull final Set<CacheNode> nodeSet) {
        final Map<String, Set<String>> childrenMap = nodeSet.stream()
                .collect(Collectors.toMap(CacheNode::getKey, node -> new HashSet<String>()));
        nodeSet.forEach(node -> node.getParents().forEach(parentKey -> childrenMap.get(parentKey)
                .add(node.getKey())));
        return childrenMap;
    }

    private static Map<String, Set<String>> getInducedTagsMap(@Nonnull final Set<CacheNode> nodeSet,
                                                              @Nonnull final Map<Integer, Set<String>> nodesByTopologicalOrderMap,
                                                              @Nonnull final Map<String, Set<String>> allParentsMap) {
        final int maxDepth = Collections.max(nodesByTopologicalOrderMap.keySet());
        /* initialize with given tags */
        final Map<String, Set<String>> allTagsMap = nodeSet.stream()
                .collect(Collectors.toMap(CacheNode::getKey, node -> new HashSet<>(node.getTags())));
        /* propagate tags to all parents */
        for (int depth = maxDepth; depth >= 0; depth--) {
            nodesByTopologicalOrderMap.get(depth)
                    .forEach(nodeKey -> allParentsMap.get(nodeKey)
                            .forEach(parentKey -> allTagsMap.get(parentKey).addAll(allTagsMap.get(nodeKey))));
        }
        return allTagsMap;
    }

    private static Map<String, Set<String>> getNodesByInducedTagMap(@Nonnull final Map<String, Set<String>> allTagsMap,
                                                                    @Nonnull final Set<String> nodeKeysSet,
                                                                    @Nonnull final Set<String> nodeTagsSet) {
        final Map<String, Set<String>> nodesByTagMap = nodeTagsSet.stream()
                .collect(Collectors.toMap(Function.identity(), tag -> new HashSet<String>()));
        nodeKeysSet.forEach(nodeKey ->
                allTagsMap.get(nodeKey).forEach(tag ->
                        nodesByTagMap.get(tag).add(nodeKey)));
        return nodesByTagMap;
    }

    private static Map<String, Set<String>> getDescendantsMap(@Nonnull final Map<String, Set<String>> childrenMap,
                                                              @Nonnull final Map<Integer, Set<String>> nodesByTopologicalOrderMap) {
        final Map<String, Set<String>> descendantsMap = new HashMap<>();
        final int maxDepth = Collections.max(nodesByTopologicalOrderMap.keySet());
        /* deepest nodes have no descendants */
        nodesByTopologicalOrderMap.get(maxDepth).forEach(node -> descendantsMap.put(node, new HashSet<>()));
        /* get all descendants by ascending the tree */
        for (int depth = maxDepth - 1; depth >= 0; depth -= 1) {
            for (final String node : nodesByTopologicalOrderMap.get(depth)) {
                final Set<String> nodeDescendants = new HashSet<>();
                nodeDescendants.addAll(childrenMap.get(node));
                for (final String child : childrenMap.get(node)) {
                    nodeDescendants.addAll(descendantsMap.get(child));
                }
                descendantsMap.put(node, nodeDescendants);
            }
        }
        return descendantsMap;
    }

    private static Map<String, Set<String>> getAncestorsMap(@Nonnull final Map<String, Set<String>> parentsMap,
                                                            @Nonnull final Map<Integer, Set<String>> nodesByTopologicalOrderMap) {
        final Map<String, Set<String>> ancestorsMap = new HashMap<>();
        final int maxDepth = Collections.max(nodesByTopologicalOrderMap.keySet());
        nodesByTopologicalOrderMap.get(0).forEach(node -> ancestorsMap.put(node, new HashSet<>()));
        for (int depth = 1; depth <= maxDepth; depth += 1) {
            for (final String node : nodesByTopologicalOrderMap.get(depth)) {
                final Set<String> nodeAncestors = new HashSet<>();
                nodeAncestors.addAll(parentsMap.get(node));
                for (final String parent : parentsMap.get(node)) {
                    nodeAncestors.addAll(ancestorsMap.get(parent));
                }
                ancestorsMap.put(node, nodeAncestors);
            }
        }
        return ancestorsMap;
    }

    private static Map<String, List<String>> getTopologicalOrderForNodeEvaluation(@Nonnull final Set<String> nodeKeysSet,
                                                                                  @Nonnull final Map<String, Integer> topologicalOrderMap,
                                                                                  @Nonnull final Map<String, Set<String>> ancestorsMap) {
        final Map<String, List<String>> topologicalOrderForNodeEvaluation = new HashMap<>();
        for (final String nodeKey : nodeKeysSet) {
            final List<String> allParentsIncludingTheNode = new ArrayList<>();
            allParentsIncludingTheNode.addAll(ancestorsMap.get(nodeKey));
            allParentsIncludingTheNode.add(nodeKey);
            /* sort by depth */
            allParentsIncludingTheNode.sort(Comparator.comparingInt(topologicalOrderMap::get));
            topologicalOrderForNodeEvaluation.put(nodeKey, allParentsIncludingTheNode);
        }
        return topologicalOrderForNodeEvaluation;
    }

    private static Map<String, List<String>> getTopologicalOrderForTagEvaluation(@Nonnull final Set<String> nodeTagsSet,
                                                                                 @Nonnull final Map<String, Integer> topologicalOrderMap,
                                                                                 @Nonnull final Map<String, Set<String>> nodesByTagMap,
                                                                                 @Nonnull final Map<String, Set<String>> ancestorsMap) {
        final Map<String, List<String>> topologicalOrderForTagEvaluation = new HashMap<>();
        for (final String tag : nodeTagsSet) {
            final Set<String> allParentsIncludingTheNodesSet = new HashSet<>();
            for (final String node : nodesByTagMap.get(tag)) {
                allParentsIncludingTheNodesSet.addAll(ancestorsMap.get(node));
                allParentsIncludingTheNodesSet.add(node);
            }
            final List<String> allParentsIncludingTheNodesList = new ArrayList<>();
            allParentsIncludingTheNodesList.addAll(allParentsIncludingTheNodesSet);
            /* sort by depth */
            allParentsIncludingTheNodesList.sort(Comparator.comparingInt(topologicalOrderMap::get));
            topologicalOrderForTagEvaluation.put(tag, allParentsIncludingTheNodesList);
        }
        return topologicalOrderForTagEvaluation;
    }

    private static List<String> getTopologicalOrderForCompleteEvaluation(@Nonnull final Set<String> nodeKeysSet,
                                                                         @Nonnull final Map<String, Integer> topologicalOrderMap) {
        return new ArrayList<>(nodeKeysSet).stream()
                .sorted(Comparator.comparingInt(topologicalOrderMap::get))
                .collect(Collectors.toList());
    }

    private static Map<String, List<String>> getTopologicalOrderForNodeMutation(@Nonnull final Set<String> nodeKeysSet,
                                                                                @Nonnull final Map<String, Set<String>> ancestorsMap,
                                                                                @Nonnull final Map<String, Set<String>> descendantsMap,
                                                                                @Nonnull final Map<String, Integer> topologicalOrderMap) {
        final Map<String, List<String>> topologicalOrderForNodeMutation = new HashMap<>();
        for (final String mutatedNodeKey : nodeKeysSet) {
            final Set<String> allInvolvedSet = new HashSet<>();
            allInvolvedSet.add(mutatedNodeKey);
            for (final String desc : descendantsMap.get(mutatedNodeKey)) {
                allInvolvedSet.add(desc);
                allInvolvedSet.addAll(ancestorsMap.get(desc));
            }
            final List<String> allInvolvedList = new ArrayList<>();
            allInvolvedList.addAll(allInvolvedSet);
            /* sort by depth */
            allInvolvedList.sort(Comparator.comparingInt(topologicalOrderMap::get));
            topologicalOrderForNodeMutation.put(mutatedNodeKey, allInvolvedList);
        }
        return topologicalOrderForNodeMutation;
    }

    /**
     * Updates the depth of a node recursively
     * */
    private static void updateDepth(@Nonnull final String nodeKey, final int recursion,
                                    @Nonnull Set<String> nodeKeysSet,
                                    @Nonnull Map<String, Set<String>> parentsMap,
                                    @Nonnull Map<String, Integer> topologicalOrderMap) {
        if (recursion > nodeKeysSet.size()) {
            throw new CyclicGraphException("The graph is not acyclic");
        }
        if (parentsMap.get(nodeKey).isEmpty()) {
            topologicalOrderMap.put(nodeKey, 0);
        } else if (topologicalOrderMap.get(nodeKey) == UNDEFINED_TOPOLOGICAL_ORDER) {
            parentsMap.get(nodeKey).forEach(parentNodeKey -> updateDepth(parentNodeKey,
                    recursion + 1, nodeKeysSet, parentsMap, topologicalOrderMap));
            final int maxParentDepth = parentsMap.get(nodeKey).stream()
                    .map(topologicalOrderMap::get)
                    .max(Integer::compareTo)
                    .get(); /* guaranteed to have a value */
            topologicalOrderMap.put(nodeKey, maxParentDepth + 1);
        }
        /* do nothing otherwise -- we already have the order for this node */
    }

    public Set<String> getNodeKeysSet() { return nodeKeysSet; }

    public Set<String> getNodeTagsSet() { return nodeTagsSet; }

    public Set<String> getInducedTagsForNode(final String nodeKey) {
        return inducedTagsMap.get(nodeKey);
    }

    public int getTopologicalOrder(@Nonnull final String nodeKey) {
        return topologicalOrderMap.get(nodeKey);
    }

    public Set<String> getChildren(@Nonnull final String nodeKey) {
        return childrenMap.get(nodeKey);
    }

    public Set<String> getParents(@Nonnull final String nodeKey) {
        return parentsMap.get(nodeKey);
    }

    public Set<String> getAncestors(@Nonnull final String nodeKey) {
        return ancestorsMap.get(nodeKey);
    }

    public Set<String> getDescendants(@Nonnull final String nodeKey) {
        return descendantsMap.get(nodeKey);
    }

    public List<String> getTopologicalOrderForNodeEvaluation(final String nodeKey) {
        return topologicalOrderForNodeEvaluation.get(nodeKey);
    }

    public List<String> getTopologicalOrderForNodeMutation(final String nodeKey) {
        return topologicalOrderForNodeMutation.get(nodeKey);
    }

    public List<String> getTopologicalOrderForTagEvaluation(final String tagKey) {
        return topologicalOrderForTagEvaluation.get(tagKey);
    }

    public List<String> getTopologicalOrderForCompleteEvaluation() {
        return topologicalOrderForCompleteEvaluation;
    }

    /**
     * This exception will be thrown if the graph has loops
     */
    public static final class CyclicGraphException extends RuntimeException {
        private static final long serialVersionUID = 5887360871425098163L;

        public CyclicGraphException(String s) {
            super(s);
        }
    }
}