package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import avro.shaded.com.google.common.collect.Sets;
import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class pre-computes a number of useful auxiliary properties for the DAG specified by
 * a set of {@link CacheNode}s. These include:
 *
 * - topological order for evaluating a computable node,
 * - topological order for mutating a primitive/externally-computed node, and
 * - topological order for evaluating all nodes associated to a tag (see {@link ImmutableComputableGraph}).
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class ComputableGraphStructure implements Serializable {

    private static final long serialVersionUID = -3124293279477371159L;

    private final Set<String> nodeKeysSet;
    private final Set<String> nodeTagsSet;
    private final Map<String, Set<String>> allTagsMap;
    private final Map<String, Set<String>> immediateDescendentsMap;
    private final Map<String, Set<String>> immediateParentsMap;
    private final Map<String, Set<String>> allDescendentsMap;
    private final Map<String, Set<String>> allParentsMap;
    private final Map<String, Integer> topologicalOrderMap;
    private final Map<Integer, Set<String>> nodesByTopologicalOrderMap;
    private final Map<String, Set<String>> nodesByTagMap;
    private final Map<String, List<String>> topologicalOrderForNodeEvaluation;
    private final Map<String, List<String>> topologicalOrderForNodeMutation;
    private final Map<String, List<String>> topologicalOrderForTagEvaluation;
    private final List<String> topologicalOrderForCompleteEvaluation;

    /**
     * An arbitrary negative number to denote the to-be-determined depth of a node
     */
    private static final int UNDEFINED_ORDER = -1;

    /**
     * Package-private constructor from a set of {@link CacheNode}s. The graph is specified by the immediate
     * parents and descendents of each node. A {@link CyclicGraphException} is thrown of the graph has a cycle.
     *
     * @param nodeSet a set of {@link CacheNode}s
     */
    ComputableGraphStructure(@Nonnull final Set<CacheNode> nodeSet) {
        /* create the set of keys */
        Utils.nonNull(nodeSet, "The given set of nodes must be non-null");
        nodeKeysSet = nodeSet.stream().map(CacheNode::getKey).collect(Collectors.toSet());
        nodeTagsSet = nodeSet.stream().map(CacheNode::getTags).flatMap(Collection::stream).collect(Collectors.toSet());

        /* assert that all parent keys exist */
        for (final CacheNode node : nodeSet) {
            Utils.validateArg(nodeKeysSet.containsAll(node.getParents()), () -> {
                final Set<String> undefinedParents = Sets.difference(new HashSet<>(node.getParents()), nodeKeysSet);
                return "Node " + ImmutableComputableGraphUtils.quote(node.getKey()) + " depends on undefined parent(s): " +
                        undefinedParents.stream().map(ImmutableComputableGraphUtils::quote).collect(Collectors.joining(", "));
            });
        }

        /* initialize the containers */
        immediateDescendentsMap = new HashMap<>();
        immediateParentsMap = new HashMap<>();
        allDescendentsMap = new HashMap<>();
        allParentsMap = new HashMap<>();
        topologicalOrderMap = new HashMap<>();
        allTagsMap = new HashMap<>();
        nodesByTagMap = new HashMap<>();
        topologicalOrderForCompleteEvaluation = new ArrayList<>();
        topologicalOrderForTagEvaluation = new HashMap<>();
        topologicalOrderForNodeEvaluation = new HashMap<>();
        topologicalOrderForNodeMutation = new HashMap<>();

        /* immediate parents, descendents and tags */
        initializeImmediateDescendentsAndParents(nodeSet);

        /* calculate the depth of each node; nodes with empty immediate parents have depth 0 (these include primitive nodes) */
        initializeTopologicalOrder();
        nodesByTopologicalOrderMap = getNodesByTopologicalOrderMap(topologicalOrderMap, nodeKeysSet);

        initializeAllDescendantsAndParents();

        /* build the full tags map; the parents inherit the tags of the descendents */
        processTags(nodeSet);

        /* topological order for evaluating a single node */
        initializeTopologicalOrderForNodeEvaluation();

        /* topological order for evaluating all nodes associated to a tag */
        initializeTopologicalOrderForTagEvaluation();

        /* topological order for evaluating all nodes */
        initializeTopologicalOrderForCompleteEvaluation();

        /* topological order for updating the descendents of a mutated node */
        initializeTopologicalOrderForNodeMutation();
    }

    private void initializeTopologicalOrder() {
        nodeKeysSet.forEach(key -> topologicalOrderMap.put(key, UNDEFINED_ORDER));
        nodeKeysSet.forEach(nodeKey -> updateDepth(nodeKey, 0));
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

    private void initializeImmediateDescendentsAndParents(@Nonnull Set<CacheNode> nodeSet) {
        nodeKeysSet.forEach(key -> {
            immediateDescendentsMap.put(key, new HashSet<>());
            immediateParentsMap.put(key, new HashSet<>());
        });
        nodeSet.forEach(node -> {
            final String nodeKey = node.getKey();
            node.getParents().forEach(parent -> immediateDescendentsMap.get(parent).add(nodeKey));
            immediateParentsMap.get(nodeKey).addAll(node.getParents());
        });
    }

    private void processTags(@Nonnull Set<CacheNode> nodeSet) {
        final Map<String, Set<String>> initialTagsMap = nodeSet.stream()
                .collect(Collectors.toMap(CacheNode::getKey, node -> new HashSet<>(node.getTags())));
        nodeKeysSet.forEach(key -> allTagsMap.put(key, new HashSet<>()));
        nodeKeysSet.forEach(node -> allTagsMap.get(node).addAll(initialTagsMap.get(node)));
        final int maxDepth = Collections.max(topologicalOrderMap.values());
        for (int depth = maxDepth - 1; depth >= 0; depth--) {
            nodesByTopologicalOrderMap.get(depth)
                    .forEach(node -> immediateDescendentsMap.get(node)
                            .forEach(desc -> allTagsMap.get(node).addAll(allTagsMap.get(desc))));
        }
        nodeTagsSet.forEach(tag -> nodesByTagMap.put(tag, new HashSet<>()));
        nodeKeysSet.forEach(node ->
                allTagsMap.get(node).forEach(tag ->
                        nodesByTagMap.get(tag).add(node)));
    }

    private void initializeAllDescendantsAndParents() {
        nodeKeysSet.forEach(key -> {
            allDescendentsMap.put(key, new HashSet<>());
            allParentsMap.put(key, new HashSet<>());
        });
        final int maxDepth = Collections.max(topologicalOrderMap.values());
        nodesByTopologicalOrderMap.get(maxDepth).forEach(node -> allDescendentsMap.put(node, new HashSet<>()));
        /* get all descendents by descending the tree */
        for (int depth = maxDepth - 1; depth >= 0; depth -= 1) {
            for (final String node : nodesByTopologicalOrderMap.get(depth)) {
                final Set<String> nodeAllDescendents = new HashSet<>();
                nodeAllDescendents.addAll(immediateDescendentsMap.get(node));
                for (final String child : immediateDescendentsMap.get(node)) {
                    nodeAllDescendents.addAll(allDescendentsMap.get(child));
                }
                allDescendentsMap.put(node, nodeAllDescendents);
            }
        }
        nodesByTopologicalOrderMap.get(0).forEach(node -> allParentsMap.put(node, new HashSet<>()));
        for (int depth = 1; depth <= maxDepth; depth += 1) {
            for (final String node : nodesByTopologicalOrderMap.get(depth)) {
                final Set<String> nodeAllParents = new HashSet<>();
                nodeAllParents.addAll(immediateParentsMap.get(node));
                for (final String parent : immediateParentsMap.get(node)) {
                    nodeAllParents.addAll(allParentsMap.get(parent));
                }
                allParentsMap.put(node, nodeAllParents);
            }
        }
    }

    private void initializeTopologicalOrderForNodeEvaluation() {
        for (final String node : nodeKeysSet) {
            final List<String> allParentsIncludingTheNode = new ArrayList<>();
            allParentsIncludingTheNode.addAll(allParentsMap.get(node));
            allParentsIncludingTheNode.add(node);
            /* sort by depth */
            allParentsIncludingTheNode.sort(Comparator.comparingInt(topologicalOrderMap::get));
            topologicalOrderForNodeEvaluation.put(node, allParentsIncludingTheNode);
        }
    }

    private void initializeTopologicalOrderForTagEvaluation() {
        for (final String tag : nodeTagsSet) {
            final Set<String> allParentsIncludingTheNodesSet = new HashSet<>();
            for (final String node : nodesByTagMap.get(tag)) {
                allParentsIncludingTheNodesSet.addAll(allParentsMap.get(node));
                allParentsIncludingTheNodesSet.add(node);
            }
            final List<String> allParentsIncludingTheNodesList = new ArrayList<>();
            allParentsIncludingTheNodesList.addAll(allParentsIncludingTheNodesSet);
            /* sort by depth */
            allParentsIncludingTheNodesList.sort(Comparator.comparingInt(topologicalOrderMap::get));
            topologicalOrderForTagEvaluation.put(tag, allParentsIncludingTheNodesList);
        }
    }

    private void initializeTopologicalOrderForCompleteEvaluation() {
        topologicalOrderForCompleteEvaluation.addAll(nodeKeysSet);
        topologicalOrderForCompleteEvaluation.sort(Comparator.comparingInt(topologicalOrderMap::get));
    }

    private void initializeTopologicalOrderForNodeMutation() {
        for (final String mutNode : nodeKeysSet) {
            final Set<String> allInvolvedSet = new HashSet<>();
            allInvolvedSet.add(mutNode);
            for (final String desc : allDescendentsMap.get(mutNode)) {
                allInvolvedSet.add(desc);
                allInvolvedSet.addAll(allParentsMap.get(desc));
            }
            final List<String> allInvolvedList = new ArrayList<>();
            allInvolvedList.addAll(allInvolvedSet);
            /* sort by depth */
            allInvolvedList.sort(Comparator.comparingInt(topologicalOrderMap::get));
            topologicalOrderForNodeMutation.put(mutNode, allInvolvedList);
        }
    }

    /**
     * Updates the depth of a node recursively
     *
     * @param nodeKey the key of the node to update
     */
    private void updateDepth(final String nodeKey, final int recursion) {
        if (recursion > nodeKeysSet.size()) {
            throw new CyclicGraphException("The graph is not acyclic");
        }
        if (immediateParentsMap.get(nodeKey).isEmpty()) {
            topologicalOrderMap.put(nodeKey, 0);
        } else if (topologicalOrderMap.get(nodeKey) == UNDEFINED_ORDER) {
            immediateParentsMap.get(nodeKey).forEach(parentNodeKey -> updateDepth(parentNodeKey, recursion + 1));
            final int maxParentDepth = immediateParentsMap.get(nodeKey).stream()
                    .map(topologicalOrderMap::get)
                    .max(Integer::compareTo)
                    .get(); /* guaranteed to have a value */
            topologicalOrderMap.put(nodeKey, maxParentDepth + 1);
        }
    }

    public Set<String> getNodeKeysSet() { return nodeKeysSet; }

    public Set<String> getNodeTagsSet() { return nodeTagsSet; }

    public Set<String> getAllTagsForNode(final String nodeKey) {
        return allTagsMap.get(nodeKey);
    }

    public Set<String> getAllDescendents(@Nonnull final String nodeKey) {
        return allDescendentsMap.get(nodeKey);
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

    @Override
    public String toString() {
        String status = "";
        for (final String nodeKey : nodeKeysSet) {
            status += "node: " + nodeKey + "\n" +
                    "\tdepth: " + topologicalOrderMap.get(nodeKey) + "\n" +
                    "\timmediate parents: " +
                    immediateParentsMap.get(nodeKey).stream().collect(Collectors.joining(", ", "[", "]\n")) +
                    "\timmediate descendents: " +
                    immediateDescendentsMap.get(nodeKey).stream().collect(Collectors.joining(", ", "[", "]\n")) +
                    "\tall parents: " +
                    allParentsMap.get(nodeKey).stream().collect(Collectors.joining(", ", "[", "]\n")) +
                    "\tall descendents: " +
                    allDescendentsMap.get(nodeKey).stream().collect(Collectors.joining(", ", "[", "]\n")) +
                    "\tall tags: " +
                    allTagsMap.get(nodeKey).stream().collect(Collectors.joining(", ", "[", "]\n"));
        }

        status += "\n";
        for (final String tag : nodeTagsSet) {
            status += "tag: " + tag + ", nodes: " +
                    nodesByTagMap.get(tag).stream().collect(Collectors.joining(", ", "[", "]\n"));
        }

        status += "\n";
        for (final String tag : nodeTagsSet) {
            status += "topological order for evaluating tag: " + tag + ", nodes:" +
                    topologicalOrderForTagEvaluation.get(tag).stream().
                            map(nodeKey -> nodeKey + "(" + topologicalOrderMap.get(nodeKey) + ")").
                            collect(Collectors.joining(", ", "[", "]\n"));
        }

        status += "\n";
        for (final String node : nodeKeysSet) {
            status += "topological order evaluating node: " + node + ", nodes: " +
                    topologicalOrderForNodeEvaluation.get(node).stream().
                            map(nodeKey -> nodeKey + "(" + topologicalOrderMap.get(nodeKey) + ")").
                            collect(Collectors.joining(", ", "[", "]\n"));
        }

        status += "\n";
        for (final String node : nodesByTopologicalOrderMap.get(0)) {
            status += "topological order for node mutation: " + node + ", nodes: " +
                    topologicalOrderForNodeMutation.get(node).stream().
                            map(nodeKey -> nodeKey + "(" + topologicalOrderMap.get(nodeKey) + ")").
                            collect(Collectors.joining(", ", "[", "]\n"));
        }
        return status;
    }

    /**
     * This exception will be thrown if the graph has loops
     */
    public final class CyclicGraphException extends RuntimeException {
        private static final long serialVersionUID = 5887360871425098163L;

        public CyclicGraphException(String s) {
            super(s);
        }
    }
}