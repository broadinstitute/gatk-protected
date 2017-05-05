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
    private final Map<String, Integer> depthsMap;
    private final Map<Integer, Set<String>> nodesByDepthMap;
    private final Map<String, Set<String>> nodesByTagMap;
    private final Map<String, List<String>> topologicalOrderForNodeEvaluation;
    private final Map<String, List<String>> topologicalOrderForNodeMutation;
    private final Map<String, List<String>> topologicalOrderForTagEvaluation;
    private final List<String> topologicalOrderForCompleteEvaluation;

    /**
     * An arbitrary negative number to denote the to-be-determined depth of a node
     */
    private static final int UNDEFINED_DEPTH = -1;

    /**
     * Package-private constructor from a set of {@link CacheNode}s. The graph is specified by the immediate
     * parents and descendents of each node. A {@link CyclicGraphException} is thrown of the graph has a cycle.
     *
     * @param nodeSet a set of {@link CacheNode}s
     */
    ComputableGraphStructure(@Nonnull final Set<CacheNode> nodeSet) {
        /* create the set of keys */
        nodeKeysSet = nodeSet.stream().map(CacheNode::getKey).collect(Collectors.toSet());

        /* assert that all parents exist */
        for (final CacheNode node : nodeSet) {
            Utils.validateArg(nodeKeysSet.containsAll(node.getParents()), () -> {
                final Set<String> undefinedParents = Sets.difference(new HashSet<>(node.getParents()), nodeKeysSet);
                return "Node " + ImmutableComputableGraphUtils.quote(node.getKey()) + " depends on undefined parent(s): " +
                        undefinedParents.stream().map(ImmutableComputableGraphUtils::quote).collect(Collectors.joining(", "));
            });
        }

        /* create maps for descendents, parents, and tags */
        immediateDescendentsMap = new HashMap<>();
        immediateParentsMap = new HashMap<>();
        final Map<String, Set<String>> initialTagsMap = new HashMap<>();
        nodeKeysSet.forEach(key -> {
            immediateDescendentsMap.put(key, new HashSet<>());
            immediateParentsMap.put(key, new HashSet<>());
            initialTagsMap.put(key, new HashSet<>());
        });

        /* immediate parents, descendents and tags */
        nodeSet.forEach(node -> {
            final String nodeKey = node.getKey();
            node.getParents().forEach(parent -> immediateDescendentsMap.get(parent).add(nodeKey));
            immediateParentsMap.get(nodeKey).addAll(node.getParents());
            initialTagsMap.get(nodeKey).addAll(node.getTags());
        });

        /* create maps for descendents, parents, and depthsMap */
        allDescendentsMap = new HashMap<>();
        allParentsMap = new HashMap<>();
        depthsMap = new HashMap<>();
        allTagsMap = new HashMap<>();
        nodeKeysSet.forEach(key -> {
            allDescendentsMap.put(key, new HashSet<>());
            allParentsMap.put(key, new HashSet<>());
            allTagsMap.put(key, new HashSet<>());
            depthsMap.put(key, UNDEFINED_DEPTH);
        });

        /* calculate the depth of each node; nodes with empty immediate parents have depth 0 (these include primitive nodes) */
        nodeKeysSet.forEach(nodeKey -> updateDepth(nodeKey, 0));
        final int maxDepth = Collections.max(depthsMap.values());

        /* list of nodes by their depthsMap */
        nodesByDepthMap = new HashMap<>();
        IntStream.range(0, maxDepth + 1).forEach(depth -> nodesByDepthMap.put(depth,
                nodeKeysSet.stream().filter(node -> depthsMap.get(node) == depth).collect(Collectors.toSet())));

        /* all descendents of the deepest nodes (empty set) */
        nodesByDepthMap.get(maxDepth).forEach(node -> allDescendentsMap.put(node, new HashSet<>()));
        /* get all descendents by descending the tree */
        for (int depth = maxDepth - 1; depth >= 0; depth -= 1) {
            for (final String node : nodesByDepthMap.get(depth)) {
                final Set<String> nodeAllDescendents = new HashSet<>();
                nodeAllDescendents.addAll(immediateDescendentsMap.get(node));
                for (final String child : immediateDescendentsMap.get(node)) {
                    nodeAllDescendents.addAll(allDescendentsMap.get(child));
                }
                allDescendentsMap.put(node, nodeAllDescendents);
            }
        }

        /* all parents of the primitive nodes (empty set) */
        nodesByDepthMap.get(0).forEach(node -> allParentsMap.put(node, new HashSet<>()));
        for (int depth = 1; depth <= maxDepth; depth += 1) {
            for (final String node : nodesByDepthMap.get(depth)) {
                final Set<String> nodeAllParents = new HashSet<>();
                nodeAllParents.addAll(immediateParentsMap.get(node));
                for (final String parent : immediateParentsMap.get(node)) {
                    nodeAllParents.addAll(allParentsMap.get(parent));
                }
                allParentsMap.put(node, nodeAllParents);
            }
        }

        /* build the full tags map; the parents inherit the tags of the descendents */
        nodeKeysSet.forEach(node -> allTagsMap.get(node).addAll(initialTagsMap.get(node)));
        for (int depth = maxDepth - 1; depth >= 0; depth--) {
            nodesByDepthMap.get(depth)
                    .forEach(node -> immediateDescendentsMap.get(node)
                            .forEach(desc -> allTagsMap.get(node).addAll(allTagsMap.get(desc))));
        }

        /* build a nodes-by-tag map and nodeTagsSet */
        nodeTagsSet = initialTagsMap.values().stream().flatMap(Set::stream).collect(Collectors.toSet());
        nodesByTagMap = new HashMap<>();
        nodeTagsSet.forEach(tag ->
                nodesByTagMap.put(tag, new HashSet<>()));
        nodeKeysSet.forEach(node ->
                allTagsMap.get(node).forEach(tag ->
                        nodesByTagMap.get(tag).add(node)));

        /* topological order for evaluating a single node */
        topologicalOrderForNodeEvaluation = new HashMap<>();
        for (final String node : nodeKeysSet) {
            final List<String> allParentsIncludingTheNode = new ArrayList<>();
            allParentsIncludingTheNode.addAll(allParentsMap.get(node));
            allParentsIncludingTheNode.add(node);
            /* sort by depth */
            allParentsIncludingTheNode.sort(Comparator.comparingInt(depthsMap::get));
            topologicalOrderForNodeEvaluation.put(node, allParentsIncludingTheNode);
        }

        /* topological order for evaluating all nodes associated to a tag */
        topologicalOrderForTagEvaluation = new HashMap<>();
        for (final String tag : nodeTagsSet) {
            final Set<String> allParentsIncludingTheNodesSet = new HashSet<>();
            for (final String node : nodesByTagMap.get(tag)) {
                allParentsIncludingTheNodesSet.addAll(allParentsMap.get(node));
                allParentsIncludingTheNodesSet.add(node);
            }
            final List<String> allParentsIncludingTheNodesList = new ArrayList<>();
            allParentsIncludingTheNodesList.addAll(allParentsIncludingTheNodesSet);
            /* sort by depth */
            allParentsIncludingTheNodesList.sort(Comparator.comparingInt(depthsMap::get));
            topologicalOrderForTagEvaluation.put(tag, allParentsIncludingTheNodesList);
        }

        /* topological order for evaluating all nodes */
        topologicalOrderForCompleteEvaluation = new ArrayList<>();
        topologicalOrderForCompleteEvaluation.addAll(nodeKeysSet);
        topologicalOrderForCompleteEvaluation.sort(Comparator.comparingInt(depthsMap::get));

        /* topological order for updating the descendents of a mutated node */
        topologicalOrderForNodeMutation = new HashMap<>();
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
            allInvolvedList.sort(Comparator.comparingInt(depthsMap::get));
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
            throw new CyclicGraphException("The graph has cycles");
        }
        if (immediateParentsMap.get(nodeKey).isEmpty()) {
            depthsMap.put(nodeKey, 0);
        } else if (depthsMap.get(nodeKey) == UNDEFINED_DEPTH) {
            immediateParentsMap.get(nodeKey).forEach(parentNodeKey -> updateDepth(parentNodeKey, recursion + 1));
            final int maxParentDepth = immediateParentsMap.get(nodeKey).stream()
                    .map(depthsMap::get)
                    .max(Integer::compareTo)
                    .get(); /* guaranteed to have a value */
            depthsMap.put(nodeKey, maxParentDepth + 1);
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
                    "\tdepth: " + depthsMap.get(nodeKey) + "\n" +
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
                            map(nodeKey -> nodeKey + "(" + depthsMap.get(nodeKey) + ")").
                            collect(Collectors.joining(", ", "[", "]\n"));
        }

        status += "\n";
        for (final String node : nodeKeysSet) {
            status += "topological order evaluating node: " + node + ", nodes: " +
                    topologicalOrderForNodeEvaluation.get(node).stream().
                            map(nodeKey -> nodeKey + "(" + depthsMap.get(nodeKey) + ")").
                            collect(Collectors.joining(", ", "[", "]\n"));
        }

        status += "\n";
        for (final String node : nodesByDepthMap.get(0)) {
            status += "topological order for node mutation: " + node + ", nodes: " +
                    topologicalOrderForNodeMutation.get(node).stream().
                            map(nodeKey -> nodeKey + "(" + depthsMap.get(nodeKey) + ")").
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