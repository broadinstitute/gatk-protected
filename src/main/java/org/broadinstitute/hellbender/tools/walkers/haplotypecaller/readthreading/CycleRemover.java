package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import org.broadinstitute.hellbender.utils.Utils;
import org.jgrapht.DirectedGraph;

import java.util.*;

/**
 * Resolves cycles detected between sources and sinks of a graph.
 */
public class CycleRemover<V, E> {
    private final DirectedGraph<V, E> graph;

    public CycleRemover(final DirectedGraph<V, E> graph) {
        this.graph = Utils.nonNull(graph);
    }

    /**
     * Removes edges that produce cycles in paths from a source to a sink vertex.
     * @param source
     * @param sink
     * @return never {@code true} iff this operation has changed the graph.
     */
    public boolean removeCycles(final V source, final V sink) {
        Utils.nonNull(source);
        Utils.nonNull(sink);
        if (!graph.containsVertex(source)) {
            throw new IllegalArgumentException("source must be part of the graph");
        } else if (!graph.containsVertex(sink)) {
            throw new IllegalArgumentException("sink must be part of the graph");
        }
        return removeCycles(Collections.singleton(source), Collections.singleton(sink));
    }

    /**
     * Removes edges that produces cycles and also dead vertices that do not lead to any sink vertex.
     *
     * @param sources considered source vertices.
     * @param sinks considered sink vertices.
     * @return never {@code null}.
     */
    public boolean removeCycles(final Collection<V> sources, final Collection<V> sinks) {
        final Set<E> edgesToRemove = new HashSet<>(graph.edgeSet().size());
        final Set<V> vertexToRemove = new HashSet<>(graph.vertexSet().size());

        boolean foundSomePath = false;
        for (final V source : sources) {
            final Set<V> parentVertices = new HashSet<>(graph.vertexSet().size());
            foundSomePath = findGuiltyVerticesAndEdgesToRemoveCycles(source, sinks, edgesToRemove, vertexToRemove, parentVertices) || foundSomePath;
        }

        if (!foundSomePath) {
            throw new IllegalStateException("could not find any path from the source vertex to the sink vertex after removing cycles: "
                    + Arrays.toString(sources.toArray()) + " => " + Arrays.toString(sinks.toArray()));
        }

        if (edgesToRemove.isEmpty() && vertexToRemove.isEmpty()) {
            return false;
        } else {
            graph.removeAllEdges(edgesToRemove);
            graph.removeAllVertices(vertexToRemove);
            return true;
        }
    }

    /**
     * Recursive call that looks for edges and vertices that need to be removed to get rid of cycles.
     *
     * @param currentVertex current search vertex.
     * @param sinks considered sink vertices.
     * @param edgesToRemove collection  of edges that need to be removed in order to get rid of cycles.
     * @param verticesToRemove collection of vertices that can be removed.
     * @param parentVertices collection of vertices that preceded the {@code currentVertex}; i.e. the it can be
     *                       reached from those vertices using edges existing in {@code graph}.
     *
     * @return {@code true} to indicate that the some sink vertex is reachable by {@code currentVertex},
     *  {@code false} otherwise.
     */
    private boolean findGuiltyVerticesAndEdgesToRemoveCycles(final V currentVertex,
                                                             final Collection<V> sinks,
                                                             final Set<E> edgesToRemove,
                                                             final Set<V> verticesToRemove,
                                                             final Set<V> parentVertices) {
        if (sinks.contains(currentVertex)) {
            return true;
        }

        final Set<E> outgoingEdges = graph.outgoingEdgesOf(currentVertex);
        parentVertices.add(currentVertex);

        boolean reachesSink = false;
        for (final E edge : outgoingEdges) {
            final V child = graph.getEdgeTarget(edge);
            if (parentVertices.contains(child)) {
                edgesToRemove.add(edge);
            } else {
                final boolean childReachSink = findGuiltyVerticesAndEdgesToRemoveCycles(child, sinks, edgesToRemove, verticesToRemove, parentVertices);
                reachesSink = reachesSink || childReachSink;
            }
        }
        parentVertices.remove(currentVertex);
        if (!reachesSink) {
            verticesToRemove.add(currentVertex);
        }
        return reachesSink;
    }

}
