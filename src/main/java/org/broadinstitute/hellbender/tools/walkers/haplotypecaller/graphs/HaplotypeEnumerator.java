package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import it.unimi.dsi.fastutil.Hash;
import it.unimi.dsi.fastutil.objects.Object2DoubleMap;
import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2LongLinkedOpenCustomHashMap;
import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Haplotype enumerator component.
 */
public class HaplotypeEnumerator {

    private final Logger logger;

    private static final Hash.Strategy<SeqVertex> VERTEX_HASH_STRATEGY = new Hash.Strategy<SeqVertex>() {
        @Override
        public int hashCode(SeqVertex seqVertex) {
            return Integer.hashCode(seqVertex.getId());
        }

        @Override
        public boolean equals(SeqVertex v1, SeqVertex v2) {
            return v1 == v2;
        }
    };

    public static final int MAXIMUM_NO_DISCOVERY_RUN_LENGTH_SHORT_NAME_MINIMUM_RECOMMENDED = 50;
    public static final int MAXIMUM_NO_DISCOVERY_RUN_LENGTH_SHORT_NAME_DEFAULT = MAXIMUM_NO_DISCOVERY_RUN_LENGTH_SHORT_NAME_MINIMUM_RECOMMENDED * 2;
    public static final int MAXIMUM_HAPLOTYPES_PER_ASSEMBLY_GRAPH_DEFAULT = 128;
    public static final String MAXIMUM_HAPLOTYPES_PER_ASSEMBLY_GRAPH_SHORT_NAME = "maximumHaplotypesPerAssemblyGraph";
    public static final String MAXIMUM_HAPLOTYPES_PER_ASSEMBLY_GRAPH_FULL_NAME = MAXIMUM_HAPLOTYPES_PER_ASSEMBLY_GRAPH_SHORT_NAME;
    public static final String MAXIMUM_NO_DISCOVERY_RUN_LENGTH_SHORT_NAME = "maximumHaplotypeSamplingWithNoDiscoveryInComplexAssemblyGraph";
    public static final String MAXIMUM_NO_DISCOVERY_RUN_LENGTH_FULL_NAME = MAXIMUM_NO_DISCOVERY_RUN_LENGTH_SHORT_NAME;

    @Argument(doc = "maximum number of haplotypes considered per assembly graph",
            shortName = MAXIMUM_HAPLOTYPES_PER_ASSEMBLY_GRAPH_SHORT_NAME,
            fullName = MAXIMUM_HAPLOTYPES_PER_ASSEMBLY_GRAPH_FULL_NAME, optional = true)
    private int maximumHaplotypeCount = MAXIMUM_HAPLOTYPES_PER_ASSEMBLY_GRAPH_DEFAULT;

    @Argument(doc = "maximum number of haplotypes samples without newly added haplotypes before stopping looking for more",
            shortName = MAXIMUM_NO_DISCOVERY_RUN_LENGTH_SHORT_NAME,
            fullName = MAXIMUM_NO_DISCOVERY_RUN_LENGTH_FULL_NAME)
    private int maximumNoDiscoveryRunLength = MAXIMUM_NO_DISCOVERY_RUN_LENGTH_SHORT_NAME_DEFAULT;

    public HaplotypeEnumerator(final int maximumHaplotypeCount, final int maximumNoDiscoveryRunLength) {
        this();
        this.maximumHaplotypeCount = maximumHaplotypeCount;
        this.maximumNoDiscoveryRunLength = maximumNoDiscoveryRunLength;
        checkUserArguments();
    }

    public HaplotypeEnumerator() {
        logger = LogManager.getLogger(this);
    }

    public List<Haplotype> list(final SeqGraph graph) {
        checkUserArguments();
        final long actualHaplotypeCount = calculateNumberOfHaplotypes(graph);
        if (actualHaplotypeCount > maximumHaplotypeCount) {
            final List<Haplotype> result = sampleHaplotypes(graph, maximumHaplotypeCount, maximumNoDiscoveryRunLength);
            logger.debug("number of haplotypes in assembly graph (%d) is larger than the maximum allowed (%d); stochastic haplotype sampling produced (%d) haplotypes with a maximum no-discovery run of (%d)",
                    actualHaplotypeCount, maximumHaplotypeCount, result.size(), maximumNoDiscoveryRunLength);
            return result;
        } else {
            return listAllHaplotypes(graph);
        }
    }

    private void checkUserArguments() {
        if (maximumHaplotypeCount < 2) {
            throw new UserException.BadArgumentValue(MAXIMUM_HAPLOTYPES_PER_ASSEMBLY_GRAPH_SHORT_NAME, "" + maximumHaplotypeCount, "it must equal or greater than 2");
        } else if (maximumNoDiscoveryRunLength < 1) {
            throw new UserException.BadArgumentValue(MAXIMUM_NO_DISCOVERY_RUN_LENGTH_SHORT_NAME, "" + maximumHaplotypeCount, "it must equal or greater than 1");
        } else if (maximumNoDiscoveryRunLength < 100) {
            logger.warn(String.format("The value give to %s (=%d) is rather small and may cause to miss plausible haplotypes in complex regions where the number of haplotypes in the assembly graph goes beyond %s (=%d) MAXIMUM_HAPLOTYPES_PER_ASSEMBLY_GRAPH_SHORT_NAME. A value of at least %d is recommendable",
                    MAXIMUM_NO_DISCOVERY_RUN_LENGTH_SHORT_NAME, maximumNoDiscoveryRunLength,
                    MAXIMUM_HAPLOTYPES_PER_ASSEMBLY_GRAPH_SHORT_NAME, maximumHaplotypeCount,
                    MAXIMUM_NO_DISCOVERY_RUN_LENGTH_SHORT_NAME_MINIMUM_RECOMMENDED));
        }
    }

    private static List<Haplotype> listAllHaplotypes(final SeqGraph graph) {
        final Map<SeqVertex, List<HaplotypeSuffix>> subProblems = new HashMap<>(graph.vertexSet().size());
        listAllHaplotypes(graph, graph.getReferenceSourceVertex(), subProblems);
        return subProblems.get(graph.getReferenceSourceVertex()).stream()
                .map(HaplotypeSuffix::toHaplotype)
                .collect(Collectors.toList());
    }

    private static List<HaplotypeSuffix> listAllHaplotypes(final SeqGraph graph, final SeqVertex vertex, final Map<SeqVertex, List<HaplotypeSuffix>> subProblems) {
        if (subProblems.containsKey(vertex)) {
            return subProblems.get(vertex);
        } else if (graph.isSink(vertex)) {
            if (vertex != graph.getReferenceSinkVertex()) {
                subProblems.put(vertex, Collections.emptyList());
                return Collections.emptyList();
            } else {
                subProblems.put(vertex, Collections.singletonList(new HaplotypeSuffix(vertex.getSequence())));
                return subProblems.get(vertex);
            }
        } else {
            final List<HaplotypeSuffix> result = new ArrayList<>();
            final byte[] head = vertex.getSequence();
            final double totalCount = outgoingEdgesMultiplicitySum(graph, vertex);
            for (final BaseEdge edge : graph.outgoingEdgesOf(vertex)) {
                final boolean isReference = edge.isRef();
                final double score = Math.log10(Math.max(edge.getMultiplicity(), 0.5) / totalCount);
                result.addAll(listAllHaplotypes(graph, graph.getEdgeTarget(edge), subProblems).stream()
                    .map(tail -> new HaplotypeSuffix(head, isReference, score, tail))
                    .collect(Collectors.toList()));
            }
            subProblems.put(vertex, result);
            return result;
        }
    }

    private static List<Haplotype> sampleHaplotypes(final SeqGraph graph, final int maximumHaplotypeCount, final int maximumNoDiscoveryRunLength) {
        final Set<Haplotype> result = new LinkedHashSet<>(maximumHaplotypeCount);
        final Object2DoubleMap<SeqVertex> log10TotalOutgoingMultiplicity = log10TotalOutgoingMultiplicityMap(graph);
        final Haplotype refHaplotype = composeReferenceHaplotype(graph, log10TotalOutgoingMultiplicity);
        result.add(refHaplotype);
        final int bufferSize = refHaplotype.length() << 1;
        final byte[] buffer = new byte[bufferSize];
        int noDiscoveryRun = 0;
        final Random rdn = Utils.getRandomGenerator();
        while (noDiscoveryRun < maximumNoDiscoveryRunLength && result.size() < maximumHaplotypeCount) {
            final Haplotype haplotype = sampleOneHaplotype(graph, rdn, buffer, log10TotalOutgoingMultiplicity);
            if (result.add(haplotype)) {
                noDiscoveryRun = 0;
            } else {
                noDiscoveryRun++;
            }
        }
        return new ArrayList<>(result);
    }

    private static Haplotype composeReferenceHaplotype(final SeqGraph graph,
                                                       final Object2DoubleMap<SeqVertex> log10TotalOutgoingMultiplicity) {
        SeqVertex next = graph.getReferenceSourceVertex();
        final Stack<SeqVertex> stack = new Stack<>();
        while (next != null) {
            stack.push(next);
            next = graph.getNextReferenceVertex(next);
        }

        SeqVertex target = stack.pop();
        HaplotypeSuffix suffix = new HaplotypeSuffix(target.getSequence());
        while (!stack.isEmpty()) {
            final SeqVertex source = stack.pop();
            final BaseEdge edge = graph.getEdge(source, target);
            suffix = new HaplotypeSuffix(source.getSequence(), edge.isRef(),
                    Math.log10(Math.max(0.5, edge.getMultiplicity())) -  log10TotalOutgoingMultiplicity.getDouble(source), suffix);
            target = source;
        }
        return suffix.toHaplotype();
    }

    private static Object2DoubleMap<SeqVertex> log10TotalOutgoingMultiplicityMap(final SeqGraph graph) {
        final Object2DoubleMap<SeqVertex> result = new Object2DoubleOpenHashMap<>(graph.vertexSet().size());
        for (final SeqVertex vertex : graph.vertexSet()) {
            result.put(vertex, Math.log10(outgoingEdgesMultiplicitySum(graph, vertex)));
        }
        return result;
    }

    private static double outgoingEdgesMultiplicitySum(final SeqGraph graph, final SeqVertex source) {
        double result = 0;
        for (final BaseEdge edge : graph.outgoingEdgesOf(source)) {
            result += Math.max(0.5, edge.getMultiplicity());
        }
        return result;
    }

    private static Haplotype sampleOneHaplotype(final SeqGraph graph, final Random rdn, byte[] basesBuffer,
                                                final Object2DoubleMap<SeqVertex> log10TotalOutgoingMultiplicity) {
        SeqVertex next = graph.getReferenceSourceVertex();
        double score = 0;
        int nextBase = 0;
        boolean isReference = true;
        while (next != null) {
            final byte[] nextVertexBases = next.getSequence();
            final int nextBaseUpdate = nextBase + nextVertexBases.length;
            if (nextBaseUpdate >= basesBuffer.length) { // just in case we run out of space; very rare occurrence.
                basesBuffer = Arrays.copyOf(basesBuffer, nextBaseUpdate << 1);
            }
            System.arraycopy(nextVertexBases, 0, basesBuffer, nextBase, nextVertexBases.length);
            nextBase = nextBaseUpdate;
            final Set<BaseEdge> outgoingEdges = graph.outgoingEdgesOf(next);
            final double log10TotalSum = log10TotalOutgoingMultiplicity.get(next);
            final double totalSum = Math.pow(10, log10TotalSum);
            double random = rdn.nextDouble() * totalSum;
            final Iterator<BaseEdge> outgoingEdgesIterator = outgoingEdges.iterator();
            BaseEdge nextChildEdge = null;
            while (random >= 0 && outgoingEdgesIterator.hasNext()) {
                nextChildEdge = outgoingEdgesIterator.next();
                random -= nextChildEdge.getMultiplicity();
            }
            if (nextChildEdge == null) {
                next = null;
            } else {
                next = graph.getEdgeTarget(nextChildEdge);
                isReference &= nextChildEdge.isRef();
                score += Math.log10(Math.max(0.5, nextChildEdge.getMultiplicity())) - log10TotalSum;
            }
        }
        final Haplotype result = new Haplotype(Arrays.copyOf(basesBuffer, nextBase), isReference);
        result.setScore(score);
        return result;
    }

    private static long calculateNumberOfHaplotypes(final SeqGraph graph) {
        final Object2LongMap<SeqVertex> subProblems = new Object2LongLinkedOpenCustomHashMap<SeqVertex>(graph.vertexSet().size(), VERTEX_HASH_STRATEGY);
        return calculateNumberOfHaplotypes(graph, graph.getReferenceSourceVertex(), subProblems);
    }

    private static long calculateNumberOfHaplotypes(final SeqGraph graph, final SeqVertex vertex, final Object2LongMap<SeqVertex> subProblems) {
        if (subProblems.containsKey(vertex)) {
            return subProblems.getLong(vertex);
        } else if (graph.isSink(vertex)) {
            if (graph.isReferenceNode(vertex)) {
                subProblems.put(vertex, 1);
                return 1;
            } else {
                subProblems.put(vertex, 0);
                return 0;
            }
        } else {
            long result = 0;
            for (final BaseEdge edge : graph.outgoingEdgesOf(vertex)) {
                final long childResult = calculateNumberOfHaplotypes(graph, graph.getEdgeTarget(edge), subProblems);
                if (Long.MAX_VALUE - result <= childResult) {
                    subProblems.put(vertex, Long.MAX_VALUE);
                    return Long.MAX_VALUE;
                } else {
                    result += childResult;
                }
            }
            subProblems.put(vertex, result);
            return result;
        }
    }

    private static class HaplotypeSuffix {
        private final double score;
        private final byte[] head;
        private final HaplotypeSuffix tail;
        private final boolean isReference;
        private final int length;

        private HaplotypeSuffix(final byte[] head) {
            this.length = head.length;
            this.score = 0;
            this.head = head;
            this.tail = null;
            this.isReference = true;
        }

        private HaplotypeSuffix(final byte[] head, final boolean isReference, final double score, final HaplotypeSuffix tail) {
            this.length = head.length + tail.length;
            this.head = head;
            this.score = score + tail.score;
            this.tail = tail;
            this.isReference = isReference && tail.isReference;
        }

        public Haplotype toHaplotype() {
            final byte[] bases = new byte[length];
            HaplotypeSuffix next = this;
            int nextBase = 0;
            while (next != null) {
                System.arraycopy(next.head, 0, bases, nextBase, next.head.length);
                nextBase += next.head.length;
                next = next.tail;
            }
            final Haplotype result = new Haplotype(bases, isReference);
            result.setScore(score);
            return result;
        }
    }
}
