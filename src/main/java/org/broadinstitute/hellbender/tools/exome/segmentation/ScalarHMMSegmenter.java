package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.apache.commons.collections4.ListUtils;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * HMM segmenter for scalar hidden data, such as allele fraction and copy ratio, but not the joint
 * allele fraction / copy ratio segmenter
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public abstract class ScalarHMMSegmenter<DATA> extends ClusteringGenomicHMMSegmenter<DATA, Double> {

    private static final double DEFAULT_MEMORY_LENGTH = 5e7;

    public ScalarHMMSegmenter(final List<SimpleInterval> positions, final List<DATA> data, final List<Double> initialNonConstantHiddenStates) {
        super(positions, data, initialNonConstantHiddenStates, DEFAULT_MEMORY_LENGTH);
    }

    @Override
    protected boolean hiddenStateValuesHaveConverged(final List<Double> oldHiddenStateValues) {
        return oldHiddenStateValues.size() == numStates() && GATKProtectedMathUtils.maxDifference(oldHiddenStateValues, getStates()) < CONVERGENCE_THRESHOLD;
    }

    @Override
    protected void relearnHiddenStateValues(final ExpectationStep eStep) {
        final ClusteringGenomicHMM<DATA, Double> model = makeModel();
        IntStream.range(0, numStates()).forEach(state -> {
            final Function<Double, Double> objective = f -> IntStream.range(0, data.size())
                    .filter(n -> eStep.pStateAtPosition(state, n) > NEGLIGIBLE_POSTERIOR_FOR_M_STEP)
                    .mapToDouble(n -> eStep.pStateAtPosition(state, n) * model.logEmissionProbability(data.get(n), f))
                    .sum();
            setState(state, OptimizationUtils.singleNewtonArgmaxUpdate(objective, minHiddenStateValue(),
                    maxHiddenStateValue(), getState(state)));
        });
    }

    protected abstract double minHiddenStateValue();
    protected abstract double maxHiddenStateValue();
}
