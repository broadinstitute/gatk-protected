package org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class VariantPloidyStatePrior {
    private final int numVariantPloidyStates;
    private final Map<VariantPloidyState, Double> logProbabilityMassFunctionMap;

    public VariantPloidyStatePrior(final Map<VariantPloidyState, Double> unnormalizedLogProbabilityMassFunctionMap) {
        Utils.nonNull(unnormalizedLogProbabilityMassFunctionMap);
        Utils.validateArg(!unnormalizedLogProbabilityMassFunctionMap.isEmpty(), "Number of variant ploidy states must be positive.");
        numVariantPloidyStates = unnormalizedLogProbabilityMassFunctionMap.size();
        this.logProbabilityMassFunctionMap = normalize(new LinkedHashMap<>(unnormalizedLogProbabilityMassFunctionMap));
    }

    public int numVariantPloidyStates() {
        return numVariantPloidyStates;
    }

    public double logProbability(final VariantPloidyState variantPloidyState) {
        Utils.nonNull(variantPloidyState);
        if (logProbabilityMassFunctionMap.containsKey(variantPloidyState)) {
            return logProbabilityMassFunctionMap.get(variantPloidyState);
        } else {
            throw new IllegalArgumentException(String.format("Could not get log probability for variant ploidy state %s.",
                    variantPloidyState.toString()));
        }
    }

    private Map<VariantPloidyState, Double> normalize(final Map<VariantPloidyState, Double> unnormalizedLogProbabilityMassFunctionMap) {
        final List<VariantPloidyState> states = new ArrayList<>(unnormalizedLogProbabilityMassFunctionMap.keySet());
        final double[] log10Probabilities = states.stream()
                .mapToDouble(s -> MathUtils.logToLog10(unnormalizedLogProbabilityMassFunctionMap.get(s))).toArray();
        final double[] probabilities = MathUtils.normalizeFromLog10(log10Probabilities);
        final LinkedHashMap<VariantPloidyState, Double> logProbabilityMassFunctionMap = new LinkedHashMap<>();
        IntStream.range(0, unnormalizedLogProbabilityMassFunctionMap.size())
                .forEach(i -> logProbabilityMassFunctionMap.put(states.get(i), Math.log(probabilities[i])));
        return logProbabilityMassFunctionMap;
    }
}
