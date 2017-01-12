package org.broadinstitute.hellbender.tools.coveragemodel.germline;

import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.coveragemodel.*;
import org.broadinstitute.hellbender.tools.coveragemodel.interfaces.CopyRatioExpectationsCalculator;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.IntegerCopyNumberHiddenMarkovModel;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.IntegerCopyNumberTransitionProbabilityCacheCollection;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;
import org.broadinstitute.hellbender.utils.hmm.HiddenMarkovModel;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class performs the following tasks:
 *
 * <dl>
 *     <dt>
 *         Calculates copy ratio posterior expectations as an instance of {@link CopyRatioExpectations}
 *         for a given list of active targets and emission data </dt>
 *     <dt>
 *         Calculates copy ratio prior expectations on a list of targets. The result is still reported as
 *         an instance of {@link CopyRatioExpectations}
 *     </dt>
 *     <dt>
 *         Generates forward-backward and Viterbi results as an instance of {@link CopyRatioHiddenMarkovModelResults}
 *     </dt>
 * </dl>
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IntegerCopyNumberExpectationsCalculator implements
        CopyRatioExpectationsCalculator<CoverageModelCopyRatioEmissionData, IntegerCopyNumberState>, Serializable {

    private static final long serialVersionUID = -6228432163615117934L;

    /**
     * This is for debugging -- disable in the future for performance gains
     */
    private final static boolean CHECK_FOR_NANS = false;

    private final Set<String> sexGenotypesSet;

    /**
     * A map from each sex genotype to a corresponding HMM
     */
    private final Map<String, IntegerCopyNumberHiddenMarkovModel<CoverageModelCopyRatioEmissionData>> hmm;

    /**
     * Integer copy number transition probability and prior cache
     */
    private final IntegerCopyNumberTransitionProbabilityCacheCollection cache;

    public IntegerCopyNumberExpectationsCalculator(
            @Nonnull final IntegerCopyNumberTransitionProbabilityCacheCollection cache,
            final int readCountThresholdPoissonSwitch) {
        this.cache = Utils.nonNull(cache, "The integer copy number transition probability cache" +
                " collection must be non-null");
        /* create an integer copy number HMM for each sex genotype in the collection */
        sexGenotypesSet = cache.getSexGenotypes();
        final CoverageModelCopyRatioEmissionProbabilityCalculator emissionProbabilityCalculator =
                new CoverageModelCopyRatioEmissionProbabilityCalculator(readCountThresholdPoissonSwitch);
        hmm = new HashMap<>();
        sexGenotypesSet.forEach(sexGenotype -> hmm.put(sexGenotype,
                new IntegerCopyNumberHiddenMarkovModel<>(emissionProbabilityCalculator, cache, sexGenotype)));
    }

    @Override
    public CopyRatioExpectations getCopyRatioPosteriorExpectations(
            @Nonnull final CopyRatioCallingMetadata copyRatioCallingMetadata,
            @Nonnull final List<Target> targetList,
            @Nonnull final List<CoverageModelCopyRatioEmissionData> emissionData) {
        verifyArgsForPosteriorExpectations(copyRatioCallingMetadata, targetList, emissionData);

        /* choose the HMM */
        final IntegerCopyNumberHiddenMarkovModel<CoverageModelCopyRatioEmissionData> genotypeSpecificHMM =
                hmm.get(copyRatioCallingMetadata.getSampleSexGenotypeData().getSexGenotype());

        /* supplement the emission data with sample metadata */
        emissionData.forEach(datum -> datum.setCopyRatioCallingMetadata(copyRatioCallingMetadata));

        /* run the forward-backward algorithm */
        final ForwardBackwardAlgorithm.Result<CoverageModelCopyRatioEmissionData, Target, IntegerCopyNumberState> result =
                ForwardBackwardAlgorithm.apply(emissionData, targetList, genotypeSpecificHMM);

        /* indices of hidden states to include in calculating posterior expectations */
        final List<IntegerCopyNumberState> hiddenStates = genotypeSpecificHMM.hiddenStates();

        /* exponentiate log copy number posteriors on all targets */
        final List<double[]> hiddenStateProbabilities = IntStream.range(0, targetList.size())
                .mapToObj(ti -> hiddenStates.stream()
                        .mapToDouble(s -> FastMath.exp(result.logProbability(ti, s)))
                        .toArray())
                .collect(Collectors.toList());

        /* copy ratio corrected for mapping error */
        final List<double[]> logCopyRatiosWithMappingError = IntStream.range(0, targetList.size())
                .mapToObj(ti -> {
                    final CoverageModelCopyRatioEmissionData currentEmissionData = emissionData.get(ti);
                    final double err = currentEmissionData.getMappingErrorProbability();
                    return hiddenStates.stream()
                            .mapToDouble(s -> FastMath.log((1 - err) * s.getCopyNumber() +
                                    err * FastMath.exp(-currentEmissionData.getMu())))
                            .toArray();
                })
                .collect(Collectors.toList());

        if (CHECK_FOR_NANS) {
            final int[] badTargets = IntStream.range(0, targetList.size())
                    .filter(ti -> Arrays.stream(hiddenStateProbabilities.get(ti))
                            .anyMatch(p -> Double.isNaN(p) || Double.isInfinite(p))).toArray();
            if (badTargets.length > 0) {
                throw new RuntimeException("Some of the copy ratio posterior probabilities are ill-defined; targets: " +
                        Arrays.stream(badTargets).mapToObj(String::valueOf).collect(Collectors.joining(", ", "[", "]")));
            }
        }

        /* calculate copy ratio posterior mean and variance */
        final int[] includedHiddenStatesIndices = IntStream.range(0, hiddenStates.size()).toArray();
        final double[] logCopyRatioPosteriorMeans = IntStream.range(0, targetList.size())
                .mapToDouble(ti -> calculateMeanDiscreteStates(
                        includedHiddenStatesIndices, hiddenStateProbabilities.get(ti),
                        logCopyRatiosWithMappingError.get(ti)))
                .toArray();
        final double[] logCopyRatioPosteriorVariances = IntStream.range(0, targetList.size())
                .mapToDouble(ti -> calculateMeanDiscreteStates(
                        includedHiddenStatesIndices, hiddenStateProbabilities.get(ti),
                        Arrays.stream(logCopyRatiosWithMappingError.get(ti)).map(x -> x * x).toArray()) -
                            FastMath.pow(logCopyRatioPosteriorMeans[ti], 2))
                .toArray();

        if (CHECK_FOR_NANS) {
            final int[] badTargets = IntStream.range(0, targetList.size())
                    .filter(ti -> Double.isNaN(logCopyRatioPosteriorMeans[ti]) ||
                            Double.isNaN(logCopyRatioPosteriorVariances[ti])).toArray();
            if (badTargets.length > 0) {
                throw new RuntimeException("Some of the copy ratio posterior expectations are ill-defined; targets: " +
                        Arrays.stream(badTargets).mapToObj(String::valueOf)
                                .collect(Collectors.joining(", ", "[", "]")));
            }
        }

        /* calculate chain posterior log probability */
        final double logChainPosteriorProbability = result.logChainPosteriorProbability();

        return new CopyRatioExpectations(logCopyRatioPosteriorMeans, logCopyRatioPosteriorVariances,
                logChainPosteriorProbability);
    }

    @Override
    public CopyRatioExpectations getCopyRatioPriorExpectations(
            @Nonnull final CopyRatioCallingMetadata copyRatioCallingMetadata,
            @Nonnull final List<Target> targetsList) {
        verifyArgsForPriorExpectations(copyRatioCallingMetadata, targetsList);

        /* choose the HMM */
        final IntegerCopyNumberHiddenMarkovModel<CoverageModelCopyRatioEmissionData> genotypeSpecificHMM =
                hmm.get(copyRatioCallingMetadata.getSampleSexGenotypeData().getSexGenotype());

        /* indices of hidden states to include in calculating posterior expectations */
        final List<IntegerCopyNumberState> hiddenStates = genotypeSpecificHMM.hiddenStates();

        /* exponentiate log copy number priors on all targets */
        final List<double[]> hiddenStatePriorProbabilities = targetsList.stream()
                .map(target -> hiddenStates.stream()
                        .mapToDouble(state -> FastMath.exp(genotypeSpecificHMM.logPriorProbability(state, target)))
                        .toArray())
                .collect(Collectors.toList());

        /* copy ratio corrected for mapping error */
        final double err = copyRatioCallingMetadata.getSampleAverageMappingErrorProbability();
        final List<double[]> logCopyRatiosWithMappingError = IntStream.range(0, targetsList.size())
                .mapToObj(ti -> hiddenStates.stream()
                            .mapToDouble(s -> FastMath.log((1 - err) * s.getCopyNumber() + err))
                            .toArray())
                .collect(Collectors.toList());

        if (CHECK_FOR_NANS) {
            final int[] badTargets = IntStream.range(0, targetsList.size())
                    .filter(ti -> Arrays.stream(hiddenStatePriorProbabilities.get(ti))
                            .anyMatch(p -> Double.isNaN(p) || Double.isInfinite(p))).toArray();
            if (badTargets.length > 0) {
                throw new RuntimeException("Some of the copy ratio prior probabilities are ill-defined; targets: " +
                        Arrays.stream(badTargets).mapToObj(String::valueOf).collect(Collectors.joining(", ", "[", "]")));
            }
        }

        /* calculate copy ratio posterior mean and variance */
        final int[] includedHiddenStatesIndices = IntStream.range(0, hiddenStates.size()).toArray();
        final double[] logCopyRatioPriorMeans = IntStream.range(0, targetsList.size())
                .mapToDouble(ti -> calculateMeanDiscreteStates(
                        includedHiddenStatesIndices, hiddenStatePriorProbabilities.get(ti),
                        logCopyRatiosWithMappingError.get(ti)))
                .toArray();
        final double[] logCopyRatioPriorVariances = IntStream.range(0, targetsList.size())
                .mapToDouble(ti -> calculateMeanDiscreteStates(
                        includedHiddenStatesIndices, hiddenStatePriorProbabilities.get(ti),
                        Arrays.stream(logCopyRatiosWithMappingError.get(ti)).map(x -> x * x).toArray()) -
                        FastMath.pow(logCopyRatioPriorMeans[ti], 2))
                .toArray();

        if (CHECK_FOR_NANS) {
            final int[] badTargets = IntStream.range(0, targetsList.size())
                    .filter(ti -> Double.isNaN(logCopyRatioPriorMeans[ti]) ||
                            Double.isNaN(logCopyRatioPriorVariances[ti])).toArray();
            if (badTargets.length > 0) {
                throw new RuntimeException("Some of the copy ratio prior expectations are ill-defined; targets: " +
                        Arrays.stream(badTargets).mapToObj(String::valueOf).collect(Collectors.joining(", ", "[", "]")));
            }
        }

        /* TODO -- this is slow and unnecessary anyway */
        /* calculate chain posterior log probability */
//        final double logChainPriorProbability = calculateLogChainPriorProbability(genotypeSpecificHMM, targetsList);

        return new CopyRatioExpectations(logCopyRatioPriorMeans, logCopyRatioPriorVariances, 0.0);
    }

    @Override
    public CopyRatioHiddenMarkovModelResults<CoverageModelCopyRatioEmissionData,
            IntegerCopyNumberState> getCopyRatioHiddenMarkovModelResults(@Nonnull final CopyRatioCallingMetadata copyRatioCallingMetadata,
                                                                         @Nonnull final List<Target> activeTargets,
                                                                         @Nonnull final List<CoverageModelCopyRatioEmissionData> emissionData) {
        verifyArgsForHiddenMarkovModelResults(copyRatioCallingMetadata, activeTargets, emissionData);

        /* choose the HMM */
        final IntegerCopyNumberHiddenMarkovModel<CoverageModelCopyRatioEmissionData> genotypeSpecificHMM =
                hmm.get(copyRatioCallingMetadata.getSampleSexGenotypeData().getSexGenotype());

        /* supplement the emission data with sample metadata */
        emissionData.forEach(datum -> datum.setCopyRatioCallingMetadata(copyRatioCallingMetadata));

        /* run the forward-backward algorithm */
        final ForwardBackwardAlgorithm.Result<CoverageModelCopyRatioEmissionData, Target, IntegerCopyNumberState> fbResult =
                ForwardBackwardAlgorithm.apply(emissionData, activeTargets, genotypeSpecificHMM);

        /* run the Viterbi algorithm */
        final List<IntegerCopyNumberState> viterbiResult = ViterbiAlgorithm.apply(emissionData, activeTargets,
                genotypeSpecificHMM);

        /* clear caches */
        genotypeSpecificHMM.clearCaches();

        return new CopyRatioHiddenMarkovModelResults<>(copyRatioCallingMetadata, fbResult, viterbiResult);
    }

    @Override
    public void initializeCaches(@Nonnull final List<Target> allTargets) {
        Utils.nonNull(allTargets, "The list of targets must be non-null");
        for (final String sexGenotype : cache.getSexGenotypes()) {
            IntStream.range(0, allTargets.size() - 1).forEach(firstTargetIndex ->
                    cache.cacheTransitionMatrix(
                            (int)Target.calculateDistance(allTargets.get(firstTargetIndex),
                                    allTargets.get(firstTargetIndex + 1),
                                    IntegerCopyNumberHiddenMarkovModel.DEFAULT_DISTANCE_BETWEEN_TARGETS),
                            sexGenotype, allTargets.get(firstTargetIndex).getContig()));
        }
    }

    @Override
    public void clearCaches() {
        cache.clearCaches();
    }

    /**
     * Calculates the mean of a function evaluated on a set of discrete states
     * (No checks are performed for speed)
     *
     * @param indices an array of indices for the states to include
     * @param pdf the probability distribution function (does not need to be normalized
     * @param func a function evaluated over the discrete states
     * @return a double value
     */
    private static double calculateMeanDiscreteStates(@Nonnull final int[] indices,
                                                      @Nonnull final double[] pdf,
                                                      @Nonnull final double[] func) {
        double prod = 0, sum = 0;
        for (final int i : indices) {
            prod += pdf[i] * func[i];
            sum += pdf[i];
        }
        return prod / sum;
    }

    /**
     * TODO
     *
     * @param hmm
     * @param positions
     * @param <D>
     * @param <T>
     * @param <S>
     * @return
     */
    private static <D, T, S> double calculateLogChainPriorProbability(final HiddenMarkovModel<D, T, S> hmm,
                                                                      final List<T> positions) {
        final List<S> states = hmm.hiddenStates();
        if (positions.isEmpty() || states.isEmpty()) {
            return 0;
        }
        final int numStates = states.size();
        final int numPositions = positions.size();
        final double[] logPriorProbabilities = states.stream()
                .mapToDouble(state -> hmm.logPriorProbability(state, positions.get(0)))
                .toArray();
        final double[] priorProbabilities = Arrays.stream(logPriorProbabilities)
                .map(FastMath::exp)
                .toArray();

        /* contribution of the first state */
        double result = IntStream.range(0, numStates)
                .mapToDouble(stateIndex -> logPriorProbabilities[stateIndex] * priorProbabilities[stateIndex])
                .sum();

        /* contribution of the rest of the chain */
        if (numPositions > 1) {
            for (int positionIndex = 0; positionIndex < numPositions - 1; positionIndex++) {
                for (int i = 0; i < numStates; i++) { /* departure state */
                    for (int j = 0; j < numStates; j++) { /* destination state */
                        final double logTransitionProbability = hmm.logTransitionProbability(
                                states.get(i), positions.get(positionIndex),
                                states.get(j), positions.get(positionIndex + 1));
                        if (logTransitionProbability != Double.NEGATIVE_INFINITY) { /* NEGATIVE_INFINITY * 0 = 0 here */
                            result += logTransitionProbability * priorProbabilities[i] * priorProbabilities[j];
                        }
                    }
                }
            }
        }
        return result;
    }

    private void verifyArgsForHiddenMarkovModelResults(@Nonnull final CopyRatioCallingMetadata copyRatioCallingMetadata,
                                                       @Nonnull List<Target> targetsList,
                                                       @Nonnull List<CoverageModelCopyRatioEmissionData> data) {
        Utils.nonNull(copyRatioCallingMetadata, "Sample metadata must be non-null");
        Utils.nonNull(targetsList, "List of targets must be non-null");
        Utils.nonNull(data, "List of copy ratio emission data must be non-null");
        Utils.validateArg(sexGenotypesSet.contains(copyRatioCallingMetadata.getSampleSexGenotypeData().getSexGenotype()),
                String.format("The sex genotype \"%s\" does not exist in the transition cache collection",
                        copyRatioCallingMetadata.getSampleSexGenotypeData().getSexGenotype()));
        Utils.validateArg(targetsList.size() >= 2, "At least two active targets are required");
        Utils.validateArg(targetsList.size() == data.size(), String.format("Number of active targets (%d) does not" +
                " must match the length of emission data list (%d)", targetsList.size(), data.size()));
        targetsList.forEach(t -> Utils.nonNull(t, "Some of the active targets are null"));
        data.forEach(d -> Utils.nonNull(d, "Some of the emission data entries are null"));
        targetsList.stream()
                .collect(Collectors.groupingBy(Target::getContig)) /* group by contig */
                .values()
                .forEach(targetList -> Utils.validateArg(IntStream.range(0, targetList.size() - 1)
                                .allMatch(i -> targetList.get(i + 1).getStart() - targetList.get(i).getEnd() >= 0),
                        "The list of active targets must be coordinate sorted and non-overlapping" +
                                " (except for their endpoints)"));
    }

    private void verifyArgsForPosteriorExpectations(@Nonnull final CopyRatioCallingMetadata copyRatioCallingMetadata,
                                                    @Nonnull final List<Target> targetList,
                                                    @Nonnull final List<CoverageModelCopyRatioEmissionData> emissionData) {
        verifyArgsForHiddenMarkovModelResults(copyRatioCallingMetadata, targetList, emissionData);
    }

    private void verifyArgsForPriorExpectations(@Nonnull final CopyRatioCallingMetadata copyRatioCallingMetadata,
                                                @Nonnull List<Target> targetsList) {
        Utils.nonNull(copyRatioCallingMetadata, "Sample metadata must be non-null");
        Utils.nonNull(targetsList, "List of targets must be non-null");
        Utils.validateArg(sexGenotypesSet.contains(copyRatioCallingMetadata.getSampleSexGenotypeData().getSexGenotype()),
                String.format("The sex genotype \"%s\" does not exist in the transition cache collection",
                        copyRatioCallingMetadata.getSampleSexGenotypeData().getSexGenotype()));
        Utils.validateArg(targetsList.size() >= 2, "At least two active targets are required");
        targetsList.forEach(t -> Utils.nonNull(t, "Some of the active targets are null"));
        targetsList.stream()
                .collect(Collectors.groupingBy(Target::getContig)) /* group by contig */
                .values()
                .forEach(targetList -> Utils.validateArg(IntStream.range(0, targetList.size() - 1)
                                .allMatch(i -> targetList.get(i + 1).getStart() - targetList.get(i).getEnd() >= 0),
                        "The list of active targets must be coordinate sorted and non-overlapping" +
                                " (except for their endpoints)"));
    }
}