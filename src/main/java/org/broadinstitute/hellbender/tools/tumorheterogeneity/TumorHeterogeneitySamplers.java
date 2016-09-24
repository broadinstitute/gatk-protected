package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Gamma;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.mcmc.ParameterSampler;
import org.broadinstitute.hellbender.utils.mcmc.SliceSampler;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class TumorHeterogeneitySamplers {
    private static final double EPSILON = 1E-10;
    private static final double EPSILON_TIMES_TWO = 2. * EPSILON;

    private static final Logger logger = LogManager.getLogger(TumorHeterogeneitySamplers.class);

    private TumorHeterogeneitySamplers() {}

    protected static final class ConcentrationSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private final double concentrationMin;
        private final double concentrationMax;
        private final double concentrationSliceSamplingWidth;

        public ConcentrationSampler(final double concentrationMin, final double concentrationMax, final double concentrationSliceSamplingWidth) {
            this.concentrationMin = concentrationMin;
            this.concentrationMax = concentrationMax;
            this.concentrationSliceSamplingWidth = concentrationSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            logger.info(String.format("Ploidy of current state: %.6f", state.calculatePopulationAndGenomicAveragedPloidy(data)));
            logger.debug("Sampling concentration.");
            final int numPopulations = state.numPopulations();
            final Function<Double, Double> logConditionalPDF = newConcentration -> {
                final double populationFractionsTerm = IntStream.range(0, numPopulations)
                        .mapToDouble(i -> (newConcentration - 1) * Math.log(state.populationFraction(i) + EPSILON)).sum();
                return (state.priors().concentrationPriorAlpha() - 1.) * Math.log(newConcentration) - state.priors().concentrationPriorBeta() * newConcentration +
                        Gamma.logGamma(newConcentration * numPopulations) - numPopulations * Gamma.logGamma(newConcentration) + populationFractionsTerm;
            };
            return new SliceSampler(rng, logConditionalPDF, concentrationMin, concentrationMax, concentrationSliceSamplingWidth).sample(state.concentration());
        }
    }

    protected static final class PopulationFractionsSampler implements ParameterSampler<TumorHeterogeneityState.PopulationFractions, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        public PopulationFractionsSampler() {
        }

        @Override
        public TumorHeterogeneityState.PopulationFractions sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            logger.debug("Sampling population fractions.");
            //sampling from Dirichlet(alpha_vec) is equivalent to sampling from individual Gamma(alpha_vec_i, 1) distributions and normalizing
            final double[] unnormalizedPopulationFractions = IntStream.range(0, state.numPopulations()).boxed()
                    .map(state::populationCount)
                    .mapToDouble(c -> new GammaDistribution(rng, state.concentration() + c.doubleValue(), 1.).sample()).toArray();
            final List<Double> populationFractions = Doubles.asList(MathUtils.normalizeFromRealSpace(unnormalizedPopulationFractions));
            return new TumorHeterogeneityState.PopulationFractions(populationFractions);
        }
    }

    protected static final class PopulationIndicatorsSampler implements ParameterSampler<TumorHeterogeneityState.PopulationIndicators, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private final Random rnd = new Random(1845);
        private final double singleCellPopulationFraction;
        private final List<Integer> populationIndices;


        public PopulationIndicatorsSampler(final int numCells, final int numPopulations) {
            singleCellPopulationFraction = 1. / numCells;
            populationIndices = Collections.unmodifiableList(IntStream.range(0, numPopulations).boxed().collect(Collectors.toList()));
        }

        @Override
        public TumorHeterogeneityState.PopulationIndicators sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            logger.debug("Sampling population indicators.");
            final List<Integer> populationIndicators = new ArrayList<>(Collections.nCopies(state.numCells(), 0));
            final List<Integer> shuffledCellIndices = IntStream.range(0, state.numCells()).boxed().collect(Collectors.toList());
            Collections.shuffle(shuffledCellIndices, rnd);
            for (int cellIndex : shuffledCellIndices) {
                final int currentPopulationIndex = state.populationIndex(cellIndex);

                final double invariantPloidyTerm = calculatePopulationAndGenomicAveragedPloidyWithExclusion(state, data, currentPopulationIndex, singleCellPopulationFraction);

                final double[] log10Probabilities = new double[state.numPopulations()];
                for (int populationIndex = 0; populationIndex < state.numPopulations(); populationIndex++) {
                    double logDensity = Math.log(state.populationFraction(populationIndex) + EPSILON);    //prior
                    for (int segmentIndex = 0; segmentIndex < state.numSegments(); segmentIndex++) {
                        final double segmentFractionalLength = state.calculateFractionalLength(data, segmentIndex);
                        final double invariantMAlleleCopyNumberTerm = calculatePopulationAveragedMAlleleCopyNumberWithExclusion(state, data, segmentIndex, currentPopulationIndex, singleCellPopulationFraction);
                        final double invariantNAlleleCopyNumberTerm = calculatePopulationAveragedNAlleleCopyNumberWithExclusion(state, data, segmentIndex, currentPopulationIndex, singleCellPopulationFraction);

                        final boolean isVariant = state.isVariant(populationIndex, segmentIndex);
                        final PloidyState ploidyState = isVariant ? state.variantPloidyState(populationIndex, segmentIndex) : state.priors().normalPloidyState();
                        logDensity += calculateLogDensityFromInvariantTerms(data, invariantPloidyTerm, invariantMAlleleCopyNumberTerm, invariantNAlleleCopyNumberTerm,
                                segmentIndex, singleCellPopulationFraction, segmentFractionalLength, ploidyState);
                    }
                    log10Probabilities[populationIndex] = MathUtils.logToLog10(logDensity);
                }
                final double[] probabilities = MathUtils.normalizeFromLog10(log10Probabilities);
                final int populationIndex = GATKProtectedMathUtils.randomSelect(populationIndices, i -> probabilities[i], rng);
                populationIndicators.set(cellIndex, populationIndex);

                //update the state as a side effect
                state.setPopulationIndex(cellIndex, populationIndex);
            }
            return new TumorHeterogeneityState.PopulationIndicators(populationIndicators);
        }

        private static double calculatePopulationAndGenomicAveragedPloidyWithExclusion(final TumorHeterogeneityState state, final TumorHeterogeneityData data, final int populationIndexToExclude, final double populationFractionToExclude) {
            final double excludedPloidy = populationFractionToExclude * IntStream.range(0, state.numSegments())
                    .mapToDouble(i -> state.calculateFractionalLength(data, i) * state.calculateCopyNumberFunction(i, populationIndexToExclude, PloidyState::total))
                    .sum();
            return state.calculatePopulationAndGenomicAveragedPloidy(data) - excludedPloidy;
        }

        private static double calculatePopulationAveragedMAlleleCopyNumberWithExclusion(final TumorHeterogeneityState state, final TumorHeterogeneityData data, final int segmentIndex, final int populationIndexToExclude, final double populationFractionToExclude) {
            final double excludedMAlleleCopyNumber = state.calculateCopyNumberFunction(segmentIndex, populationIndexToExclude, PloidyState::m) * populationFractionToExclude;
            return state.calculatePopulationAveragedMAlleleCopyNumber(data, segmentIndex) - excludedMAlleleCopyNumber;
        }

        private static double calculatePopulationAveragedNAlleleCopyNumberWithExclusion(final TumorHeterogeneityState state, final TumorHeterogeneityData data, final int segmentIndex, final int populationIndexToExclude, final double populationFractionToExclude) {
            final double excludedNAlleleCopyNumber = state.calculateCopyNumberFunction(segmentIndex, populationIndexToExclude, PloidyState::n) * populationFractionToExclude;
            return state.calculatePopulationAveragedNAlleleCopyNumber(data, segmentIndex) - excludedNAlleleCopyNumber;
        }
    }

    /**
     * Samples genomic profiles for a collection of variant populations.
     */
    protected static final class VariantProfileCollectionSampler implements ParameterSampler<TumorHeterogeneityState.VariantProfileCollection, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private final List<VariantProfileSampler> variantProfileSamplers;

        public VariantProfileCollectionSampler(final int numVariantPopulations, final PloidyStatePrior variantPloidyStatePrior) {
            variantProfileSamplers = IntStream.range(0, numVariantPopulations).boxed().map(n -> new VariantProfileSampler(n, variantPloidyStatePrior)).collect(Collectors.toList());
        }

        public TumorHeterogeneityState.VariantProfileCollection sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            logger.debug("Sampling variant profile collection.");
            final List<TumorHeterogeneityState.VariantProfile> variantProfiles = variantProfileSamplers.stream().map(sampler -> sampler.sample(rng, state, data)).collect(Collectors.toList());
            return new TumorHeterogeneityState.VariantProfileCollection(variantProfiles);
        }
    }

    /**
     * Samples a genomic profile for a single variant population.
     */
    protected static final class VariantProfileSampler implements ParameterSampler<TumorHeterogeneityState.VariantProfile, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private int populationIndex;
        private final VariantSegmentFractionSampler variantSegmentFractionSampler;
        private final VariantIndicatorsSampler variantIndicatorsSampler;
        private final VariantPloidyStateIndicatorsSampler variantPloidyStateIndicatorsSampler;

        public VariantProfileSampler(final int populationIndex, final PloidyStatePrior variantPloidyStatePrior) {
            this.populationIndex = populationIndex;
            variantSegmentFractionSampler = new VariantSegmentFractionSampler();
            variantIndicatorsSampler = new VariantIndicatorsSampler();
            variantPloidyStateIndicatorsSampler = new VariantPloidyStateIndicatorsSampler(variantPloidyStatePrior);
        }

        public TumorHeterogeneityState.VariantProfile sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            logger.debug("Sampling variant profile for population " + populationIndex);
            final double variantSegmentFraction = variantSegmentFractionSampler.sample(rng, state, data);
            final TumorHeterogeneityState.VariantProfile.VariantIndicators variantIndicators = variantIndicatorsSampler.sample(rng, state, data);
            final TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators variantPloidyStateIndicators = variantPloidyStateIndicatorsSampler.sample(rng, state, data);
            return new TumorHeterogeneityState.VariantProfile(variantSegmentFraction, variantIndicators, variantPloidyStateIndicators);
        }

        /**
         * Samples variant-segment fraction for this variant population.
         */
        private final class VariantSegmentFractionSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
            private VariantSegmentFractionSampler() {}

            @Override
            public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
                logger.debug("Sampling variant-segment fraction for population " + populationIndex);
                final int numSegmentsVariant = (int) IntStream.range(0, state.numSegments()).filter(i -> state.isVariant(populationIndex, i)).count();
                return new BetaDistribution(rng,
                        state.priors().variantSegmentFractionPriorAlpha() + numSegmentsVariant,
                        state.priors().variantSegmentFractionPriorBeta() + data.numSegments() - numSegmentsVariant).sample();
            }
        }

        /**
         * Samples per-segment variant indicators for this variant population (these indicate whether or not population is variant in each segment).
         * Updates the variant indicators held by the {@link TumorHeterogeneityState}.
         */
        private final class VariantIndicatorsSampler implements ParameterSampler<TumorHeterogeneityState.VariantProfile.VariantIndicators, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
            private final Random rnd = new Random(1425);

            private VariantIndicatorsSampler() {}

            @Override
            public TumorHeterogeneityState.VariantProfile.VariantIndicators sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
                logger.debug("Sampling variant indicators for population " + populationIndex);
                final List<Boolean> variantIndicators = new ArrayList<>(Collections.nCopies(state.numSegments(), false));
                final List<Integer> shuffledSegmentIndices = IntStream.range(0, state.numSegments()).boxed().collect(Collectors.toList());
                Collections.shuffle(shuffledSegmentIndices, rnd);
                for (int segmentIndex : shuffledSegmentIndices) {
                    final double populationFraction = state.calculatePopulationFractionFromCounts(populationIndex);
                    final double segmentFractionalLength = state.calculateFractionalLength(data, segmentIndex);
                    final double invariantPloidyTerm = calculatePopulationAndGenomicAveragedPloidyWithExclusion(state, data, segmentIndex, populationIndex);
                    final double invariantMAlleleCopyNumberTerm = calculatePopulationAveragedMAlleleCopyNumberWithExclusion(state, data, segmentIndex, populationIndex);
                    final double invariantNAlleleCopyNumberTerm = calculatePopulationAveragedNAlleleCopyNumberWithExclusion(state, data, segmentIndex, populationIndex);

                    //approximation: ignore coupling of copy-ratio posteriors in different segments due to ploidy term

                    final double variantSegmentFraction = Math.max(EPSILON, Math.min(state.variantSegmentFraction(populationIndex), 1. - EPSILON));

                    //calculate unnormalized probability for isVariant = true
                    final double logDensityVariant = Math.log(variantSegmentFraction) +
                            calculateLogDensityFromInvariantTerms(
                                    data, invariantPloidyTerm, invariantMAlleleCopyNumberTerm, invariantNAlleleCopyNumberTerm,
                                    segmentIndex, populationFraction, segmentFractionalLength, state.variantPloidyState(populationIndex, segmentIndex));
                    //calculate unnormalized probability for isVariant = false
                    final double logDensityNormal = Math.log(1. - variantSegmentFraction) +
                            calculateLogDensityFromInvariantTerms(
                                    data, invariantPloidyTerm, invariantMAlleleCopyNumberTerm, invariantNAlleleCopyNumberTerm,
                                    segmentIndex, populationFraction, segmentFractionalLength, state.priors().normalPloidyState());

                    final double[] log10Probabilities = new double[]{MathUtils.logToLog10(logDensityVariant), MathUtils.logToLog10(logDensityNormal)};
                    final double isVariantProbability = MathUtils.normalizeFromLog10(log10Probabilities)[0];
                    final boolean isVariant = rng.nextDouble() < isVariantProbability;
                    variantIndicators.set(segmentIndex, isVariant);

                    //update the state as a side effect
                    state.setIsVariant(populationIndex, segmentIndex, isVariant);
                }
                return new TumorHeterogeneityState.VariantProfile.VariantIndicators(variantIndicators);
            }
        }

        /**
         * Samples per-segment variant-ploidy-state indicators for this variant population (these indicate the ploidy state if the population is variant in a segment).
         * Updates the variant-ploidy-state indicators held by the {@link TumorHeterogeneityState}.
         */
        private final class VariantPloidyStateIndicatorsSampler implements ParameterSampler<TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
            private final Random rnd = new Random(4615);
            private final List<Integer> variantPloidyStateIndices;
            private final double[] variantPloidyStatePriorLog10Probabilities;

            private VariantPloidyStateIndicatorsSampler(final PloidyStatePrior variantPloidyStatePrior) {
                final int numVariantPloidyStates = variantPloidyStatePrior.numPloidyStates();
                variantPloidyStateIndices = Collections.unmodifiableList(IntStream.range(0, numVariantPloidyStates).boxed().collect(Collectors.toList()));
                variantPloidyStatePriorLog10Probabilities = variantPloidyStatePrior.ploidyStates().stream()
                        .mapToDouble(vps -> MathUtils.logToLog10(variantPloidyStatePrior.logProbability(vps)))
                        .toArray();
            }

            @Override
            public TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
                logger.debug("Sampling variant-ploidy-state indicators for population " + populationIndex);
                final List<Integer> variantPloidyStateIndicators = new ArrayList<>(Collections.nCopies(state.numSegments(), 0));
                final List<Integer> shuffledSegmentIndices = IntStream.range(0, state.numSegments()).boxed().collect(Collectors.toList());
                Collections.shuffle(shuffledSegmentIndices, rnd);
                for (int segmentIndex : shuffledSegmentIndices) {
                    final boolean isVariant = state.isVariant(populationIndex, segmentIndex);
                    final double segmentFractionalLength = state.calculateFractionalLength(data, segmentIndex);
                    final double populationFraction = state.calculatePopulationFractionFromCounts(populationIndex);
                    final double invariantPloidyTerm = calculatePopulationAndGenomicAveragedPloidyWithExclusion(state, data, segmentIndex, populationIndex);
                    final double invariantMAlleleCopyNumberTerm = calculatePopulationAveragedMAlleleCopyNumberWithExclusion(state, data, segmentIndex, populationIndex);
                    final double invariantNAlleleCopyNumberTerm = calculatePopulationAveragedNAlleleCopyNumberWithExclusion(state, data, segmentIndex, populationIndex);

                    //approximation: ignore coupling of copy-ratio posteriors in different segments due to ploidy term

                    //calculate unnormalized probabilities for all variant ploidy states
                    final List<PloidyState> variantPloidyStates = state.priors().variantPloidyStatePrior().ploidyStates();
                    final int si = segmentIndex;
                    final double[] log10Probabilities = isVariant ?
                            variantPloidyStateIndices.stream()
                                    .mapToDouble(i -> variantPloidyStatePriorLog10Probabilities[i] +
                                            MathUtils.logToLog10(calculateLogDensityFromInvariantTerms(
                                                    data, invariantPloidyTerm, invariantMAlleleCopyNumberTerm, invariantNAlleleCopyNumberTerm,
                                                    si, populationFraction, segmentFractionalLength, variantPloidyStates.get(i))))
                                    .toArray() :
                            variantPloidyStatePriorLog10Probabilities;
                    final double[] probabilities = MathUtils.normalizeFromLog10(log10Probabilities);
                    final int variantPloidyStateIndex = GATKProtectedMathUtils.randomSelect(variantPloidyStateIndices, i -> probabilities[i], rng);
                    variantPloidyStateIndicators.set(segmentIndex, variantPloidyStateIndex);

                    //update the state as a side effect
                    state.setVariantPloidyStateIndex(populationIndex, segmentIndex, variantPloidyStateIndex);
                }
                return new TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators(variantPloidyStateIndicators);
            }
        }

        private static double calculatePopulationAndGenomicAveragedPloidyWithExclusion(final TumorHeterogeneityState state, final TumorHeterogeneityData data, final int segmentIndexToExclude, final int populationIndexToExclude) {
            final double excludedPloidy = state.populationFraction(populationIndexToExclude) * state.calculateFractionalLength(data, segmentIndexToExclude) * state.calculateCopyNumberFunction(segmentIndexToExclude, populationIndexToExclude, PloidyState::total);
            return state.calculatePopulationAndGenomicAveragedPloidy(data) - excludedPloidy;
        }

        private static double calculatePopulationAveragedMAlleleCopyNumberWithExclusion(final TumorHeterogeneityState state, final TumorHeterogeneityData data, final int segmentIndex, final int populationIndexToExclude) {
            return state.calculatePopulationAveragedCopyNumberFunction(data, segmentIndex, populationIndexToExclude, PloidyState::m, false);
        }

        private static double calculatePopulationAveragedNAlleleCopyNumberWithExclusion(final TumorHeterogeneityState state, final TumorHeterogeneityData data, final int segmentIndex, final int populationIndexToExclude) {
            return state.calculatePopulationAveragedCopyNumberFunction(data, segmentIndex, populationIndexToExclude, PloidyState::n, false);
        }
    }

    private static double calculateLogDensityFromInvariantTerms(final TumorHeterogeneityData data,
                                                                final double invariantPloidyTerm,
                                                                final double invariantMAlleleCopyNumberTerm,
                                                                final double invariantNAlleleCopyNumberTerm,
                                                                final int segmentIndex,
                                                                final double populationFraction,
                                                                final double segmentFractionalLength,
                                                                final PloidyState ploidyState) {
        final double ploidy = invariantPloidyTerm + populationFraction * segmentFractionalLength * ploidyState.total();
        final double copyRatio = (invariantMAlleleCopyNumberTerm + invariantNAlleleCopyNumberTerm + populationFraction * ploidyState.total()) / (ploidy + EPSILON);
        final double minorAlleleFraction = calculateMinorAlleleFraction(
                invariantMAlleleCopyNumberTerm + populationFraction * ploidyState.m(),
                invariantNAlleleCopyNumberTerm + populationFraction * ploidyState.n());
        return data.logDensity(segmentIndex, copyRatio, minorAlleleFraction);
    }

    private static double calculateMinorAlleleFraction(final double m, final double n) {
        return Math.min(m + EPSILON, n + EPSILON) / (m + n + EPSILON_TIMES_TWO);
    }
}
