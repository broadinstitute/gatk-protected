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

import java.util.*;
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

    protected static final class DoMetropolisStepSampler implements ParameterSampler<Boolean, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        protected DoMetropolisStepSampler() {}

        @Override
        public Boolean sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            final boolean doMetropolisStep = rng.nextDouble() < state.priors().metropolisIterationFraction();
            if (doMetropolisStep) {
                logger.info("Performing Metropolis step.");
                //if Metropolis step is accepted, update the relevant parameters blocks in the state;
                //when these blocks are sampled during this iteration, these updated values will be returned
                performMetropolisStep(rng, state, data);
            }
            return doMetropolisStep;
        }

        private void performMetropolisStep(final RandomGenerator rng,
                                           final TumorHeterogeneityState state,
                                           final TumorHeterogeneityData data) {
            final TumorHeterogeneityState proposedState = TumorHeterogeneityStateInitializationUtils.proposeState(rng, state, data);
            final double proposedLogPosterior = calculateLogPosterior(proposedState, data);
            final double currentLogPosterior = calculateLogPosterior(state, data);
            final double acceptanceProbability = Math.min(1., Math.exp(proposedLogPosterior - currentLogPosterior));
            logger.debug("Log posterior of current state: " + currentLogPosterior);
            logger.debug("Log posterior of proposed state: " + proposedLogPosterior);
            if (rng.nextDouble() < acceptanceProbability) {
                logger.info("Proposed state accepted.");
                state.set(proposedState);
            }
        }

        private static double calculateLogPosterior(final TumorHeterogeneityState state,
                                                    final TumorHeterogeneityData data) {
            final int numPopulations = state.numPopulations();
            final int numSegments = data.numSegments();

            //concentration prior
            final double concentrationPriorAlpha = state.priors().concentrationPriorAlpha();
            final double concentrationPriorBeta = state.priors().concentrationPriorBeta();
            final double concentration = state.concentration();
            final double logPriorConcentration =
                    concentrationPriorAlpha * Math.log(concentrationPriorBeta + EPSILON)
                            + (concentrationPriorAlpha - 1.) * Math.log(concentration)
                            - concentrationPriorBeta * concentration
                            - Gamma.logGamma(concentrationPriorAlpha);

            //population-fractions prior
            final double logPriorPopulationFractionsSum = IntStream.range(0, numPopulations)
                    .mapToDouble(i -> (concentration - 1.) * Math.log(state.populationFraction(i) + EPSILON))
                    .sum();
            final double logPriorPopulationFractions =
                    Gamma.logGamma(concentration * numPopulations)
                            - numPopulations * Gamma.logGamma(concentration)
                            + logPriorPopulationFractionsSum;
            //population-fractions likelihood
            final double logLikelihoodPopulationFractions = IntStream.range(0, numPopulations)
                    .mapToDouble(i -> state.populationCount(i) * Math.log(state.populationFraction(i) + EPSILON))
                    .sum();

            //variant-profiles prior
            double logPriorVariantProfiles = 0.;
            for (int populationIndex = 0; populationIndex < numPopulations - 1; populationIndex++) {
                for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
                    final int ploidyStateIndex = state.ploidyStateIndex(populationIndex, segmentIndex);
                    final PloidyState ploidyState = state.priors().ploidyStatePrior().ploidyStates().get(ploidyStateIndex);
                    logPriorVariantProfiles += state.priors().ploidyStatePrior().logProbability(ploidyState);
                }
            }

            //copy-ratio--minor-allele-fraction likelihood
            double logLikelihoodSegments = 0.;
            final double ploidy = state.calculatePopulationAndGenomicAveragedPloidy(data);
            for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
                final double totalCopyNumber = state.calculatePopulationAveragedTotalCopyNumber(data, segmentIndex);
                final double mAlleleCopyNumber = state.calculatePopulationAveragedMAlleleCopyNumber(data, segmentIndex);
                final double nAlleleCopyNumber = state.calculatePopulationAveragedNAlleleCopyNumber(data, segmentIndex);
                final double copyRatio = totalCopyNumber / (ploidy + EPSILON);
                final double minorAlleleFraction = calculateMinorAlleleFraction(mAlleleCopyNumber, nAlleleCopyNumber);
                logLikelihoodSegments += data.logDensity(segmentIndex, copyRatio, minorAlleleFraction);
            }

            return logPriorConcentration +
                    logPriorPopulationFractions + logLikelihoodPopulationFractions +
                    logPriorVariantProfiles +
                    logLikelihoodSegments;
        }
    }

    protected static final class ConcentrationSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private final double concentrationMin;
        private final double concentrationMax;
        private final double concentrationSliceSamplingWidth;

        protected ConcentrationSampler(final double concentrationMin, final double concentrationMax, final double concentrationSliceSamplingWidth) {
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
        protected PopulationFractionsSampler() {
        }

        @Override
        public TumorHeterogeneityState.PopulationFractions sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            logger.debug("Sampling population fractions.");
            return state.doMetropolisStep() ? new TumorHeterogeneityState.PopulationFractions(state.populationFractions()) : sampleGibbs(rng, state);
        }

        private static TumorHeterogeneityState.PopulationFractions sampleGibbs(final RandomGenerator rng, final TumorHeterogeneityState state) {
            //sampling from Dirichlet(alpha_vec) is equivalent to sampling from individual Gamma(alpha_vec_i, 1) distributions and normalizing
            logger.debug("Current population counts: " + IntStream.range(0, state.numPopulations()).boxed()
                    .map(state::populationCount).collect(Collectors.toList()));
            final double[] unnormalizedPopulationFractions = IntStream.range(0, state.numPopulations()).boxed()
                    .mapToDouble(i -> new GammaDistribution(rng, state.concentration() + state.populationCount(i), 1.).sample()).toArray();
            final List<Double> populationFractions = Doubles.asList(MathUtils.normalizeFromRealSpace(unnormalizedPopulationFractions));
            logger.debug("Sampled population fractions: " + populationFractions);
            return new TumorHeterogeneityState.PopulationFractions(populationFractions);
        }
    }

    protected static final class PopulationIndicatorsSampler implements ParameterSampler<TumorHeterogeneityState.PopulationIndicators, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private final Random rnd = new Random(1845);
        private final double singleCellPopulationFraction;
        private final List<Integer> populationIndices;

        protected PopulationIndicatorsSampler(final int numCells, final int numPopulations) {
            singleCellPopulationFraction = 1. / numCells;
            populationIndices = Collections.unmodifiableList(IntStream.range(0, numPopulations).boxed().collect(Collectors.toList()));
        }

        @Override
        public TumorHeterogeneityState.PopulationIndicators sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            logger.debug("Sampling population indicators.");
            return state.doMetropolisStep() ? new TumorHeterogeneityState.PopulationIndicators(state.populationIndicators()) : sampleGibbs(rng, state, data);
        }

        private TumorHeterogeneityState.PopulationIndicators sampleGibbs(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            final List<Integer> populationIndicators = new ArrayList<>(Collections.nCopies(state.numCells(), 0));
            final List<Integer> shuffledCellIndices = IntStream.range(0, state.numCells()).boxed().collect(Collectors.toList());
            Collections.shuffle(shuffledCellIndices, rnd);
            for (int cellIndex : shuffledCellIndices) {
                final int currentPopulationIndex = state.populationIndex(cellIndex);

                final double inploidyTerm = calculatePopulationAndGenomicAveragedPloidyExcludingPopulation(state, data, currentPopulationIndex, singleCellPopulationFraction);

                final double[] log10Probabilities = new double[state.numPopulations()];
                for (int populationIndex = 0; populationIndex < state.numPopulations(); populationIndex++) {
                    double logDensity = Math.log(state.populationFraction(populationIndex) + EPSILON);    //prior
                    for (int segmentIndex = 0; segmentIndex < state.numSegments(); segmentIndex++) {
                        final double segmentFractionalLength = state.calculateFractionalLength(data, segmentIndex);
                        final double invariantMAlleleCopyNumberTerm = calculatePopulationAveragedMAlleleCopyNumberExcludingPopulation(state, data, segmentIndex, currentPopulationIndex, singleCellPopulationFraction);
                        final double invariantNAlleleCopyNumberTerm = calculatePopulationAveragedNAlleleCopyNumberExcludingPopulation(state, data, segmentIndex, currentPopulationIndex, singleCellPopulationFraction);

                        final PloidyState ploidyState = state.ploidyState(populationIndex, segmentIndex);
                        logDensity += calculateSegmentLogLikelihoodFromInvariantTerms(data, inploidyTerm, invariantMAlleleCopyNumberTerm, invariantNAlleleCopyNumberTerm,
                                segmentIndex, singleCellPopulationFraction, segmentFractionalLength, ploidyState);
                    }
                    log10Probabilities[populationIndex] = MathUtils.logToLog10(logDensity);
                }
                final double[] probabilities = MathUtils.normalizeFromLog10(log10Probabilities);
                final int populationIndex = GATKProtectedMathUtils.randomSelect(populationIndices, i -> probabilities[i], rng);
                populationIndicators.set(cellIndex, populationIndex);

                //update the state as a side effect
                state.setPopulationIndexAndIncrementCounts(cellIndex, populationIndex);
            }
            logger.debug("Sampled population indicators: " + populationIndicators);
            return new TumorHeterogeneityState.PopulationIndicators(populationIndicators);
        }

        private static double calculatePopulationAndGenomicAveragedPloidyExcludingPopulation(final TumorHeterogeneityState state, final TumorHeterogeneityData data,
                                                                                             final int populationIndexToExclude, final double populationFractionToExclude) {
            final double excludedPloidy = populationFractionToExclude * IntStream.range(0, state.numSegments())
                    .mapToDouble(i -> state.calculateFractionalLength(data, i) * state.calculateCopyNumberFunction(i, populationIndexToExclude, PloidyState::total))
                    .sum();
            return state.calculatePopulationAndGenomicAveragedPloidy(data) - excludedPloidy;
        }

        private static double calculatePopulationAveragedMAlleleCopyNumberExcludingPopulation(final TumorHeterogeneityState state, final TumorHeterogeneityData data,
                                                                                              final int segmentIndex, final int populationIndexToExclude, final double populationFractionToExclude) {
            final double excludedMAlleleCopyNumber = state.calculateCopyNumberFunction(segmentIndex, populationIndexToExclude, PloidyState::m) * populationFractionToExclude;
            return state.calculatePopulationAveragedMAlleleCopyNumber(data, segmentIndex) - excludedMAlleleCopyNumber;
        }

        private static double calculatePopulationAveragedNAlleleCopyNumberExcludingPopulation(final TumorHeterogeneityState state, final TumorHeterogeneityData data,
                                                                                              final int segmentIndex, final int populationIndexToExclude, final double populationFractionToExclude) {
            final double excludedNAlleleCopyNumber = state.calculateCopyNumberFunction(segmentIndex, populationIndexToExclude, PloidyState::n) * populationFractionToExclude;
            return state.calculatePopulationAveragedNAlleleCopyNumber(data, segmentIndex) - excludedNAlleleCopyNumber;
        }
    }

    /**
     * Samples genomic profiles for a collection of variant populations.
     */
    protected static final class VariantProfileCollectionSampler implements ParameterSampler<TumorHeterogeneityState.VariantProfileCollection, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private final Random rnd = new Random(1562);
        private final int numVariantPopulations;
        private final List<VariantProfileSampler> variantProfileSamplers;

        protected VariantProfileCollectionSampler(final int numVariantPopulations, final PloidyStatePrior ploidyStatePrior) {
            this.numVariantPopulations = numVariantPopulations;
            variantProfileSamplers = IntStream.range(0, numVariantPopulations).boxed().map(n -> new VariantProfileSampler(n, ploidyStatePrior)).collect(Collectors.toList());
        }

        public TumorHeterogeneityState.VariantProfileCollection sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            logger.debug("Sampling variant profile collection.");
            return state.doMetropolisStep() ? new TumorHeterogeneityState.VariantProfileCollection(state.variantProfiles()) : sampleGibbs(rng, state, data);
        }

        private TumorHeterogeneityState.VariantProfileCollection sampleGibbs(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            final List<TumorHeterogeneityState.VariantProfile> shuffledVariantProfiles = new ArrayList<>(numVariantPopulations);
            final List<Integer> shuffledVariantPopulationIndices = IntStream.range(0, numVariantPopulations).boxed().collect(Collectors.toList());
            Collections.shuffle(shuffledVariantPopulationIndices, rnd);
            shuffledVariantPopulationIndices.stream().forEach(i -> shuffledVariantProfiles.add(variantProfileSamplers.get(i).sample(rng, state, data)));
            final List<TumorHeterogeneityState.VariantProfile> variantProfiles = IntStream.range(0, numVariantPopulations).boxed()
                    .map(i -> shuffledVariantProfiles.get(shuffledVariantPopulationIndices.indexOf(i))).collect(Collectors.toList());
            return new TumorHeterogeneityState.VariantProfileCollection(variantProfiles);
        }
    }

    /**
     * Samples a genomic profile for a single variant population.
     */
    private static final class VariantProfileSampler implements ParameterSampler<TumorHeterogeneityState.VariantProfile, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private int populationIndex;
        private final PloidyStateIndicatorsSampler ploidyStateIndicatorsSampler;

        private VariantProfileSampler(final int populationIndex, final PloidyStatePrior ploidyStatePrior) {
            this.populationIndex = populationIndex;
            ploidyStateIndicatorsSampler = new PloidyStateIndicatorsSampler(ploidyStatePrior);
        }

        @Override
        public TumorHeterogeneityState.VariantProfile sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            logger.debug("Sampling variant profile for population " + populationIndex);
            final TumorHeterogeneityState.VariantProfile.PloidyStateIndicators ploidyStateIndicators = ploidyStateIndicatorsSampler.sample(rng, state, data);
            return new TumorHeterogeneityState.VariantProfile(ploidyStateIndicators);
        }

        /**
         * Samples per-segment ploidy-state indicators for this variant population (these indicate the ploidy state if the population is variant in a segment).
         * Updates the ploidy-state indicators held by the {@link TumorHeterogeneityState}.
         */
        private final class PloidyStateIndicatorsSampler implements ParameterSampler<TumorHeterogeneityState.VariantProfile.PloidyStateIndicators, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
            private final Random rnd = new Random(4615);
            private final List<Integer> ploidyStateIndices;
            private final double[] ploidyStatePriorLog10Probabilities;

            private PloidyStateIndicatorsSampler(final PloidyStatePrior ploidyStatePrior) {
                final int numPloidyStates = ploidyStatePrior.numPloidyStates();
                ploidyStateIndices = Collections.unmodifiableList(IntStream.range(0, numPloidyStates).boxed().collect(Collectors.toList()));
                ploidyStatePriorLog10Probabilities = ploidyStatePrior.ploidyStates().stream()
                        .mapToDouble(vps -> MathUtils.logToLog10(ploidyStatePrior.logProbability(vps)))
                        .toArray();
            }

            @Override
            public TumorHeterogeneityState.VariantProfile.PloidyStateIndicators sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
                logger.debug("Sampling ploidy-state indicators for population " + populationIndex);
                final List<Integer> ploidyStateIndicators = new ArrayList<>(Collections.nCopies(state.numSegments(), 0));
                final List<Integer> shuffledSegmentIndices = IntStream.range(0, state.numSegments()).boxed().collect(Collectors.toList());
                Collections.shuffle(shuffledSegmentIndices, rnd);
                final List<PloidyState> ploidyStates = state.priors().ploidyStatePrior().ploidyStates();
                for (int segmentIndex : shuffledSegmentIndices) {
                    final double segmentFractionalLength = state.calculateFractionalLength(data, segmentIndex);
                    final double populationFraction = state.calculatePopulationFractionFromCounts(populationIndex);
                    final double inploidyTerm = calculatePopulationAndGenomicAveragedPloidyExcludingPopulationInSegment(state, data, segmentIndex, populationIndex);
                    final double invariantMAlleleCopyNumberTerm = calculatePopulationAveragedMAlleleCopyNumberExcludingPopulationInSegment(state, data, segmentIndex, populationIndex);
                    final double invariantNAlleleCopyNumberTerm = calculatePopulationAveragedNAlleleCopyNumberExcludingPopulationInSegment(state, data, segmentIndex, populationIndex);

                    //approximation: ignore coupling of copy-ratio posteriors in different segments due to ploidy term

                    //calculate unnormalized probabilities for all ploidy states
                    final int si = segmentIndex;
                    final double[] log10Probabilities = ploidyStateIndices.stream()
                                    .mapToDouble(i -> ploidyStatePriorLog10Probabilities[i] +
                                            MathUtils.logToLog10(calculateSegmentLogLikelihoodFromInvariantTerms(
                                                    data, inploidyTerm, invariantMAlleleCopyNumberTerm, invariantNAlleleCopyNumberTerm,
                                                    si, populationFraction, segmentFractionalLength, ploidyStates.get(i))))
                                    .toArray();
                    final double[] probabilities = MathUtils.normalizeFromLog10(log10Probabilities);
                    logger.debug("Population " + populationIndex + " segment " + segmentIndex + " ploidy-state probabilities: " + Doubles.asList(probabilities));
                    final int ploidyStateIndex = GATKProtectedMathUtils.randomSelect(ploidyStateIndices, i -> probabilities[i], rng);
                    ploidyStateIndicators.set(segmentIndex, ploidyStateIndex);

                    //update the state as a side effect
                    state.setPloidyStateIndex(populationIndex, segmentIndex, ploidyStateIndex);
                }
                logger.debug("Sampled ploidy states for population " + populationIndex + ": " + ploidyStateIndicators.stream().map(i -> ploidyStates.get(i).toString()).collect(Collectors.joining(" ")));
                return new TumorHeterogeneityState.VariantProfile.PloidyStateIndicators(ploidyStateIndicators);
            }
        }

        private static double calculatePopulationAndGenomicAveragedPloidyExcludingPopulationInSegment(final TumorHeterogeneityState state, final TumorHeterogeneityData data,
                                                                                                      final int segmentIndexToExclude, final int populationIndexToExclude) {
            final double excludedPloidy = state.calculatePopulationFractionFromCounts(populationIndexToExclude) * state.calculateFractionalLength(data, segmentIndexToExclude) * state.calculateCopyNumberFunction(segmentIndexToExclude, populationIndexToExclude, PloidyState::total);
            return state.calculatePopulationAndGenomicAveragedPloidy(data) - excludedPloidy;
        }

        private static double calculatePopulationAveragedMAlleleCopyNumberExcludingPopulationInSegment(final TumorHeterogeneityState state, final TumorHeterogeneityData data,
                                                                                                       final int segmentIndex, final int populationIndexToExclude) {
            return state.calculatePopulationAveragedCopyNumberFunction(data, segmentIndex, populationIndexToExclude, PloidyState::m, false);
        }

        private static double calculatePopulationAveragedNAlleleCopyNumberExcludingPopulationInSegment(final TumorHeterogeneityState state, final TumorHeterogeneityData data,
                                                                                                       final int segmentIndex, final int populationIndexToExclude) {
            return state.calculatePopulationAveragedCopyNumberFunction(data, segmentIndex, populationIndexToExclude, PloidyState::n, false);
        }
    }

    private static double calculateSegmentLogLikelihoodFromInvariantTerms(final TumorHeterogeneityData data,
                                                                          final double inploidyTerm,
                                                                          final double invariantMAlleleCopyNumberTerm,
                                                                          final double invariantNAlleleCopyNumberTerm,
                                                                          final int segmentIndex,
                                                                          final double populationFraction,
                                                                          final double segmentFractionalLength,
                                                                          final PloidyState ploidyState) {
        final double ploidy = inploidyTerm + populationFraction * segmentFractionalLength * ploidyState.total();
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
