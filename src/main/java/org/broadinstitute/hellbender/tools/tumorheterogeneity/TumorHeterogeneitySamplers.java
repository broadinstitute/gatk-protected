package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.Dirichlet;
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

    private static final Logger logger = LogManager.getLogger(TumorHeterogeneitySamplers.class);

    private TumorHeterogeneitySamplers() {}

    static final class ConcentrationSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private final double concentrationMin;
        private final double concentrationMax;
        private final double concentrationSliceSamplingWidth;

        ConcentrationSampler(final double concentrationMin, final double concentrationMax, final double concentrationSliceSamplingWidth) {
            this.concentrationMin = concentrationMin;
            this.concentrationMax = concentrationMax;
            this.concentrationSliceSamplingWidth = concentrationSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            logger.info("Ploidy of current state: " + state.ploidy(data));
            final int numPopulations = state.numPopulations();
            final Function<Double, Double> logConditionalPDF = newConcentration -> {
                final double populationFractionsTerm = IntStream.range(0, numPopulations)
                        .mapToDouble(i -> (newConcentration - 1) * Math.log(state.populationFraction(i) + EPSILON)).sum();
                return (state.priors().concentrationPriorAlpha() - 1.) * Math.log(newConcentration) - state.priors().concentrationPriorBeta() * newConcentration +
                        Gamma.logGamma(newConcentration * numPopulations) - numPopulations * Gamma.logGamma(newConcentration) + populationFractionsTerm;
            };
            final double concentration = new SliceSampler(rng, logConditionalPDF, concentrationMin, concentrationMax, concentrationSliceSamplingWidth).sample(state.concentration());
            logger.debug("Sampled concentration: " + concentration);
            return concentration;
        }
    }

    static final class CopyRatioNoiseFactorSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private static final double MIN = 1.;
        private static final double MAX = Double.POSITIVE_INFINITY;
        private final double copyRatioNoiseFactorSliceSamplingWidth;

        CopyRatioNoiseFactorSampler(final double copyRatioNoiseFactorSliceSamplingWidth) {
            this.copyRatioNoiseFactorSliceSamplingWidth = copyRatioNoiseFactorSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            final Function<Double, Double> logConditionalPDF = newCopyRatioNoiseFactor -> {
                final TumorHeterogeneityState newState = new TumorHeterogeneityState(
                        state.concentration(), newCopyRatioNoiseFactor, state.minorAlleleFractionNoiseFactor(), state.populationFractions(), state.variantProfiles(), state.priors());
                return calculateLogPosterior(newState, data);
            };
            final double copyRatioNoiseFactor = new SliceSampler(rng, logConditionalPDF, MIN, MAX, copyRatioNoiseFactorSliceSamplingWidth).sample(state.copyRatioNoiseFactor());
            logger.debug("Sampled CR noise factor: " + copyRatioNoiseFactor);
            return copyRatioNoiseFactor;
        }
    }

    static final class MinorAlleleFractionNoiseFactorSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private static final double MIN = 1.;
        private static final double MAX = Double.POSITIVE_INFINITY;
        private final double minorAlleleFractionNoiseFactorSliceSamplingWidth;

        MinorAlleleFractionNoiseFactorSampler(final double minorAlleleFractionNoiseFactorSliceSamplingWidth) {
            this.minorAlleleFractionNoiseFactorSliceSamplingWidth = minorAlleleFractionNoiseFactorSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            final Function<Double, Double> logConditionalPDF = newMinorAlleleFractionNoiseFactor -> {
                final TumorHeterogeneityState newState = new TumorHeterogeneityState(
                        state.concentration(), state.copyRatioNoiseFactor(), newMinorAlleleFractionNoiseFactor, state.populationFractions(), state.variantProfiles(), state.priors());
                return calculateLogPosterior(newState, data);
            };
            final double minorAlleleFractionNoiseFactor = new SliceSampler(rng, logConditionalPDF, MIN, MAX, minorAlleleFractionNoiseFactorSliceSamplingWidth).sample(state.minorAlleleFractionNoiseFactor());
            logger.debug("Sampled MAF noise factor: " + minorAlleleFractionNoiseFactor);
            return minorAlleleFractionNoiseFactor;
        }
    }

    static final class PopulationFractionsSampler implements ParameterSampler<TumorHeterogeneityState.PopulationFractions, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private double max = -1E10;
        private static int numSamples = 0;
        private static int numAccepted = 0;
        PopulationFractionsSampler() {}

        @Override
        public TumorHeterogeneityState.PopulationFractions sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            final TumorHeterogeneityState.PopulationFractions populationFractions = sampleMetropolis(rng, state, data);
            logger.debug("Sampled population fractions: " + populationFractions);
            return populationFractions;
        }

        private TumorHeterogeneityState.PopulationFractions sampleMetropolis(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            final TumorHeterogeneityState proposedState = TumorHeterogeneityStateInitializationUtils.proposeState(rng, state, data);
            final double proposedLogPosterior = calculateLogPosterior(proposedState, data);
            final double currentLogPosterior = calculateLogPosterior(state, data);
            final double logProposalRatio = calculateLogProposalRatio(state, proposedState);
            final double acceptanceProbability = Math.min(1., Math.exp(proposedLogPosterior - currentLogPosterior + logProposalRatio));
            logger.debug("Log posterior of current state: " + currentLogPosterior);
            logger.debug("Log posterior of proposed state: " + proposedLogPosterior);
            numSamples += 1;
            if (rng.nextDouble() < acceptanceProbability) {
                numAccepted += 1;
                logger.info("Proposed state accepted.");
                state.set(proposedState);
            }
            if (currentLogPosterior > max) {
                max = currentLogPosterior;
                logger.info("New maximum: " + max);
            }
            logger.info("Acceptance rate: " + (double) numAccepted / numSamples);
            return new TumorHeterogeneityState.PopulationFractions(state.populationFractions());
        }

        private static double calculateLogProposalRatio(final TumorHeterogeneityState currentState,
                                                        final TumorHeterogeneityState proposedState) {
            final double currentLogProposalProbability = calculateLogProposalProbability(currentState, proposedState);
            final double proposedLogProposalProbability = calculateLogProposalProbability(proposedState, currentState);
            logger.debug("Proposal log probability of current state: " + currentLogProposalProbability);
            logger.debug("Proposal log probability of proposed state: " + proposedLogProposalProbability);
            return currentLogProposalProbability - proposedLogProposalProbability;
        }

        private static double calculateLogProposalProbability(final TumorHeterogeneityState newState,
                                                              final TumorHeterogeneityState conditionalState) {
            final int numPopulations = conditionalState.numPopulations();
            final double proposalWidthFactor = conditionalState.priors().proposalWidthFactor();
            final double concentration = conditionalState.concentration();
            final Dirichlet prior = Dirichlet.symmetricDirichlet(numPopulations, concentration * numPopulations);   //concentration convention differs from that used in Dirichlet class
            final double[] proposedPopulationFractions = Doubles.toArray(newState.populationFractions());
            final double[] currentEffectiveCounts = conditionalState.populationFractions().stream().mapToDouble(f -> proposalWidthFactor * f).toArray();
            final Dirichlet proposal = new Dirichlet(prior, currentEffectiveCounts);
            return proposal.logDensity(proposedPopulationFractions);
        }
    }

    /**
     * Samples genomic profiles for a collection of variant populations.
     */
    static final class VariantProfileCollectionSampler implements ParameterSampler<TumorHeterogeneityState.VariantProfileCollection, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private final int numVariantPopulations;
        private final List<Integer> ploidyStateIndices;
        private final double[] ploidyStatePriorLog10Probabilities;

        VariantProfileCollectionSampler(final int numVariantPopulations, final PloidyStatePrior ploidyStatePrior) {
            this.numVariantPopulations = numVariantPopulations;
            ploidyStateIndices = Collections.unmodifiableList(IntStream.range(0, ploidyStatePrior.numPloidyStates()).boxed().collect(Collectors.toList()));
            ploidyStatePriorLog10Probabilities = ploidyStatePrior.ploidyStates().stream()
                    .mapToDouble(vps -> MathUtils.logToLog10(ploidyStatePrior.logProbability(vps)))
                    .toArray();
        }

        public TumorHeterogeneityState.VariantProfileCollection sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            return new TumorHeterogeneityState.VariantProfileCollection(state.variantProfiles());
        }

        TumorHeterogeneityState.VariantProfileCollection sampleGibbs(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            final List<TumorHeterogeneityState.VariantProfile> variantProfiles = new ArrayList<>(Collections.nCopies(numVariantPopulations,
                    new TumorHeterogeneityState.VariantProfile(new TumorHeterogeneityState.VariantProfile(Collections.nCopies(state.numSegments(), 0)))));
            final List<PloidyState> ploidyStates = state.priors().ploidyStatePrior().ploidyStates();

            final List<Integer> shuffledSegmentIndices = IntStream.range(0, state.numSegments()).boxed().collect(Collectors.toList());
            final List<Integer> shuffledPopulationIndices = IntStream.range(0, state.numPopulations() - 1).boxed().collect(Collectors.toList());
            final List<Pair<Integer, Integer>> shuffledPopulationAndSegmentIndices = new ArrayList<>();
//            Collections.shuffle(shuffledPopulationIndices, new Random(rng.nextLong()));
            for (final int populationIndex : shuffledPopulationIndices) {
//                Collections.shuffle(shuffledSegmentIndices, new Random(rng.nextLong()));
                for (final int segmentIndex : shuffledSegmentIndices) {
                    shuffledPopulationAndSegmentIndices.add(new Pair<>(populationIndex, segmentIndex));
                }
            }
            Collections.shuffle(shuffledPopulationAndSegmentIndices, new Random(rng.nextLong()));

            double ploidy = state.ploidy(data);// + new NormalDistribution(rng, 0, 0.5).sample();
            for (final Pair<Integer, Integer> populationAndSegmentIndices : shuffledPopulationAndSegmentIndices) {
                final int populationIndex = populationAndSegmentIndices.getFirst();
                final int segmentIndex = populationAndSegmentIndices.getSecond();
                final double segmentFractionalLength = data.fractionalLength(segmentIndex);
                final double populationFraction = state.populationFraction(populationIndex);

                ploidy -= populationFraction * segmentFractionalLength * state.calculateCopyNumberFunction(segmentIndex, populationIndex, PloidyState::total);
                final double invariantPloidyTerm = ploidy;
                final double invariantMAlleleCopyNumberTerm = state.calculatePopulationAveragedCopyNumberFunctionExcludingPopulation(populationIndex, segmentIndex, PloidyState::m);
                final double invariantNAlleleCopyNumberTerm = state.calculatePopulationAveragedCopyNumberFunctionExcludingPopulation(populationIndex, segmentIndex, PloidyState::n);

                //approximation: ignore coupling of copy-ratio posteriors in different segments due to ploidy term

                //calculate unnormalized probabilities for all ploidy states
                final double[] log10Probabilities = ploidyStateIndices.stream()
                        .mapToDouble(i -> ploidyStatePriorLog10Probabilities[i] +
                                MathUtils.logToLog10(calculateSegmentLogLikelihoodFromInvariantTerms(
                                        data, invariantPloidyTerm, invariantMAlleleCopyNumberTerm, invariantNAlleleCopyNumberTerm,
                                        segmentIndex, populationFraction, segmentFractionalLength, ploidyStates.get(i),
                                        state.copyRatioNoiseFactor(), state.minorAlleleFractionNoiseFactor())))
                        .toArray();
                final double[] probabilities = MathUtils.normalizeFromLog10(log10Probabilities);
                final int ploidyStateIndex = GATKProtectedMathUtils.randomSelect(ploidyStateIndices, i -> probabilities[i], rng);
                variantProfiles.get(populationIndex).set(segmentIndex, ploidyStateIndex);

                //update the current state as a side effect
                state.setPloidyStateIndex(populationIndex, segmentIndex, ploidyStateIndex);
                ploidy += populationFraction * segmentFractionalLength * state.calculateCopyNumberFunction(segmentIndex, populationIndex, PloidyState::total);
            }
            return new TumorHeterogeneityState.VariantProfileCollection(variantProfiles);
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
                        + (concentrationPriorAlpha - 1.) * Math.log(concentration + EPSILON)
                        - concentrationPriorBeta * concentration
                        - Gamma.logGamma(concentrationPriorAlpha);

        //copy-ratio noise-factor prior
        final double copyRatioNoiseFactorPriorAlpha = state.priors().copyRatioNoiseFactorPriorAlpha();
        final double copyRatioNoiseFactorPriorBeta = state.priors().copyRatioNoiseFactorPriorBeta();
        final double copyRatioNoiseFactor = state.copyRatioNoiseFactor();
        final double logPriorCopyRatioNoiseFactor =
                concentrationPriorAlpha * Math.log(copyRatioNoiseFactorPriorBeta + EPSILON)
                        + (copyRatioNoiseFactorPriorAlpha - 1.) * Math.log(copyRatioNoiseFactor - 1. + EPSILON)
                        - copyRatioNoiseFactorPriorBeta * (copyRatioNoiseFactor - 1.)
                        - Gamma.logGamma(copyRatioNoiseFactorPriorAlpha);

        //minor-allele-fraction noise-factor prior
        final double minorAlleleFractionNoiseFactorPriorAlpha = state.priors().minorAlleleFractionNoiseFactorPriorAlpha();
        final double minorAlleleFractionNoiseFactorPriorBeta = state.priors().minorAlleleFractionNoiseFactorPriorBeta();
        final double minorAlleleFractionNoiseFactor = state.minorAlleleFractionNoiseFactor();
        final double logPriorMinorAlleleFractionNoiseFactor =
                concentrationPriorAlpha * Math.log(minorAlleleFractionNoiseFactorPriorBeta + EPSILON)
                        + (minorAlleleFractionNoiseFactorPriorAlpha - 1.) * Math.log(minorAlleleFractionNoiseFactor - 1. + EPSILON)
                        - minorAlleleFractionNoiseFactorPriorBeta * (minorAlleleFractionNoiseFactor - 1.)
                        - Gamma.logGamma(minorAlleleFractionNoiseFactorPriorAlpha);

        //population-fractions prior
        final double logPriorPopulationFractionsSum = IntStream.range(0, numPopulations)
                .mapToDouble(i -> (concentration - 1.) * Math.log(state.populationFraction(i) + EPSILON))
                .sum();
        final double logPriorPopulationFractions =
                Gamma.logGamma(concentration * numPopulations)
                        - numPopulations * Gamma.logGamma(concentration)
                        + logPriorPopulationFractionsSum;

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
        final double ploidy = state.ploidy(data);
        for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
            final double totalCopyNumber = state.calculatePopulationAveragedCopyNumberFunction(segmentIndex, PloidyState::total);
            final double mAlleleCopyNumber = state.calculatePopulationAveragedCopyNumberFunction(segmentIndex, PloidyState::m);
            final double nAlleleCopyNumber = state.calculatePopulationAveragedCopyNumberFunction(segmentIndex, PloidyState::n);
            final double copyRatio = totalCopyNumber / (ploidy + EPSILON);
            final double minorAlleleFraction = calculateMinorAlleleFraction(mAlleleCopyNumber, nAlleleCopyNumber);
            logLikelihoodSegments += data.logDensity(segmentIndex, copyRatio, minorAlleleFraction, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor);
        }
        return logPriorConcentration + logPriorCopyRatioNoiseFactor + logPriorMinorAlleleFractionNoiseFactor +
                logPriorPopulationFractions + logPriorVariantProfiles + logLikelihoodSegments;
    }

    private static double calculateSegmentLogLikelihoodFromInvariantTerms(final TumorHeterogeneityData data,
                                                                          final double invariantPloidyTerm,
                                                                          final double invariantMAlleleCopyNumberTerm,
                                                                          final double invariantNAlleleCopyNumberTerm,
                                                                          final int segmentIndex,
                                                                          final double populationFraction,
                                                                          final double segmentFractionalLength,
                                                                          final PloidyState ploidyState,
                                                                          final double copyRatioNoiseFactor,
                                                                          final double minorAlleleFractionNoiseFactor) {
        final double ploidy = invariantPloidyTerm + populationFraction * segmentFractionalLength * ploidyState.total();
        final double copyRatio = (invariantMAlleleCopyNumberTerm + invariantNAlleleCopyNumberTerm + populationFraction * ploidyState.total()) / (ploidy + EPSILON);
        final double minorAlleleFraction = calculateMinorAlleleFraction(
                invariantMAlleleCopyNumberTerm + populationFraction * ploidyState.m(),
                invariantNAlleleCopyNumberTerm + populationFraction * ploidyState.n());
        return data.logDensity(segmentIndex, copyRatio, minorAlleleFraction, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor);
    }

    private static double calculateMinorAlleleFraction(final double m, final double n) {
        return Math.min(m, n) / (m + n + EPSILON);
    }
}
