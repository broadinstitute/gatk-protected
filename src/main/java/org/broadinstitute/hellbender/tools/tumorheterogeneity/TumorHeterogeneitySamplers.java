package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Sets;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Gamma;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
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
    static final double EPSILON = 1E-10;
    private static final double COPY_RATIO_NOISE_FLOOR_MAX = 1E-2;

    static final Logger logger = LogManager.getLogger(TumorHeterogeneitySamplers.class);

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
            final int numPopulations = state.populationMixture().numPopulations();
            final double concentrationPriorAlpha = state.priors().concentrationPriorHyperparameterValues().getAlpha();
            final double concentrationPriorBeta = state.priors().concentrationPriorHyperparameterValues().getBeta();
            final Function<Double, Double> logConditionalPDF = newConcentration -> {
                final double populationFractionsTerm = IntStream.range(0, numPopulations)
                        .mapToDouble(i -> (newConcentration - 1) * Math.log(state.populationMixture().populationFraction(i) + EPSILON)).sum();
                return (concentrationPriorAlpha - 1.) * Math.log(newConcentration) - concentrationPriorBeta * newConcentration +
                        Gamma.logGamma(newConcentration * numPopulations) - numPopulations * Gamma.logGamma(newConcentration) + populationFractionsTerm;
            };
            final double concentration = new SliceSampler(rng, logConditionalPDF, concentrationMin, concentrationMax, concentrationSliceSamplingWidth).sample(state.concentration());
            logger.debug("Sampled concentration: " + concentration);
            return concentration;
        }
    }

    static final class CopyRatioNoiseFloorSampler implements ParameterSampler<Double, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private static final double MIN = 0.;
        private static final double MAX = COPY_RATIO_NOISE_FLOOR_MAX;
        private final double copyRatioNoiseFloorSliceSamplingWidth;

        CopyRatioNoiseFloorSampler(final double copyRatioNoiseFloorSliceSamplingWidth) {
            this.copyRatioNoiseFloorSliceSamplingWidth = copyRatioNoiseFloorSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            final Function<Double, Double> logConditionalPDF = newCopyRatioNoiseFloor -> {
                final TumorHeterogeneityState newState = new TumorHeterogeneityState(
                        state.concentration(),
                        newCopyRatioNoiseFloor,
                        state.copyRatioNoiseFactor(),
                        state.minorAlleleFractionNoiseFactor(),
                        state.populationMixture(),
                        state.priors());
                return TumorHeterogeneityUtils.calculateLogPosterior(newState, data);
            };
            final double copyRatioNoiseFloor = new SliceSampler(rng, logConditionalPDF, MIN, MAX, copyRatioNoiseFloorSliceSamplingWidth).sample(state.copyRatioNoiseFloor());
            logger.debug("Sampled CR noise floor: " + copyRatioNoiseFloor);
            return copyRatioNoiseFloor;
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
                        state.concentration(),
                        state.copyRatioNoiseFloor(),
                        newCopyRatioNoiseFactor,
                        state.minorAlleleFractionNoiseFactor(),
                        state.populationMixture(),
                        state.priors());
                return TumorHeterogeneityUtils.calculateLogPosterior(newState, data);
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
                        state.concentration(),
                        state.copyRatioNoiseFloor(),
                        state.copyRatioNoiseFactor(),
                        newMinorAlleleFractionNoiseFactor,
                        state.populationMixture(),
                        state.priors());
                return TumorHeterogeneityUtils.calculateLogPosterior(newState, data);
            };
            final double minorAlleleFractionNoiseFactor = new SliceSampler(rng, logConditionalPDF, MIN, MAX, minorAlleleFractionNoiseFactorSliceSamplingWidth).sample(state.minorAlleleFractionNoiseFactor());
            logger.debug("Sampled MAF noise factor: " + minorAlleleFractionNoiseFactor);
            return minorAlleleFractionNoiseFactor;
        }
    }

    static final class PopulationMixtureSampler implements ParameterSampler<PopulationMixture, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private static double maxLogPosterior = -1E100;
        private static final double scaleParameter = 2.;
        private static final int numWalkers = 10;

        private int numSamples = 0;
        private int numAccepted = 0;

        private final int maxTotalCopyNumber;
        private final List<List<Integer>> totalCopyNumberProductStates;
        private final Map<Integer, Set<PloidyState>> ploidyStateSetsMap = new HashMap<>();

        private final List<List<Double>> walkerPositions;

        PopulationMixtureSampler(final RandomGenerator rng, final int numPopulations, final List<PloidyState> ploidyStates, final PloidyState normalPloidyState) {
            final Set<Integer> totalCopyNumberStates = ploidyStates.stream().map(PloidyState::total).collect(Collectors.toSet());
            maxTotalCopyNumber = Collections.max(totalCopyNumberStates);
            totalCopyNumberProductStates = new ArrayList<>(Sets.cartesianProduct(Collections.nCopies(numPopulations, totalCopyNumberStates)));
            for (final int totalCopyNumber : totalCopyNumberStates) {
                final Set<PloidyState> ploidyStateSet = ploidyStates.stream().filter(ps -> ps.total() == totalCopyNumber).collect(Collectors.toSet());
                ploidyStateSetsMap.put(totalCopyNumber, ploidyStateSet);
            }

            walkerPositions = new ArrayList<>(numWalkers);
            for (int walkerIndex = 0; walkerIndex < numWalkers; walkerIndex++) {
                final List<Double> walkerPosition = IntStream.range(0, numPopulations - 1).boxed()
                        .map(i -> rng.nextGaussian()).collect(Collectors.toList());
                walkerPosition.add(Math.max(EPSILON, Math.min(normalPloidyState.total() + rng.nextGaussian(), maxTotalCopyNumber)));
                walkerPositions.add(walkerPosition);
            }
        }

        @Override
        public PopulationMixture sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            final int numPopulations = state.populationMixture().numPopulations();
            final PloidyState normalPloidyState = state.priors().normalPloidyState();

            final int currentWalkerIndex = numSamples % numWalkers;
            final int selectedWalkerIndex = IntStream.range(0, numWalkers).boxed().filter(i -> i != currentWalkerIndex).collect(Collectors.toList()).get(rng.nextInt(numWalkers - 1));
            final List<Double> currentWalkerPosition = walkerPositions.get(currentWalkerIndex);
            final List<Double> selectedWalkerPosition = walkerPositions.get(selectedWalkerIndex);

            final double z = Math.pow((scaleParameter - 1.) * rng.nextDouble() + 1, 2.) / scaleParameter;
            final List<Double> proposedWalkerPosition = IntStream.range(0, numPopulations).boxed()
                    .map(i -> selectedWalkerPosition.get(i) + z * (currentWalkerPosition.get(i) - selectedWalkerPosition.get(i)))
                    .collect(Collectors.toList());
            proposedWalkerPosition.set(numPopulations - 1, Math.max(EPSILON, Math.min(normalPloidyState.total(), proposedWalkerPosition.get(numPopulations - 1))));

            final List<Double> currentTransformedPopulationFractions = currentWalkerPosition.subList(0, numPopulations - 1);
            final PopulationMixture.PopulationFractions currentPopulationFractions =
                    TumorHeterogeneityUtils.calculatePopulationFractionsFromTransformedPopulationFractions(currentTransformedPopulationFractions);
            final double currentInitialPloidy = currentWalkerPosition.get(numPopulations - 1);
            final PopulationMixture.VariantProfileCollection currentVariantProfileCollection =
                    TumorHeterogeneityUtils.proposeVariantProfileCollection(rng, state, data,
                            currentPopulationFractions, currentInitialPloidy, totalCopyNumberProductStates, ploidyStateSetsMap);
            final PopulationMixture currentPopulationMixture = new PopulationMixture(
                    currentPopulationFractions, currentVariantProfileCollection, normalPloidyState);
            final double currentResultingPloidy = currentPopulationMixture.ploidy(data);
            logger.debug("Current population fractions: " + currentPopulationFractions);
            logger.debug("Current ploidy: " + currentResultingPloidy);

            final TumorHeterogeneityState currentState = new TumorHeterogeneityState(
                    state.concentration(),
                    state.copyRatioNoiseFloor(),
                    state.copyRatioNoiseFactor(),
                    state.minorAlleleFractionNoiseFactor(),
                    currentPopulationMixture,
                    state.priors());

            final List<Double> proposedTransformedPopulationFractions = proposedWalkerPosition.subList(0, numPopulations - 1);
            final PopulationMixture.PopulationFractions proposedPopulationFractions =
                    TumorHeterogeneityUtils.calculatePopulationFractionsFromTransformedPopulationFractions(proposedTransformedPopulationFractions);
            final double proposedInitialPloidy = proposedWalkerPosition.get(numPopulations - 1);
            final PopulationMixture.VariantProfileCollection proposedVariantProfileCollection =
                    TumorHeterogeneityUtils.proposeVariantProfileCollection(rng, state, data,
                            proposedPopulationFractions, proposedInitialPloidy, totalCopyNumberProductStates, ploidyStateSetsMap);
            final PopulationMixture proposedPopulationMixture = new PopulationMixture(
                    proposedPopulationFractions, proposedVariantProfileCollection, normalPloidyState);
            final double proposedResultingPloidy = proposedPopulationMixture.ploidy(data);
            logger.debug("Proposed population fractions: " + proposedPopulationFractions);
            logger.debug("Proposed initial ploidy: " + proposedInitialPloidy);
            logger.debug("Proposed resulting ploidy: " + proposedResultingPloidy);

            final TumorHeterogeneityState proposedState = new TumorHeterogeneityState(
                    state.concentration(),
                    state.copyRatioNoiseFloor(),
                    state.copyRatioNoiseFactor(),
                    state.minorAlleleFractionNoiseFactor(),
                    proposedPopulationMixture,
                    state.priors());

            final double proposedLogPosterior = TumorHeterogeneityUtils.calculateLogPriorVariantProfiles(proposedVariantProfileCollection, state.priors().ploidyStatePrior())
                    + TumorHeterogeneityUtils.calculateLogLikelihoodSegments(proposedState, data)
                    + TumorHeterogeneityUtils.calculateLogJacobianFactor(proposedPopulationFractions);
            final double currentLogPosterior = TumorHeterogeneityUtils.calculateLogPriorVariantProfiles(currentVariantProfileCollection, state.priors().ploidyStatePrior())
                    + TumorHeterogeneityUtils.calculateLogLikelihoodSegments(currentState, data)
                    + TumorHeterogeneityUtils.calculateLogJacobianFactor(currentPopulationFractions);
            final double acceptanceProbability = Math.min(1., Math.exp((numPopulations - 1.) * Math.log(z) + proposedLogPosterior - currentLogPosterior));
            logger.debug("Log posterior of current state: " + currentLogPosterior);
            logger.debug("Log posterior of proposed state: " + proposedLogPosterior);
            numSamples += 1;
            if (proposedLogPosterior > maxLogPosterior) {
                maxLogPosterior = currentLogPosterior;
                logger.debug("New maximum: " + maxLogPosterior);
            }
            if (rng.nextDouble() < acceptanceProbability) {
                numAccepted += 1;
                logger.debug("Proposed state accepted.");
                walkerPositions.set(currentWalkerIndex, proposedWalkerPosition);
                return proposedPopulationMixture;
            }

            logger.debug("Acceptance rate: " + (double) numAccepted / numSamples);
            return currentPopulationMixture;
        }

//        @Override
//        public PopulationMixture sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
//            final PloidyState normalPloidyState = state.priors().normalPloidyState();
//            final PopulationMixture currentPopulationMixture = state.populationMixture();
//            final PopulationMixture.PopulationFractions currentPopulationFractions = currentPopulationMixture.populationFractions();
//            final List<Double> currentTransformedPopulationFractions =
//                    TumorHeterogeneityUtils.calculateTransformedPopulationFractionsFromPopulationFractions(currentPopulationFractions);
//            final double currentPloidy = currentPopulationMixture.ploidy(data);
//            logger.debug("Current population fractions: " + currentPopulationFractions);
//            logger.debug("Current ploidy: " + currentPloidy);
//
//            final List<Double> proposedTransformedPopulationFractions = currentTransformedPopulationFractions.stream()
//                    .map(x -> x + new NormalDistribution(rng, 0., transformedPopulationFractionProposalWidth).sample())
//                    .collect(Collectors.toList());
//            final PopulationMixture.PopulationFractions proposedPopulationFractions =
//                    TumorHeterogeneityUtils.calculatePopulationFractionsFromTransformedPopulationFractions(proposedTransformedPopulationFractions);
//
//            final double normalPloidy = state.priors().normalPloidyState().total();
//            final double currentTumorPloidy = currentPloidy - currentPopulationFractions.normalFraction() * normalPloidy;
//            final double proposedTumorPloidy = currentTumorPloidy + new NormalDistribution(rng, 0., ploidyProposalWidth).sample();
//            final double proposedInitialPloidy = Math.min(maxTotalCopyNumber, Math.max(EPSILON,
//                    proposedPopulationFractions.normalFraction() * normalPloidy + proposedPopulationFractions.tumorFraction() * proposedTumorPloidy));
//            logger.debug("Proposed initial ploidy: " + proposedInitialPloidy);
//
//            final PopulationMixture.VariantProfileCollection proposedVariantProfileCollection =
//                    TumorHeterogeneityUtils.proposeVariantProfileCollection(rng, state, data,
//                            proposedPopulationFractions, proposedInitialPloidy, totalCopyNumberProductStates, ploidyStateSetsMap);
//            final PopulationMixture proposedPopulationMixture = new PopulationMixture(
//                    proposedPopulationFractions, proposedVariantProfileCollection, normalPloidyState);
//            final double proposedResultingPloidy = proposedPopulationMixture.ploidy(data);
//            logger.debug("Proposed population fractions: " + proposedPopulationFractions);
//            logger.debug("Proposed resulting ploidy: " + proposedResultingPloidy);
//
//            final TumorHeterogeneityState proposedState = new TumorHeterogeneityState(
//                    state.concentration(),
//                    state.copyRatioNoiseFloor(),
//                    state.copyRatioNoiseFactor(),
//                    state.minorAlleleFractionNoiseFactor(),
//                    proposedPopulationMixture,
//                    state.priors());
//
//            final double proposedLogPosterior = TumorHeterogeneityUtils.calculateLogPosterior(proposedState, data)
//                    + TumorHeterogeneityUtils.calculateLogJacobianFactor(proposedPopulationFractions);
//            final double currentLogPosterior = TumorHeterogeneityUtils.calculateLogPosterior(state, data)
//                    + TumorHeterogeneityUtils.calculateLogJacobianFactor(currentPopulationFractions);
//            final double acceptanceProbability = Math.min(1., Math.exp(proposedLogPosterior - currentLogPosterior));
//            logger.debug("Log posterior of current state: " + currentLogPosterior);
//            logger.debug("Log posterior of proposed state: " + proposedLogPosterior);
//            numSamples += 1;
//            if (proposedLogPosterior > maxLogPosterior) {
//                maxLogPosterior = currentLogPosterior;
//                logger.debug("New maximum: " + maxLogPosterior);
//            }
//            if (rng.nextDouble() < acceptanceProbability) {
//                numAccepted += 1;
//                logger.debug("Proposed state accepted.");
//                return proposedPopulationMixture;
//            }
//
//            logger.debug("Acceptance rate: " + (double) numAccepted / numSamples);
//            return new PopulationMixture(
//                    currentPopulationMixture.populationFractions(),
//                    currentPopulationMixture.variantProfileCollection(),
//                    state.priors().normalPloidyState()
//            );
//        }
    }
}
