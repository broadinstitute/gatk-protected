package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Sets;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Gamma;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
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
        private static final double MIN = 0.;
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
        private static final int maxNumIterationsPloidyStep = 25;
        private static final double transformedPopulationFractionProposalWidth = 0.1;
        private static final double ploidyProposalWidth = 0.05;

        private static double maxLogPosterior = -1E100;
        private static int numSamples = 0;
        private static int numAccepted = 0;

        private final int maxTotalCopyNumber;
        private final List<List<Integer>> totalCopyNumberProductStates;
        private final Map<Integer, Set<PloidyState>> ploidyStateSetsMap = new HashMap<>();

        PopulationMixtureSampler(final int numPopulations, final int maxTotalCopyNumber) {
            this.maxTotalCopyNumber = maxTotalCopyNumber;
            final Set<Integer> copyNumberStates = IntStream.range(0, maxTotalCopyNumber + 1).boxed().collect(Collectors.toSet());
            totalCopyNumberProductStates = new ArrayList<>(Sets.cartesianProduct(Collections.nCopies(numPopulations, copyNumberStates)));
            for (int totalCopyNumber = 0; totalCopyNumber <= maxTotalCopyNumber; totalCopyNumber++) {
                final Set<PloidyState> ploidyStateSet = new HashSet<>();
                for (int m = 0; m <= totalCopyNumber / 2; m++) {
                    ploidyStateSet.add(new PloidyState(m, totalCopyNumber - m));
                }
                ploidyStateSetsMap.put(totalCopyNumber, ploidyStateSet);
            }
        }

        @Override
        public PopulationMixture sample(final RandomGenerator rng, final TumorHeterogeneityState state, final TumorHeterogeneityData data) {
            final PloidyState normalPloidyState = state.priors().normalPloidyState();
            final PopulationMixture currentPopulationMixture = state.populationMixture();
            final PopulationMixture.PopulationFractions currentPopulationFractions = currentPopulationMixture.populationFractions();
            final List<Double> currentTransformedPopulationFractions =
                    TumorHeterogeneityUtils.calculateTransformedPopulationFractionsFromPopulationFractions(currentPopulationFractions);

            final List<Double> proposedTransformedPopulationFractions = currentTransformedPopulationFractions.stream()
                    .map(x -> TumorHeterogeneityUtils.proposeTransformedPopulationFraction(rng, x, transformedPopulationFractionProposalWidth))
                    .collect(Collectors.toList());
            final PopulationMixture.PopulationFractions proposedPopulationFractions =
                    TumorHeterogeneityUtils.calculatePopulationFractionsFromTransformedPopulationFractions(proposedTransformedPopulationFractions);

            final PopulationMixture.VariantProfileCollection proposedVariantProfileCollection =
                    TumorHeterogeneityUtils.proposeVariantProfileCollection(rng, state, data,
                            proposedPopulationFractions, transformedPopulationFractionProposalWidth,
                            ploidyProposalWidth, maxTotalCopyNumber, maxNumIterationsPloidyStep, totalCopyNumberProductStates, ploidyStateSetsMap);
            final PopulationMixture proposedPopulationMixture = new PopulationMixture(
                    proposedPopulationFractions, proposedVariantProfileCollection, normalPloidyState);
            logger.info("Proposed population fractions: " + proposedPopulationFractions);
            logger.info("Proposed ploidy: " + proposedPopulationMixture.ploidy(data));

            final TumorHeterogeneityState proposedState = new TumorHeterogeneityState(
                    state.concentration(),
                    state.copyRatioNoiseFloor(),
                    state.copyRatioNoiseFactor(),
                    state.minorAlleleFractionNoiseFactor(),
                    proposedPopulationMixture,
                    state.priors());

            final double proposedLogPosterior = TumorHeterogeneityUtils.calculateLogPosterior(proposedState, data) + TumorHeterogeneityUtils.calculateLogJacobianFactor(proposedPopulationFractions);
            final double currentLogPosterior = TumorHeterogeneityUtils.calculateLogPosterior(state, data) + TumorHeterogeneityUtils.calculateLogJacobianFactor(currentPopulationFractions);
            final double acceptanceProbability = Math.min(1., Math.exp(proposedLogPosterior - currentLogPosterior));
            logger.debug("Log posterior of current state: " + currentLogPosterior);
            logger.debug("Log posterior of proposed state: " + proposedLogPosterior);
            numSamples += 1;
            if (proposedLogPosterior > maxLogPosterior) {
                maxLogPosterior = currentLogPosterior;
                logger.info("New maximum: " + maxLogPosterior);
            }
            if (rng.nextDouble() < acceptanceProbability) {
                numAccepted += 1;
                logger.info("Proposed state accepted.");
                return proposedPopulationMixture;
            }

            logger.info("Acceptance rate: " + (double) numAccepted / numSamples);
            return new PopulationMixture(
                    currentPopulationMixture.populationFractions(),
                    currentPopulationMixture.variantProfileCollection(),
                    state.priors().normalPloidyState()
            );
        }
    }
}
