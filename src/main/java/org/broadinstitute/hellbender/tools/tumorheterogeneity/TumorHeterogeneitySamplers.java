package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import autovalue.shaded.com.google.common.common.collect.Sets;
import org.apache.commons.math3.analysis.function.Logit;
import org.apache.commons.math3.analysis.function.Sigmoid;
import org.apache.commons.math3.distribution.NormalDistribution;
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
    private static final double COPY_RATIO_NOISE_FLOOR_MAX = 1E-2;

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
                return TumorHeterogeneityLikelihoods.calculateLogPosterior(newState, data);
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
                return TumorHeterogeneityLikelihoods.calculateLogPosterior(newState, data);
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
                return TumorHeterogeneityLikelihoods.calculateLogPosterior(newState, data);
            };
            final double minorAlleleFractionNoiseFactor = new SliceSampler(rng, logConditionalPDF, MIN, MAX, minorAlleleFractionNoiseFactorSliceSamplingWidth).sample(state.minorAlleleFractionNoiseFactor());
            logger.debug("Sampled MAF noise factor: " + minorAlleleFractionNoiseFactor);
            return minorAlleleFractionNoiseFactor;
        }
    }

    static final class PopulationMixtureSampler implements ParameterSampler<PopulationMixture, TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> {
        private static final int MAX_NUM_PLOIDY_STEP_ITERATIONS = 25;
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
                    calculateTransformedPopulationFractionsFromPopulationFractions(currentPopulationFractions);
            final double currentPloidy = currentPopulationMixture.ploidy(data);
            logger.info("Current population fractions: " + currentPopulationFractions);
            logger.info("Current ploidy: " + currentPloidy);

            final List<Double> proposedTransformedPopulationFractions = currentTransformedPopulationFractions.stream()
                    .map(x -> proposeTransformedPopulationFraction(rng, x))
                    .collect(Collectors.toList());
            final PopulationMixture.PopulationFractions proposedPopulationFractions =
                    calculatePopulationFractionsFromTransformedPopulationFractions(proposedTransformedPopulationFractions);

            final double proposedPloidy = proposePloidy(rng, currentPloidy, maxTotalCopyNumber);
            logger.info("Proposed initial ploidy: " + proposedPloidy);

            final PopulationMixture.VariantProfileCollection proposedVariantProfileCollection =
                    proposeVariantProfileCollection(rng, state, data, proposedPopulationFractions, proposedPloidy,
                            totalCopyNumberProductStates, ploidyStateSetsMap);
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

            final double proposedLogPosterior = TumorHeterogeneityLikelihoods.calculateLogPosterior(proposedState, data) + calculateLogJacobianFactor(proposedPopulationFractions);
            final double currentLogPosterior = TumorHeterogeneityLikelihoods.calculateLogPosterior(state, data) + calculateLogJacobianFactor(currentPopulationFractions);
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

        private static PopulationMixture.VariantProfileCollection proposeVariantProfileCollection(final RandomGenerator rng,
                                                                                                  final TumorHeterogeneityState currentState,
                                                                                                  final TumorHeterogeneityData data,
                                                                                                  final PopulationMixture.PopulationFractions proposedPopulationFractions,
                                                                                                  final double proposedPloidy,
                                                                                                  final List<List<Integer>> totalCopyNumberProductStates,
                                                                                                  final Map<Integer, Set<PloidyState>> ploidyStateSetsMap) {
            final int numPopulations = currentState.populationMixture().numPopulations();
            final int numSegments = data.numSegments();
            final PloidyState normalPloidyState = currentState.priors().normalPloidyState();
            final List<PopulationMixture.VariantProfile> variantProfiles = new ArrayList<>(Collections.nCopies(numPopulations - 1,
                    new PopulationMixture.VariantProfile(Collections.nCopies(numSegments, normalPloidyState))));

            for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
                final int si = segmentIndex;
                final double[] log10ProbabilitiesCopyRatio = totalCopyNumberProductStates.stream()
                        .mapToDouble(tcnps -> calculateTotalCopyNumber(proposedPopulationFractions, tcnps, normalPloidyState) / proposedPloidy)
                        .map(cr -> data.copyRatioLogDensity(si, cr, currentState.copyRatioNoiseFloor(), currentState.copyRatioNoiseFactor()))
                        .map(MathUtils::logToLog10)
                        .toArray();
                final double[] probabilitiesCopyRatio = MathUtils.normalizeFromLog10ToLinearSpace(log10ProbabilitiesCopyRatio);
                final Function<List<Integer>, Double> probabilityFunctionCopyRatio = totalCopyNumberProductState ->
                        probabilitiesCopyRatio[totalCopyNumberProductStates.indexOf(totalCopyNumberProductState)];
                final List<Integer> totalCopyNumberProductState = GATKProtectedMathUtils.randomSelect(totalCopyNumberProductStates, probabilityFunctionCopyRatio, rng);
                final double totalCopyRatio = calculateTotalCopyNumber(proposedPopulationFractions, totalCopyNumberProductState, normalPloidyState) / proposedPloidy;

                final List<List<PloidyState>> ploidyStateProductStates =
                        new ArrayList<>(Sets.cartesianProduct(totalCopyNumberProductState.stream().map(ploidyStateSetsMap::get).collect(Collectors.toList())));
                final double[] log10Probabilities = ploidyStateProductStates.stream()
                        .mapToDouble(ps -> calculateMinorAlleleFraction(proposedPopulationFractions, ps, normalPloidyState))
                        .map(maf -> data.logDensity(si, totalCopyRatio, maf, currentState.copyRatioNoiseFloor(), currentState.copyRatioNoiseFactor(), currentState.minorAlleleFractionNoiseFactor()))
                        .map(MathUtils::logToLog10)
                        .toArray();
                final double[] probabilities = MathUtils.normalizeFromLog10ToLinearSpace(log10Probabilities);
                final Function<List<PloidyState>, Double> probabilityFunction = ploidyStateProductState ->
                        probabilities[ploidyStateProductStates.indexOf(ploidyStateProductState)];
                final List<PloidyState> ploidyStateProductState = GATKProtectedMathUtils.randomSelect(ploidyStateProductStates, probabilityFunction, rng);

                IntStream.range(0, numPopulations - 1).forEach(i -> variantProfiles.get(i).set(si, ploidyStateProductState.get(i)));
            }
            return new PopulationMixture.VariantProfileCollection(variantProfiles);
        }

        private static double calculateTotalCopyNumber(final PopulationMixture.PopulationFractions populationFractions,
                                                       final List<Integer> totalCopyNumberProductState,
                                                       final PloidyState normalPloidyState) {
            final int numPopulations = populationFractions.size();
            return IntStream.range(0, numPopulations - 1).boxed()
                    .mapToDouble(i -> totalCopyNumberProductState.get(i) * populationFractions.get(i))
                    .sum() + normalPloidyState.total() * populationFractions.get(numPopulations - 1);
        }

        private static double calculateMinorAlleleFraction(final PopulationMixture.PopulationFractions populationFractions,
                                                           final List<PloidyState> ploidyStateProductState,
                                                           final PloidyState normalPloidyState) {
            final int numPopulations = populationFractions.size();
            final double mAlleleCopyNumber = IntStream.range(0, numPopulations - 1).boxed()
                    .mapToDouble(i -> ploidyStateProductState.get(i).m() * populationFractions.get(i))
                    .sum() + normalPloidyState.m() * populationFractions.get(numPopulations - 1);
            final double nAlleleCopyNumber = IntStream.range(0, numPopulations - 1).boxed()
                    .mapToDouble(i -> ploidyStateProductState.get(i).n() * populationFractions.get(i))
                    .sum() + normalPloidyState.n() * populationFractions.get(numPopulations - 1);
            return Math.min(mAlleleCopyNumber, nAlleleCopyNumber) / (mAlleleCopyNumber + nAlleleCopyNumber + EPSILON);
        }

        private static double calculateLogJacobianFactor(final PopulationMixture.PopulationFractions populationFractions) {
            final List<Double> breakProportions = calculateBreakProportionsFromPopulationFractions(populationFractions);
            return IntStream.range(0, populationFractions.size() - 1).boxed()
                    .mapToDouble(i -> Math.log(populationFractions.get(i)) + Math.log(1. - breakProportions.get(i))).sum();

        }

        private static List<Double> calculateTransformedPopulationFractionsFromPopulationFractions(final PopulationMixture.PopulationFractions populationFractions) {
            final List<Double> breakProportions = calculateBreakProportionsFromPopulationFractions(populationFractions);
            return calculateTransformedPopulationFractionsFromBreakProportions(breakProportions);
        }

        private static PopulationMixture.PopulationFractions calculatePopulationFractionsFromTransformedPopulationFractions(final List<Double> transformedPopulationFractions) {
            final List<Double> breakProportions = calculateBreakProportionsFromTransformedPopulationFractions(transformedPopulationFractions);
            final List<Double> populationFractions = calculatePopulationFractionsFromBreakProportions(breakProportions);
            return new PopulationMixture.PopulationFractions(populationFractions);
        }

        private static List<Double> calculatePopulationFractionsFromBreakProportions(final List<Double> breakProportions) {
            final int numPopulations = breakProportions.size() + 1;
            final List<Double> populationFractions = new ArrayList<>();
            double cumulativeSum = 0.;
            for (int populationIndex = 0; populationIndex < numPopulations - 1; populationIndex++) {
                final double populationFraction = (1. - cumulativeSum) * breakProportions.get(populationIndex);
                populationFractions.add(populationFraction);
                cumulativeSum += populationFraction;
            }
            populationFractions.add(1. - cumulativeSum);
            return new PopulationMixture.PopulationFractions(populationFractions);
        }

        private static List<Double> calculateBreakProportionsFromPopulationFractions(final PopulationMixture.PopulationFractions populationFractions) {
            final int numPopulations = populationFractions.size();
            final List<Double> breakProportions = new ArrayList<>();
            double cumulativeSum = 0.;
            for (int populationIndex = 0; populationIndex < numPopulations - 1; populationIndex++) {
                final double breakProportion = populationFractions.get(populationIndex) / (1. - cumulativeSum);
                breakProportions.add(breakProportion);
                cumulativeSum += populationFractions.get(populationIndex);
            }
            return breakProportions;
        }

        private static List<Double> calculateBreakProportionsFromTransformedPopulationFractions(final List<Double> transformedPopulationFractions) {
            final int numPopulations = transformedPopulationFractions.size() + 1;
            return IntStream.range(0, numPopulations - 1).boxed()
                    .map(i -> new Sigmoid().value(transformedPopulationFractions.get(i) + Math.log(1. / (numPopulations - (i + 1)))))
                    .collect(Collectors.toList());
        }

        private static List<Double> calculateTransformedPopulationFractionsFromBreakProportions(final List<Double> breakProportions) {
            final int numPopulations = breakProportions.size() + 1;
            return IntStream.range(0, numPopulations - 1).boxed()
                    .map(i -> new Logit().value(breakProportions.get(i)) - Math.log(1. / (numPopulations - (i + 1))))
                    .collect(Collectors.toList());
        }

        private static double proposeTransformedPopulationFraction(final RandomGenerator rng, final Double currentTransformedPopulationFraction) {
            return rng.nextDouble() < 0.5
                    ? currentTransformedPopulationFraction + new NormalDistribution(rng, 0., transformedPopulationFractionProposalWidth).sample()
                    : currentTransformedPopulationFraction + new NormalDistribution(rng, 0., 10 * transformedPopulationFractionProposalWidth).sample();
        }

        private static double proposePloidy(final RandomGenerator rng,
                                            final double currentPloidy,
                                            final int maxTotalCopyNumber) {
            int numIterations = 0;
            final NormalDistribution normal = rng.nextDouble() < 0.5
                    ? new NormalDistribution(rng, 0., ploidyProposalWidth)
                    : new NormalDistribution(rng, 0., 10 * ploidyProposalWidth);
            while (numIterations < MAX_NUM_PLOIDY_STEP_ITERATIONS) {
                final double proposedPloidy = currentPloidy + normal.sample();
                if (0 < proposedPloidy && proposedPloidy <= maxTotalCopyNumber) {
                    return proposedPloidy;
                }
                numIterations++;
            }
            return currentPloidy;
        }
    }
}
