package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Sets;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Gamma;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.CoordinateUtils;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.SimplexPosition;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.WalkerPosition;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class TumorHeterogeneityUtils {
    private static final Logger logger = LogManager.getLogger(TumorHeterogeneityUtils.class);

    static final double EPSILON = 1E-10;

    //fixes concentration to practically unity for clonal-only version
    static final double CONCENTRATION_MIN = 1. - EPSILON;
    static final double CONCENTRATION_MAX = 1. + EPSILON;


    static final double COPY_RATIO_NOISE_CONSTANT_MIN = EPSILON;
    static final double COPY_RATIO_NOISE_CONSTANT_MAX = 5E-2;

    static final double COPY_RATIO_NOISE_FACTOR_MIN = EPSILON;
    static final double COPY_RATIO_NOISE_FACTOR_MAX = 1.;

    static final double MINOR_ALLELE_FRACTION_NOISE_FACTOR_MIN = EPSILON;
    static final double MINOR_ALLELE_FRACTION_NOISE_FACTOR_MAX = 1.;

    static final double PLOIDY_MIN = EPSILON;

    private static final int CONCENTRATION_WALKER_DIMENSION_INDEX = 0;
    private static final int COPY_RATIO_NOISE_CONSTANT_WALKER_DIMENSION_INDEX = 1;
    private static final int COPY_RATIO_NOISE_FACTOR_WALKER_DIMENSION_INDEX = 2;
    private static final int MINOR_ALLELE_FRACTION_NOISE_FACTOR_WALKER_DIMENSION_INDEX = 3;
    private static final int INITIAL_PLOIDY_WALKER_DIMENSION_INDEX = 4;
    private static final int POPULATION_FRACTIONS_WALKER_DIMENSION_START_INDEX = 5;
    static final int NUM_GLOBAL_PARAMETERS = 5;

    private TumorHeterogeneityUtils() {}

    static double calculateLogPosterior(final TumorHeterogeneityState state,
                                        final TumorHeterogeneityData data) {
        if (isOutsideBounds(state, data)) {
            return Double.NEGATIVE_INFINITY;
        }

        //copy-ratio noise-constant prior
        final double copyRatioNoiseConstantPriorAlpha = data.priors().copyRatioNoiseConstantPriorHyperparameterValues().getAlpha();
        final double copyRatioNoiseConstantPriorBeta = data.priors().copyRatioNoiseConstantPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseConstant = state.copyRatioNoiseConstant();
        final double logPriorCopyRatioNoiseConstant =
                copyRatioNoiseConstantPriorAlpha * Math.log(copyRatioNoiseConstantPriorBeta)
                        + (copyRatioNoiseConstantPriorAlpha - 1.) * Math.log(copyRatioNoiseConstant)
                        - copyRatioNoiseConstantPriorBeta * copyRatioNoiseConstant
                        - Gamma.logGamma(copyRatioNoiseConstantPriorAlpha);

        //copy-ratio noise-factor prior
        final double copyRatioNoiseFactorPriorAlpha = data.priors().copyRatioNoiseFactorPriorHyperparameterValues().getAlpha();
        final double copyRatioNoiseFactorPriorBeta = data.priors().copyRatioNoiseFactorPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseFactor = state.copyRatioNoiseFactor();
        final double logPriorCopyRatioNoiseFactor =
                (copyRatioNoiseFactorPriorAlpha - 1.) * Math.log(Math.max(EPSILON, copyRatioNoiseFactor))
                        + (copyRatioNoiseFactorPriorBeta - 1.) * Math.log(Math.max(EPSILON, 1. - copyRatioNoiseFactor));

        //minor-allele-fraction noise-factor prior
        final double minorAlleleFractionNoiseFactorPriorAlpha = data.priors().minorAlleleFractionNoiseFactorPriorHyperparameterValues().getAlpha();
        final double minorAlleleFractionNoiseFactorPriorBeta = data.priors().minorAlleleFractionNoiseFactorPriorHyperparameterValues().getBeta();
        final double minorAlleleFractionNoiseFactor = state.minorAlleleFractionNoiseFactor();
        final double logPriorMinorAlleleFractionNoiseFactor =
                (minorAlleleFractionNoiseFactorPriorAlpha - 1.) * Math.log(Math.max(EPSILON, minorAlleleFractionNoiseFactor))
                        + (minorAlleleFractionNoiseFactorPriorBeta - 1.) * Math.log(Math.max(EPSILON, 1. - minorAlleleFractionNoiseFactor));

        //variant-profiles prior
        final PopulationMixture.VariantProfileCollection variantProfileCollection = state.populationMixture().variantProfileCollection();
        double logPriorVariantProfiles = 0.;
        for (int populationIndex = 0; populationIndex < variantProfileCollection.numVariantPopulations(); populationIndex++) {
            for (int segmentIndex = 0; segmentIndex < variantProfileCollection.numSegments(); segmentIndex++) {
                final PloidyState ploidyState = variantProfileCollection.ploidyState(populationIndex, segmentIndex);
                logPriorVariantProfiles += data.length(segmentIndex) * data.priors().ploidyStatePrior().logProbability(ploidyState);
            }
        }

        //copy-ratio--minor-allele-fraction likelihood
        double logLikelihoodSegments = 0.;
        final double ploidy = state.ploidy();
        for (int segmentIndex = 0; segmentIndex < data.numSegments(); segmentIndex++) {
            final double mAlleleCopyNumber = state.populationMixture().calculatePopulationAveragedCopyNumberFunction(segmentIndex, PloidyState::m);
            final double nAlleleCopyNumber = state.populationMixture().calculatePopulationAveragedCopyNumberFunction(segmentIndex, PloidyState::n);
            final double totalCopyNumber = mAlleleCopyNumber + nAlleleCopyNumber;
            final double copyRatio = totalCopyNumber / (ploidy + EPSILON);
            final double minorAlleleFraction = calculateMinorAlleleFraction(mAlleleCopyNumber, nAlleleCopyNumber);
            logLikelihoodSegments += data.logDensity(
                    segmentIndex, copyRatio, minorAlleleFraction,
                    state.copyRatioNoiseConstant(), state.copyRatioNoiseFactor(), state.minorAlleleFractionNoiseFactor());
        }

        double logPloidyMismatchPenalty = -data.priors().ploidyMismatchPenalty() * Math.abs(state.initialPloidy() - state.ploidy());

        logger.debug("Log-posterior components:"
                + " " + logPriorCopyRatioNoiseConstant
                + " " + logPriorCopyRatioNoiseFactor
                + " " + logPriorMinorAlleleFractionNoiseFactor
                + " " + logPriorVariantProfiles
                + " " + logLikelihoodSegments
                + " " + logPloidyMismatchPenalty);

        return logPriorCopyRatioNoiseConstant + logPriorCopyRatioNoiseFactor  + logPriorMinorAlleleFractionNoiseFactor
                + logPriorVariantProfiles + logLikelihoodSegments + logPloidyMismatchPenalty;
    }

    static double calculateLogJacobianFactor(final TumorHeterogeneityState state,
                                             final TumorHeterogeneityData data) {
        final double ploidyMax = data.priors().ploidyStatePrior().maxCopyNumber();
        return CoordinateUtils.calculateLogJacobianFactor(state.concentration(), CONCENTRATION_MIN, CONCENTRATION_MAX)
                + CoordinateUtils.calculateLogJacobianFactor(state.copyRatioNoiseConstant(), COPY_RATIO_NOISE_CONSTANT_MIN, COPY_RATIO_NOISE_CONSTANT_MAX)
                + CoordinateUtils.calculateLogJacobianFactor(state.copyRatioNoiseFactor(), COPY_RATIO_NOISE_FACTOR_MIN, COPY_RATIO_NOISE_FACTOR_MAX)
                + CoordinateUtils.calculateLogJacobianFactor(state.minorAlleleFractionNoiseFactor(), MINOR_ALLELE_FRACTION_NOISE_FACTOR_MIN, MINOR_ALLELE_FRACTION_NOISE_FACTOR_MAX)
                + CoordinateUtils.calculateLogJacobianFactor(state.populationMixture().ploidy(data), PLOIDY_MIN, ploidyMax)
                + SimplexPosition.calculateLogJacobianFactor(state.populationMixture().populationFractions());
    }

    static TumorHeterogeneityState transformWalkerPositionToState(final WalkerPosition walkerPosition,
                                                                  final RandomGenerator rng,
                                                                  final TumorHeterogeneityData data,
                                                                  final List<List<Integer>> totalCopyNumberProductStates,
                                                                  final Map<Integer, Set<PloidyState>> ploidyStateSetsMap) {
        final double concentration = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(CONCENTRATION_WALKER_DIMENSION_INDEX), CONCENTRATION_MIN, CONCENTRATION_MAX);
        final double copyRatioNoiseConstant = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(COPY_RATIO_NOISE_CONSTANT_WALKER_DIMENSION_INDEX), COPY_RATIO_NOISE_CONSTANT_MIN, COPY_RATIO_NOISE_CONSTANT_MAX);
        final double copyRatioNoiseFactor = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(COPY_RATIO_NOISE_FACTOR_WALKER_DIMENSION_INDEX), COPY_RATIO_NOISE_FACTOR_MIN, COPY_RATIO_NOISE_FACTOR_MAX);
        final double minorAlleleFractionNoiseFactor = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(MINOR_ALLELE_FRACTION_NOISE_FACTOR_WALKER_DIMENSION_INDEX), MINOR_ALLELE_FRACTION_NOISE_FACTOR_MIN, MINOR_ALLELE_FRACTION_NOISE_FACTOR_MAX);
        final double ploidyMax = data.priors().ploidyStatePrior().maxCopyNumber();
        final double initialPloidy = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(INITIAL_PLOIDY_WALKER_DIMENSION_INDEX), PLOIDY_MIN, ploidyMax);

        final PopulationMixture.PopulationFractions populationFractions = new PopulationMixture.PopulationFractions(
                SimplexPosition.calculateSimplexPositionFromWalkerPosition(
                        new WalkerPosition(walkerPosition.subList(POPULATION_FRACTIONS_WALKER_DIMENSION_START_INDEX, walkerPosition.numDimensions()))));

        final PopulationMixture.VariantProfileCollection variantProfileCollection = proposeVariantProfileCollection(
                rng, copyRatioNoiseConstant, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, initialPloidy,
                populationFractions, data, totalCopyNumberProductStates, ploidyStateSetsMap);

        final PopulationMixture populationMixture = new PopulationMixture(populationFractions, variantProfileCollection, data.priors().normalPloidyState());
        final double ploidy = populationMixture.ploidy(data);

        return new TumorHeterogeneityState(concentration, copyRatioNoiseConstant, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, initialPloidy, ploidy, populationMixture);
    }

    static WalkerPosition transformStateToWalkerPosition(final TumorHeterogeneityState state,
                                                         final TumorHeterogeneityData data) {
        final double concentrationWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.concentration(), CONCENTRATION_MIN, CONCENTRATION_MAX);
        final double copyRatioNoiseConstantWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.copyRatioNoiseConstant(), COPY_RATIO_NOISE_CONSTANT_MIN, COPY_RATIO_NOISE_CONSTANT_MAX);
        final double copyRatioNoiseFactorWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.copyRatioNoiseFactor(), COPY_RATIO_NOISE_FACTOR_MIN, COPY_RATIO_NOISE_FACTOR_MAX);
        final double minorAlleleFractionNoiseFactorWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.minorAlleleFractionNoiseFactor(), MINOR_ALLELE_FRACTION_NOISE_FACTOR_MIN, MINOR_ALLELE_FRACTION_NOISE_FACTOR_MAX);
        final double ploidyMax = data.priors().ploidyStatePrior().maxCopyNumber();
        final double initialPloidyWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.initialPloidy(), PLOIDY_MIN, ploidyMax);
        final WalkerPosition populationFractionsWalkerCoordinates = SimplexPosition.calculateWalkerPositionFromSimplexPosition(state.populationMixture().populationFractions());

        return new WalkerPosition(ListUtils.union(
                Arrays.asList(
                        concentrationWalkerCoordinate,
                        copyRatioNoiseConstantWalkerCoordinate,
                        copyRatioNoiseFactorWalkerCoordinate,
                        minorAlleleFractionNoiseFactorWalkerCoordinate,
                        initialPloidyWalkerCoordinate),
                populationFractionsWalkerCoordinates));
    }

    private static PopulationMixture.VariantProfileCollection proposeVariantProfileCollection(final RandomGenerator rng,
                                                                                              final double copyRatioNoiseConstant,
                                                                                              final double copyRatioNoiseFactor,
                                                                                              final double minorAlleleFractionNoiseFactor,
                                                                                              final double initialPloidy,
                                                                                              final PopulationMixture.PopulationFractions populationFractions,
                                                                                              final TumorHeterogeneityData data,
                                                                                              final List<List<Integer>> totalCopyNumberProductStates,
                                                                                              final Map<Integer, Set<PloidyState>> ploidyStateSetsMap) {
        final int numPopulations = populationFractions.numPopulations();
        final int numSegments = data.numSegments();
        final PloidyState normalPloidyState = data.priors().normalPloidyState();
        final PloidyStatePrior ploidyStatePrior = data.priors().ploidyStatePrior();
        final List<PopulationMixture.VariantProfile> variantProfiles = new ArrayList<>(
                Collections.nCopies(numPopulations - 1, new PopulationMixture.VariantProfile(Collections.nCopies(numSegments, normalPloidyState))));

        for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
            final int si = segmentIndex;
            final double[] log10ProbabilitiesCopyRatio = totalCopyNumberProductStates.stream()
                    .mapToDouble(tcnps -> calculateTotalCopyNumber(populationFractions, tcnps, normalPloidyState) / (initialPloidy + EPSILON))
                    .map(cr -> data.copyRatioLogDensity(si, cr, copyRatioNoiseConstant, copyRatioNoiseFactor))
                    .map(MathUtils::logToLog10)
                    .toArray();
            final double[] probabilitiesCopyRatio = MathUtils.normalizeFromLog10ToLinearSpace(log10ProbabilitiesCopyRatio);
            final Map<List<Integer>, Double> probabilityMapCopyRatio = new HashMap<>(totalCopyNumberProductStates.size());
            IntStream.range(0, totalCopyNumberProductStates.size()).forEach(i -> probabilityMapCopyRatio.put(totalCopyNumberProductStates.get(i), probabilitiesCopyRatio[i]));
            final Function<List<Integer>, Double> probabilityFunctionCopyRatio = probabilityMapCopyRatio::get;
            final List<Integer> totalCopyNumberProductState = GATKProtectedMathUtils.randomSelect(totalCopyNumberProductStates, probabilityFunctionCopyRatio, rng);
            final double totalCopyRatio = calculateTotalCopyNumber(populationFractions, totalCopyNumberProductState, normalPloidyState) / (initialPloidy + EPSILON);

            final List<List<PloidyState>> ploidyStateProductStates =
                    new ArrayList<>(Sets.cartesianProduct(totalCopyNumberProductState.stream().map(ploidyStateSetsMap::get).collect(Collectors.toList())));
            final double[] log10ProbabilitiesPloidyStateProductStates = ploidyStateProductStates.stream()
                    .mapToDouble(psps ->
                            data.logDensity(si, totalCopyRatio, calculateMinorAlleleFraction(populationFractions, psps, normalPloidyState), copyRatioNoiseConstant, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor)
                                    + psps.stream().mapToDouble(ps -> data.length(si) * ploidyStatePrior.logProbability(ps)).sum())
                    .map(MathUtils::logToLog10)
                    .toArray();
            final double[] probabilitiesPloidyStateProductStates = MathUtils.normalizeFromLog10ToLinearSpace(log10ProbabilitiesPloidyStateProductStates);
            final Map<List<PloidyState>, Double> probabilityMapPloidyStateProductStates = new HashMap<>(ploidyStateProductStates.size());
            IntStream.range(0, ploidyStateProductStates.size()).forEach(i -> probabilityMapPloidyStateProductStates.put(ploidyStateProductStates.get(i), probabilitiesPloidyStateProductStates[i]));
            final Function<List<PloidyState>, Double> probabilityFunctionPloidyStateProductStates = probabilityMapPloidyStateProductStates::get;
            final List<PloidyState> ploidyStateProductState = GATKProtectedMathUtils.randomSelect(ploidyStateProductStates, probabilityFunctionPloidyStateProductStates, rng);

            IntStream.range(0, numPopulations - 1).forEach(i -> variantProfiles.get(i).set(si, ploidyStateProductState.get(i)));
        }
        return new PopulationMixture.VariantProfileCollection(variantProfiles);
    }

    private static double calculateTotalCopyNumber(final PopulationMixture.PopulationFractions populationFractions,
                                                   final List<Integer> totalCopyNumberProductState,
                                                   final PloidyState normalPloidyState) {
        final int numPopulations = populationFractions.numPopulations();
        return IntStream.range(0, numPopulations - 1).boxed()
                .mapToDouble(i -> totalCopyNumberProductState.get(i) * populationFractions.get(i))
                .sum() + normalPloidyState.total() * populationFractions.normalFraction();
    }

    private static double calculateMinorAlleleFraction(final PopulationMixture.PopulationFractions populationFractions,
                                                       final List<PloidyState> ploidyStateProductState,
                                                       final PloidyState normalPloidyState) {
        final int numPopulations = populationFractions.numPopulations();
        final double mAlleleCopyNumber = IntStream.range(0, numPopulations - 1).boxed()
                .mapToDouble(i -> ploidyStateProductState.get(i).m() * populationFractions.get(i))
                .sum() + normalPloidyState.m() * populationFractions.normalFraction();
        final double nAlleleCopyNumber = IntStream.range(0, numPopulations - 1).boxed()
                .mapToDouble(i -> ploidyStateProductState.get(i).n() * populationFractions.get(i))
                .sum() + normalPloidyState.n() * populationFractions.normalFraction();
        return Math.min(mAlleleCopyNumber, nAlleleCopyNumber) / (mAlleleCopyNumber + nAlleleCopyNumber + EPSILON);
    }

    private static double calculateMinorAlleleFraction(final double m, final double n) {
        return Math.min(m, n) / (m + n + EPSILON);
    }

    private static boolean isOutsideBounds(final TumorHeterogeneityState state,
                                           final TumorHeterogeneityData data) {
        final double ploidyMax = data.priors().ploidyStatePrior().maxCopyNumber();
        return state.concentration() < CONCENTRATION_MIN || state.concentration() > CONCENTRATION_MAX ||
                state.copyRatioNoiseConstant() < COPY_RATIO_NOISE_CONSTANT_MIN || state.copyRatioNoiseConstant() > COPY_RATIO_NOISE_CONSTANT_MAX ||
                state.copyRatioNoiseFactor() < COPY_RATIO_NOISE_FACTOR_MIN || state.copyRatioNoiseFactor() > COPY_RATIO_NOISE_FACTOR_MAX ||
                state.minorAlleleFractionNoiseFactor() < MINOR_ALLELE_FRACTION_NOISE_FACTOR_MIN || state.minorAlleleFractionNoiseFactor() > MINOR_ALLELE_FRACTION_NOISE_FACTOR_MAX ||
                state.ploidy() < PLOIDY_MIN || state.ploidy() > ploidyMax ||
                state.populationMixture().variantProfileCollection().stream().anyMatch(PopulationMixture.VariantProfile::isCompleteDeletion);
    }
}
