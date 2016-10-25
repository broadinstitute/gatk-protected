package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Iterables;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

final class TumorHeterogeneityStateInitializationUtils {
    private static final int NUM_POPULATIONS_CLONAL = 2;
    private static final int MAX_POPULATION_FRACTIONS_ITERATIONS = 50;

    static TumorHeterogeneityState initializeState(final double initialConcentration,
                                                   final TumorHeterogeneityPriorCollection priors,
                                                   final int numSegments,
                                                   final int numPopulations,
                                                   final int numCells,
                                                   final RandomGenerator rng) {
        //start with Gibbs step
        final boolean doMetropolisStep = false;
        //initialize population fractions to be evenly distributed
        final TumorHeterogeneityState.PopulationFractions initialPopulationFractions =
                new TumorHeterogeneityState.PopulationFractions(Collections.nCopies(numPopulations, 1. / numPopulations));
        //randomly initialize population indicators for each cell
        final TumorHeterogeneityState.PopulationIndicators initialPopulationIndicators =
                initializePopulationIndicators(numCells, initialPopulationFractions, rng);
        //initialize variant profiles to normal
        final int numVariantPopulations = numPopulations - 1;
        final TumorHeterogeneityState.VariantProfileCollection initialVariantProfileCollection =
                initializeNormalProfiles(numVariantPopulations, numSegments, priors.normalPloidyStateIndex());
        return new TumorHeterogeneityState(
                doMetropolisStep, initialConcentration, initialPopulationFractions, initialPopulationIndicators, initialVariantProfileCollection, priors);
    }

    static TumorHeterogeneityState proposeState(final RandomGenerator rng,
                                                final TumorHeterogeneityState state,
                                                final TumorHeterogeneityData data) {
        final int numPopulations = state.numPopulations();
        final int numCells = state.numCells();
        final int numSegments = state.numSegments();
        final boolean doMetropolisStep = true;
        final double concentration = state.concentration();
        //randomly initialize population fractions from prior
        final TumorHeterogeneityState.PopulationFractions populationFractions =
                initializePopulationFractions(numPopulations, concentration, rng);
        //randomly initialize population indicators for each cell
        final TumorHeterogeneityState.PopulationIndicators populationIndicators =
                initializePopulationIndicators(numCells, populationFractions, rng);
        //randomly initialize variant profiles
        final int numVariantPopulations = numPopulations - 1;
        final TumorHeterogeneityPriorCollection priors = state.priors();
        final TumorHeterogeneityState.VariantProfileCollection variantProfileCollection =
                initializeRandomProfiles(numVariantPopulations, numSegments, priors, rng);
        final TumorHeterogeneityState proposedState = new TumorHeterogeneityState(
                doMetropolisStep,concentration, populationFractions, populationIndicators, variantProfileCollection, priors);
        new TumorHeterogeneitySamplers.VariantProfileCollectionSampler(numVariantPopulations, priors.ploidyStatePrior()).sample(rng, proposedState, data);
        return proposedState;
    }

    static TumorHeterogeneityState initializeStateFromClonalResult(final TumorHeterogeneityData data,
                                                                   final TumorHeterogeneityPriorCollection priors,
                                                                   final TumorHeterogeneityModeller clonalModeller,
                                                                   final int maxNumPopulations,
                                                                   final int numCells) {
        final boolean doMetropolisStep = false;
        final double clonalConcentration = Iterables.getLast(clonalModeller.getConcentrationSamples());
        final double clonalNormalFraction = Iterables.getLast(Iterables.getLast(clonalModeller.getPopulationFractionsSamples()));
        final TumorHeterogeneityState.PopulationIndicators initialPopulationIndicators =
                new TumorHeterogeneityState.PopulationIndicators(Iterables.getLast(clonalModeller.getPopulationIndicatorsSamples()));
        IntStream.range(0, numCells).filter(i -> initialPopulationIndicators.get(i) == 1).forEach(i -> initialPopulationIndicators.set(i, maxNumPopulations - 1));
        final TumorHeterogeneityState.VariantProfileCollection clonalVariantProfileCollection = Iterables.getLast(clonalModeller.getVariantProfileCollectionSamples());
        final List<Double> initialFractions = new ArrayList<>();
        initialFractions.add(1. - clonalNormalFraction);
        final List<TumorHeterogeneityState.VariantProfile> initialVariantProfiles = new ArrayList<>();
        initialVariantProfiles.addAll(clonalVariantProfileCollection);
        for (int i = 0; i < maxNumPopulations - NUM_POPULATIONS_CLONAL; i++) {
            initialFractions.add(1, 0.);
            initialVariantProfiles.add(1, TumorHeterogeneityStateInitializationUtils.initializeNormalProfile(data.numSegments(), priors.normalPloidyStateIndex()));
        }
        initialFractions.add(clonalNormalFraction);
        final TumorHeterogeneityState.PopulationFractions initialPopulationFractions =
                new TumorHeterogeneityState.PopulationFractions(initialFractions);
        final TumorHeterogeneityState.VariantProfileCollection initialVariantProfileCollection =
                new TumorHeterogeneityState.VariantProfileCollection(initialVariantProfiles);
        return new TumorHeterogeneityState(doMetropolisStep, clonalConcentration, initialPopulationFractions, initialPopulationIndicators, initialVariantProfileCollection, priors);
    }

    private static TumorHeterogeneityState.PopulationFractions initializePopulationFractions(final int numPopulations,
                                                                                             final double concentration,
                                                                                             final RandomGenerator rng) {
        //sampling from Dirichlet(alpha_vec) is equivalent to sampling from individual Gamma(alpha_vec_i, beta) distributions and normalizing
        final GammaDistribution gammaDistribution = new GammaDistribution(rng, concentration, 1.);
        for (int i = 0; i < MAX_POPULATION_FRACTIONS_ITERATIONS; i++) {
            final double[] unnormalizedPopulationFractions = IntStream.range(0, numPopulations).boxed().mapToDouble(pi -> gammaDistribution.sample()).toArray();
            final double sum = DoubleStream.of(unnormalizedPopulationFractions).sum();
            if (sum > 0) {
                final List<Double> populationFractions = Doubles.asList(MathUtils.normalizeFromRealSpace(unnormalizedPopulationFractions));
                System.out.println(populationFractions);
                return new TumorHeterogeneityState.PopulationFractions(populationFractions);
            }
        }
        System.out.println("Failed to sample population fractions.");
        return new TumorHeterogeneityState.PopulationFractions(Collections.nCopies(numPopulations, 1. / numPopulations));
//        final List<Double> betaSamples = IntStream.range(0, numPopulations).boxed().map(i -> new BetaDistribution(rng, concentration, concentration * (numPopulations - i - 1)).sample()).collect(Collectors.toList());
//        final List<Double> unnormalizedPopulationFractions = new ArrayList<>(numPopulations);
//        for (int populationIndex = 0; populationIndex < numPopulations - 1; populationIndex++) {
//            final double cumSum = unnormalizedPopulationFractions.stream().mapToDouble(Double::doubleValue).sum();
//            unnormalizedPopulationFractions.add((1. - cumSum) * betaSamples.get(populationIndex));
//        }
//        System.out.println(unnormalizedPopulationFractions);
//        final List<Double> populationFractions = Doubles.asList(MathUtils.normalizeFromRealSpace(Doubles.toArray(unnormalizedPopulationFractions)));
    }

    private static TumorHeterogeneityState.PopulationIndicators initializePopulationIndicators(final int numCells,
                                                                                               final TumorHeterogeneityState.PopulationFractions populationFractions,
                                                                                               final RandomGenerator rng) {
        final List<Integer> populationIndices = IntStream.range(0, populationFractions.size()).boxed().collect(Collectors.toList());
        final Function<Integer, Double> probabilityFunction = populationFractions::get;
        return new TumorHeterogeneityState.PopulationIndicators(IntStream.range(0, numCells).boxed()
                .map(p -> GATKProtectedMathUtils.randomSelect(populationIndices, probabilityFunction, rng))
                .collect(Collectors.toList()));
    }

    private static TumorHeterogeneityState.VariantProfileCollection initializeNormalProfiles(final int numVariantPopulations,
                                                                                             final int numSegments,
                                                                                             final int normalPloidyStateIndex) {
        return new TumorHeterogeneityState.VariantProfileCollection(Collections.nCopies(numVariantPopulations, initializeNormalProfile(numSegments, normalPloidyStateIndex)));
    }

    private static TumorHeterogeneityState.VariantProfileCollection initializeRandomProfiles(final int numVariantPopulations,
                                                                                             final int numSegments,
                                                                                             final TumorHeterogeneityPriorCollection priors,
                                                                                             final RandomGenerator rng) {
        return new TumorHeterogeneityState.VariantProfileCollection(
                IntStream.range(0, numVariantPopulations).boxed().map(i -> initializeRandomProfile(numSegments, priors, rng)).collect(Collectors.toList()));
    }

    private static TumorHeterogeneityState.VariantProfile initializeNormalProfile(final int numSegments,
                                                                                  final int normalPloidyStateIndex) {
        final TumorHeterogeneityState.VariantProfile.PloidyStateIndicators ploidyStateIndicators =
                new TumorHeterogeneityState.VariantProfile.PloidyStateIndicators(Collections.nCopies(numSegments, normalPloidyStateIndex));
        return new TumorHeterogeneityState.VariantProfile(ploidyStateIndicators);
    }

    private static TumorHeterogeneityState.VariantProfile initializeRandomProfile(final int numSegments,
                                                                                  final TumorHeterogeneityPriorCollection priors,
                                                                                  final RandomGenerator rng) {
        final PloidyStatePrior ploidyStatePrior = priors.ploidyStatePrior();
        final List<Integer> ploidyStateIndices = IntStream.range(0, ploidyStatePrior.numPloidyStates()).boxed().collect(Collectors.toList());
        final List<Double> ploidyStateProbabilities = ploidyStatePrior.ploidyStates().stream().map(vps -> Math.exp(ploidyStatePrior.logProbability(vps))).collect(Collectors.toList());
        final Function<Integer, Double> probabilityFunction = ploidyStateProbabilities::get;
        final TumorHeterogeneityState.VariantProfile.PloidyStateIndicators ploidyStateIndicators =
                new TumorHeterogeneityState.VariantProfile.PloidyStateIndicators(IntStream.range(0, numSegments).boxed()
                        .map(p -> GATKProtectedMathUtils.randomSelect(ploidyStateIndices, probabilityFunction, rng))
                        .collect(Collectors.toList()));
        return new TumorHeterogeneityState.VariantProfile(ploidyStateIndicators);
    }
}
