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
import java.util.stream.IntStream;

final class TumorHeterogeneityStateInitializationUtils {
    private static final int NUM_POPULATIONS_CLONAL = 2;

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
                initializeNormalProfiles(numVariantPopulations, numSegments);
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
        new TumorHeterogeneitySamplers.VariantProfileCollectionSampler(numVariantPopulations, priors.variantPloidyStatePrior()).sample(rng, proposedState, data);
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
            initialVariantProfiles.add(1, TumorHeterogeneityStateInitializationUtils.initializeNormalProfile(data.numSegments()));
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
        final double[] unnormalizedPopulationFractions = IntStream.range(0, numPopulations).boxed().mapToDouble(i -> gammaDistribution.sample()).toArray();
        final List<Double> populationFractions = Doubles.asList(MathUtils.normalizeFromRealSpace(unnormalizedPopulationFractions));
//        final List<Double> betaSamples = IntStream.range(0, numPopulations).boxed().map(i -> new BetaDistribution(rng, concentration, concentration * (numPopulations - i - 1)).sample()).collect(Collectors.toList());
//        final List<Double> unnormalizedPopulationFractions = new ArrayList<>(numPopulations);
//        for (int populationIndex = 0; populationIndex < numPopulations - 1; populationIndex++) {
//            final double cumSum = unnormalizedPopulationFractions.stream().mapToDouble(Double::doubleValue).sum();
//            unnormalizedPopulationFractions.add((1. - cumSum) * betaSamples.get(populationIndex));
//        }
//        System.out.println(unnormalizedPopulationFractions);
//        final List<Double> populationFractions = Doubles.asList(MathUtils.normalizeFromRealSpace(Doubles.toArray(unnormalizedPopulationFractions)));
        System.out.println(populationFractions);
        return new TumorHeterogeneityState.PopulationFractions(populationFractions);
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
                                                                                             final int numSegments) {
        return new TumorHeterogeneityState.VariantProfileCollection(Collections.nCopies(numVariantPopulations, initializeNormalProfile(numSegments)));
    }

    private static TumorHeterogeneityState.VariantProfileCollection initializeRandomProfiles(final int numVariantPopulations,
                                                                                             final int numSegments,
                                                                                             final TumorHeterogeneityPriorCollection priors,
                                                                                             final RandomGenerator rng) {
        return new TumorHeterogeneityState.VariantProfileCollection(
                IntStream.range(0, numVariantPopulations).boxed().map(i -> initializeRandomProfile(numSegments, priors, rng)).collect(Collectors.toList()));
    }

    private static TumorHeterogeneityState.VariantProfile initializeNormalProfile(final int numSegments) {
        final double variantSegmentFraction = 0.;
        final TumorHeterogeneityState.VariantProfile.VariantIndicators variantIndicators =
                new TumorHeterogeneityState.VariantProfile.VariantIndicators(Collections.nCopies(numSegments, false));
        final TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators variantPloidyStateIndicators =
                new TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators(Collections.nCopies(numSegments, 0));
        return new TumorHeterogeneityState.VariantProfile(variantSegmentFraction, variantIndicators, variantPloidyStateIndicators);
    }

    private static TumorHeterogeneityState.VariantProfile initializeRandomProfile(final int numSegments,
                                                                                  final TumorHeterogeneityPriorCollection priors,
                                                                                  final RandomGenerator rng) {
        final double variantSegmentFraction = new BetaDistribution(rng, priors.variantSegmentFractionPriorAlpha(), priors.variantSegmentFractionPriorBeta()).sample();
        final TumorHeterogeneityState.VariantProfile.VariantIndicators variantIndicators = new TumorHeterogeneityState.VariantProfile.VariantIndicators(
                IntStream.range(0, numSegments).boxed().map(i -> rng.nextDouble() < variantSegmentFraction).collect(Collectors.toList()));
        final PloidyStatePrior variantPloidyStatePrior = priors.variantPloidyStatePrior();
        final List<Integer> variantPloidyStateIndices = IntStream.range(0, variantPloidyStatePrior.numPloidyStates()).boxed().collect(Collectors.toList());
        final List<Double> variantPloidyStateProbabilities = variantPloidyStatePrior.ploidyStates().stream().map(vps -> Math.exp(variantPloidyStatePrior.logProbability(vps))).collect(Collectors.toList());
        final Function<Integer, Double> probabilityFunction = variantPloidyStateProbabilities::get;
        final TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators variantPloidyStateIndicators =
                new TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators(IntStream.range(0, numSegments).boxed()
                        .map(p -> GATKProtectedMathUtils.randomSelect(variantPloidyStateIndices, probabilityFunction, rng))
                        .collect(Collectors.toList()));
        return new TumorHeterogeneityState.VariantProfile(variantSegmentFraction, variantIndicators, variantPloidyStateIndicators);
    }
}
