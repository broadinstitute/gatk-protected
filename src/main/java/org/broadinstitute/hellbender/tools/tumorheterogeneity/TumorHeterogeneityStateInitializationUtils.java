package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Iterables;
import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;

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
        //initialize population fractions to be evenly distributed
        final TumorHeterogeneityState.PopulationFractions initialPopulationFractions =
                new TumorHeterogeneityState.PopulationFractions(Collections.nCopies(numPopulations, 1. / numPopulations));
        //randomly initialize population indicators for each cell
        final TumorHeterogeneityState.PopulationIndicators initialPopulationIndicators =
                initializePopulationIndicators(initialPopulationFractions, numPopulations, numCells, rng);
        //initialize variant profiles
        final int numVariantPopulations = numPopulations - 1;
        final TumorHeterogeneityState.VariantProfileCollection initialVariantProfileCollection =
                initializeProfiles(numVariantPopulations, numSegments, priors.variantPloidyStatePrior(), rng);
        return new TumorHeterogeneityState(
                initialConcentration, initialPopulationFractions, initialPopulationIndicators, initialVariantProfileCollection, priors);
    }

    private static TumorHeterogeneityState.PopulationIndicators initializePopulationIndicators(final TumorHeterogeneityState.PopulationFractions initialPopulationFractions,
                                                                                               final int numPopulations,
                                                                                               final int numCells,
                                                                                               final RandomGenerator rng) {
        final List<Integer> populationIndices = IntStream.range(0, numPopulations).boxed().collect(Collectors.toList());
        final Function<Integer, Double> probabilityFunction = initialPopulationFractions::get;
        return new TumorHeterogeneityState.PopulationIndicators(IntStream.range(0, numCells).boxed()
                .map(p -> GATKProtectedMathUtils.randomSelect(populationIndices, probabilityFunction, rng))
                .collect(Collectors.toList()));
    }

    private static TumorHeterogeneityState.VariantProfileCollection initializeProfiles(final int numVariantPopulations,
                                                                                       final int numSegments,
                                                                                       final PloidyStatePrior variantPloidyStatePrior,
                                                                                       final RandomGenerator rng) {
//        return new TumorHeterogeneityState.VariantProfileCollection(Collections.nCopies(numVariantPopulations, initializeNormalProfile(numSegments)));
        return new TumorHeterogeneityState.VariantProfileCollection(Collections.nCopies(numVariantPopulations, initializeRandomProfile(numSegments, variantPloidyStatePrior, rng)));
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
                                                                                  final PloidyStatePrior variantPloidyStatePrior,
                                                                                  final RandomGenerator rng) {
        final double variantSegmentFraction = 1.;
        final TumorHeterogeneityState.VariantProfile.VariantIndicators variantIndicators =
                new TumorHeterogeneityState.VariantProfile.VariantIndicators(Collections.nCopies(numSegments, true));
        final List<Integer> variantPloidyStateIndices = IntStream.range(0, variantPloidyStatePrior.numPloidyStates()).boxed().collect(Collectors.toList());
        final List<Double> variantPloidyStateProbabilities = variantPloidyStatePrior.ploidyStates().stream().map(vps -> Math.exp(variantPloidyStatePrior.logProbability(vps))).collect(Collectors.toList());
        final Function<Integer, Double> probabilityFunction = variantPloidyStateProbabilities::get;
        final TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators variantPloidyStateIndicators =
                new TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators(IntStream.range(0, numSegments).boxed()
                        .map(p -> GATKProtectedMathUtils.randomSelect(variantPloidyStateIndices, probabilityFunction, rng))
                        .collect(Collectors.toList()));
        return new TumorHeterogeneityState.VariantProfile(variantSegmentFraction, variantIndicators, variantPloidyStateIndicators);
    }

    static TumorHeterogeneityState initializeStateFromClonalResult(final TumorHeterogeneityData data,
                                                                   final TumorHeterogeneityPriorCollection priors,
                                                                   final TumorHeterogeneityModeller clonalModeller,
                                                                   final int maxNumPopulations,
                                                                   final int numCells) {
        final double clonalConcentration = Iterables.getLast(clonalModeller.getConcentrationSamples());
        final double clonalNormalFraction = Iterables.getLast(Iterables.getLast(clonalModeller.getPopulationFractionsSamples()));
        final TumorHeterogeneityState.PopulationIndicators initialPopulationIndicators = new TumorHeterogeneityState.PopulationIndicators(Iterables.getLast(clonalModeller.getPopulationIndicatorsSamples()));
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
        return new TumorHeterogeneityState(clonalConcentration, initialPopulationFractions, initialPopulationIndicators, initialVariantProfileCollection, priors);
    }
}
