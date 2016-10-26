package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Iterables;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

final class TumorHeterogeneityStateInitializationUtils {
    private static final int NUM_POPULATIONS_CLONAL = 2;
    private static final int MAX_POPULATION_FRACTIONS_ITERATIONS = 50;

    private static final Logger logger = LogManager.getLogger(TumorHeterogeneityStateInitializationUtils.class);

    /**
     * Initialize state with evenly distributed population fractions and normal variant profiles.
     */
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

    /**
     * Propose a state for a Metropolis step.  Population fractions are sampled from the Dirichlet prior distribution
     * and corresponding variant profiles are then sampled in the same way as in the Gibbs steps.
     */
    static TumorHeterogeneityState proposeState(final RandomGenerator rng,
                                                final TumorHeterogeneityState state,
                                                final TumorHeterogeneityData data) {
        final int numPopulations = state.numPopulations();
        final int numCells = state.numCells();
        final int numSegments = state.numSegments();
        //record that we are performing a Metropolis step
        final boolean doMetropolisStep = true;
        final double concentration = state.concentration();
        //randomly initialize population fractions from prior
        final TumorHeterogeneityState.PopulationFractions populationFractions =
                initializePopulationFractions(numPopulations, concentration, rng);
        //randomly initialize population indicators for each cell
        final TumorHeterogeneityState.PopulationIndicators populationIndicators =
                initializePopulationIndicators(numCells, populationFractions, rng);
        //initialize variant profiles using sampler
        final int numVariantPopulations = numPopulations - 1;
        final TumorHeterogeneityPriorCollection priors = state.priors();
        final TumorHeterogeneityState.VariantProfileCollection variantProfileCollection =
                initializeNormalProfiles(numVariantPopulations, numSegments, priors.normalPloidyStateIndex());
        final TumorHeterogeneityState proposedState = new TumorHeterogeneityState(
                doMetropolisStep,concentration, populationFractions, populationIndicators, variantProfileCollection, priors);
        new TumorHeterogeneitySamplers.VariantProfileCollectionSampler(numVariantPopulations, priors.ploidyStatePrior()).sample(rng, proposedState, data);
        return proposedState;
    }

    /**
     * Initialize state from a modeller run that assumes one clonal population and one normal population.
     */
    static TumorHeterogeneityState initializeStateFromClonalResult(final TumorHeterogeneityData data,
                                                                   final TumorHeterogeneityPriorCollection priors,
                                                                   final TumorHeterogeneityModeller clonalModeller,
                                                                   final int maxNumPopulations,
                                                                   final int numCells) {
        //double check some of the parameters
        Utils.validateArg(clonalModeller.getConcentrationSamples().size() > 0,
                "Clonal modeller must have a non-zero number of samples.");
        Utils.validateArg(clonalModeller.getPopulationFractionsSamples().get(0).size() == 2,
                "Clonal modeller must have two populations (clone + normal).");
        //start with Gibbs step
        final boolean doMetropolisStep = false;
        //initialize concentration to posterior mean of clonal result
        final double clonalConcentration = new Mean().evaluate(Doubles.toArray(clonalModeller.getConcentrationSamples()));

        //initialize normal fraction to posterior mean of clonal result
        final double[] normalFractionSamples = clonalModeller.getPopulationFractionsSamples().stream()
                .mapToDouble(pfs -> pfs.get(NUM_POPULATIONS_CLONAL - 1)).toArray();
        final double clonalNormalFraction = new Mean().evaluate(normalFractionSamples);

        //initialize population indicators to last state of clonal result
        final TumorHeterogeneityState.PopulationIndicators initialPopulationIndicators =
                new TumorHeterogeneityState.PopulationIndicators(Iterables.getLast(clonalModeller.getPopulationIndicatorsSamples()));
        //reset population index of normal-population indicators to maxNumPopulations - 1 (instead of NUM_POPULATIONS_CLONAL - 1 = 1)
        IntStream.range(0, numCells).filter(i -> initialPopulationIndicators.get(i) == 1).forEach(i -> initialPopulationIndicators.set(i, maxNumPopulations - 1));

        //initialize one variant profile to posterior mode of clonal result
        final List<TumorHeterogeneityState.VariantProfile> clonalProfileSamples = clonalModeller.getVariantProfileCollectionSamples().stream()
                .map(vpc -> vpc.get(0)).collect(Collectors.toList());
        final TumorHeterogeneityState.VariantProfile initialClonalProfile = calculateVariantProfilePosteriorMode(clonalProfileSamples);

        //build new population fractions and variant-profile collection with additional variant populations
        final List<Double> initialFractions = new ArrayList<>();
        final List<TumorHeterogeneityState.VariantProfile> initialVariantProfiles = new ArrayList<>();
        //add clonal population
        initialFractions.add(1. - clonalNormalFraction);
        initialVariantProfiles.add(initialClonalProfile);
        //initialize additional variant profiles with zero population fraction and normal profile
        for (int i = 0; i < maxNumPopulations - NUM_POPULATIONS_CLONAL; i++) {
            initialFractions.add(1, 0.);
            initialVariantProfiles.add(1, TumorHeterogeneityStateInitializationUtils.initializeNormalProfile(data.numSegments(), priors.normalPloidyStateIndex()));
        }
        //add normal population fraction
        initialFractions.add(clonalNormalFraction);
        final TumorHeterogeneityState.PopulationFractions initialPopulationFractions =
                new TumorHeterogeneityState.PopulationFractions(initialFractions);
        final TumorHeterogeneityState.VariantProfileCollection initialVariantProfileCollection =
                new TumorHeterogeneityState.VariantProfileCollection(initialVariantProfiles);

        return new TumorHeterogeneityState(doMetropolisStep, clonalConcentration, initialPopulationFractions, initialPopulationIndicators, initialVariantProfileCollection, priors);
    }

    /**
     * Sample population fractions from a symmetric Dirichlet distribution.
     */
    private static TumorHeterogeneityState.PopulationFractions initializePopulationFractions(final int numPopulations,
                                                                                             final double concentration,
                                                                                             final RandomGenerator rng) {
        //sampling from Dirichlet(alpha_vec) is equivalent to sampling from individual Gamma(alpha_vec_i, 1) distributions and normalizing
        final GammaDistribution gammaDistribution = new GammaDistribution(rng, concentration, 1.);
        //for small values of concentration, all Gamma(alpha_vec_i, 1) samples may be zero,
        //so we iterate until this is not the case (unless a maximum number of iterations is reached)
        for (int i = 0; i < MAX_POPULATION_FRACTIONS_ITERATIONS; i++) {
            final double[] unnormalizedPopulationFractions = IntStream.range(0, numPopulations).boxed().mapToDouble(pi -> gammaDistribution.sample()).toArray();
            final double sum = DoubleStream.of(unnormalizedPopulationFractions).sum();
            if (sum > 0) {
                final List<Double> populationFractions = Doubles.asList(MathUtils.normalizeFromRealSpace(unnormalizedPopulationFractions));
                logger.debug("Proposed population fractions: " + populationFractions);
                return new TumorHeterogeneityState.PopulationFractions(populationFractions);
            }
        }
        //if non-zero population fractions were not sampled, emit warning and return evenly distributed fractions
        logger.warn("Failed to sample non-zero population fractions within " + MAX_POPULATION_FRACTIONS_ITERATIONS + " attempts. " +
                "Proposing evenly distributed population fractions instead. " +
                "Consider adjusting concentration-prior hyperparameters to increase prior mean.");
        return new TumorHeterogeneityState.PopulationFractions(Collections.nCopies(numPopulations, 1. / numPopulations));
    }

    /**
     * Calculate the per-segment ploidy-state modes
     */
    private static TumorHeterogeneityState.VariantProfile calculateVariantProfilePosteriorMode(final List<TumorHeterogeneityState.VariantProfile> variantProfiles) {
        final int numSegments = variantProfiles.get(0).numSegments();
        final List<Integer> ploidyStateIndicators = new ArrayList<>(numSegments);
        for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
            final int si = segmentIndex;
            final List<Integer> ploidyStateIndicatorSamples = variantProfiles.stream()
                    .map(vp -> vp.ploidyStateIndex(si)).collect(Collectors.toList());
            //get the mode of the samples; if there is more than one, return the first
            final int ploidyStateIndicatorMode = ploidyStateIndicatorSamples.stream()
                    .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()))
                    .entrySet().stream()
                    .collect(Collectors.groupingBy(Map.Entry::getValue, Collectors.mapping(Map.Entry::getKey, Collectors.toList())))
                    .entrySet().stream().max((o1, o2) -> o1.getKey().compareTo(o2.getKey())).map(Map.Entry::getValue)
                    .get().get(0);
            ploidyStateIndicators.add(ploidyStateIndicatorMode);
        }
        return new TumorHeterogeneityState.VariantProfile(
                new TumorHeterogeneityState.VariantProfile.PloidyStateIndicators(ploidyStateIndicators));
    }

    /**
     * Given population fractions, initialize corresponding population indicators for each cell;
     * the probability that a cell belongs to a population is given by the population fraction.
     */
    private static TumorHeterogeneityState.PopulationIndicators initializePopulationIndicators(final int numCells,
                                                                                               final TumorHeterogeneityState.PopulationFractions populationFractions,
                                                                                               final RandomGenerator rng) {
        final List<Integer> populationIndices = IntStream.range(0, populationFractions.size()).boxed().collect(Collectors.toList());
        final Function<Integer, Double> probabilityFunction = populationFractions::get;
        return new TumorHeterogeneityState.PopulationIndicators(IntStream.range(0, numCells).boxed()
                .map(p -> GATKProtectedMathUtils.randomSelect(populationIndices, probabilityFunction, rng))
                .collect(Collectors.toList()));
    }

    /**
     * Initialize variant profiles to normal.
     */
    private static TumorHeterogeneityState.VariantProfileCollection initializeNormalProfiles(final int numVariantPopulations,
                                                                                             final int numSegments,
                                                                                             final int normalPloidyStateIndex) {
        return new TumorHeterogeneityState.VariantProfileCollection(
                Collections.nCopies(numVariantPopulations, initializeNormalProfile(numSegments, normalPloidyStateIndex)));
    }

    /**
     * Initialize a variant profile to normal.
     */
    private static TumorHeterogeneityState.VariantProfile initializeNormalProfile(final int numSegments,
                                                                                  final int normalPloidyStateIndex) {
        final TumorHeterogeneityState.VariantProfile.PloidyStateIndicators ploidyStateIndicators =
                new TumorHeterogeneityState.VariantProfile.PloidyStateIndicators(Collections.nCopies(numSegments, normalPloidyStateIndex));
        return new TumorHeterogeneityState.VariantProfile(ploidyStateIndicators);
    }
}
