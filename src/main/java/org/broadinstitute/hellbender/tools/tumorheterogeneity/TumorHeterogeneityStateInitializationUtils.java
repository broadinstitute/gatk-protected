package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
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
    static TumorHeterogeneityState initializeState(final TumorHeterogeneityPriorCollection priors,
                                                   final int numSegments,
                                                   final int numPopulations) {
        final double concentration = priors.concentrationPriorAlpha() / priors.concentrationPriorBeta();
        final double copyRatioNoiseFloor = priors.copyRatioNoiseFloorPriorAlpha() / priors.copyRatioNoiseFloorPriorBeta();
        final double copyRatioNoiseFactor = 1. + priors.copyRatioNoiseFactorPriorAlpha() / priors.copyRatioNoiseFactorPriorBeta();
        final double minorAlleleFractionNoiseFactor = 1. + priors.minorAlleleFractionNoiseFactorPriorAlpha() / priors.minorAlleleFractionNoiseFactorPriorBeta();
        //initialize population fractions to be evenly distributed
        final PopulationMixture.PopulationFractions populationFractions =
                new PopulationMixture.PopulationFractions(Collections.nCopies(numPopulations, 1. / numPopulations));
        //initialize variant profiles to normal
        final int numVariantPopulations = numPopulations - 1;
        final PopulationMixture.VariantProfileCollection variantProfileCollection =
                initializeNormalProfiles(numVariantPopulations, numSegments, priors.normalPloidyState());
        final PopulationMixture populationMixture = new PopulationMixture(populationFractions, variantProfileCollection, priors.normalPloidyState());
        return new TumorHeterogeneityState(
                concentration, copyRatioNoiseFloor, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationMixture, priors);
    }

    /**
     * Initialize state from a modeller run that assumes one clonal population and one normal population.
     */
    static TumorHeterogeneityState initializeStateFromClonalResult(final TumorHeterogeneityPriorCollection priors,
                                                                   final TumorHeterogeneityModeller clonalModeller,
                                                                   final int maxNumPopulations) {
        //double check some of the parameters
        Utils.validateArg(clonalModeller.getConcentrationSamples().size() > 0,
                "Clonal modeller must have a non-zero number of samples.");
        Utils.validateArg(clonalModeller.getPopulationFractionsSamples().get(0).size() == 2,
                "Clonal modeller must have two populations (clone + normal).");

        //initialize global parameters to prior mean
        final double concentration = priors.concentrationPriorAlpha() / priors.concentrationPriorBeta();
        final double copyRatioNoiseFloor = priors.copyRatioNoiseFloorPriorAlpha() / priors.copyRatioNoiseFloorPriorBeta();
        final double copyRatioNoiseFactor = 1. + priors.copyRatioNoiseFactorPriorAlpha() / priors.copyRatioNoiseFactorPriorBeta();
        final double minorAlleleFractionNoiseFactor = 1. + priors.minorAlleleFractionNoiseFactorPriorAlpha() / priors.minorAlleleFractionNoiseFactorPriorBeta();

        //initialize normal fraction to posterior mean of clonal result
        final double[] normalFractionSamples = clonalModeller.getPopulationFractionsSamples().stream()
                .mapToDouble(pfs -> pfs.get(NUM_POPULATIONS_CLONAL - 1)).toArray();
        final double normalFraction = new Mean().evaluate(normalFractionSamples);

        //initialize one variant profile to posterior mode of clonal result
        final List<PopulationMixture.VariantProfile> clonalProfileSamples = clonalModeller.getVariantProfileCollectionSamples().stream()
                .map(vpc -> vpc.get(0)).collect(Collectors.toList());
        final PopulationMixture.VariantProfile initialClonalProfile = calculateVariantProfilePosteriorMode(clonalProfileSamples);

        //build new population fractions and variant-profile collection with additional variant populations
        //split clonal fraction evenly among new populations and initialize with clonal profile
        final List<Double> initialFractions = new ArrayList<>();
        final List<PopulationMixture.VariantProfile> initialVariantProfiles = new ArrayList<>();
        //add clonal population
        initialFractions.add((1. - normalFraction) / (maxNumPopulations - 1));
        initialVariantProfiles.add(initialClonalProfile);
        //initialize additional variant profiles
        for (int i = 0; i < maxNumPopulations - NUM_POPULATIONS_CLONAL; i++) {
            initialFractions.add((1. - normalFraction) / (maxNumPopulations - 1));
            initialVariantProfiles.add(1, new PopulationMixture.VariantProfile(initialClonalProfile));
        }
        //add normal population fraction
        initialFractions.add(normalFraction);
        final PopulationMixture populationMixture = new PopulationMixture(initialFractions, initialVariantProfiles, priors.normalPloidyState());

        return new TumorHeterogeneityState(
                concentration, copyRatioNoiseFloor, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationMixture, priors);
    }

    /**
     * Sample population fractions from a symmetric Dirichlet distribution.
     */
    private static PopulationMixture.PopulationFractions proposePopulationFractions(final int numPopulations,
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
                logger.debug("Proposing population fractions using prior.");
                logger.debug("Proposed population fractions: " + populationFractions);
                return new PopulationMixture.PopulationFractions(populationFractions);
            }
        }
        //if non-zero population fractions were not sampled, emit warning and return evenly distributed fractions
        logger.warn("Failed to sample non-zero population fractions within " + MAX_POPULATION_FRACTIONS_ITERATIONS + " attempts. " +
                "Proposing evenly distributed population fractions instead. " +
                "Consider adjusting concentration-prior hyperparameters to increase prior mean.");
        return new PopulationMixture.PopulationFractions(Collections.nCopies(numPopulations, 1. / numPopulations));
    }

    private static PopulationMixture.PopulationFractions proposePopulationFractions(final PopulationMixture.PopulationFractions populationFractions,
                                                                                    final double concentration,
                                                                                    final double proposalWidthFactor,
                                                                                    final RandomGenerator rng) {
        final int numPopulations = populationFractions.size();
        //sampling from Dirichlet(alpha_vec) is equivalent to sampling from individual Gamma(alpha_vec_i, 1) distributions and normalizing
        final double[] unnormalizedPopulationFractions = IntStream.range(0, numPopulations).boxed().mapToDouble(pi -> new GammaDistribution(rng, concentration + proposalWidthFactor * populationFractions.get(pi), 1.).sample()).toArray();
        final List<Double> proposedPopulationFractions = Doubles.asList(MathUtils.normalizeFromRealSpace(unnormalizedPopulationFractions));
        logger.debug("Proposing population fractions using posterior.");
        logger.debug("Proposed population fractions: " + proposedPopulationFractions);
        return new PopulationMixture.PopulationFractions(proposedPopulationFractions);
    }

    /**
     * Calculate the per-segment ploidy-state modes
     */
    private static PopulationMixture.VariantProfile calculateVariantProfilePosteriorMode(final List<PopulationMixture.VariantProfile> variantProfiles) {
        final int numSegments = variantProfiles.get(0).numSegments();
        final List<PloidyState> ploidyStates = new ArrayList<>(numSegments);
        for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
            final int si = segmentIndex;
            final List<PloidyState> ploidyStateSamples = variantProfiles.stream()
                    .map(vp -> vp.ploidyState(si)).collect(Collectors.toList());
            //get the mode of the samples; if there is more than one, return the first
            final PloidyState ploidyStateMode = ploidyStateSamples.stream()
                    .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()))
                    .entrySet().stream()
                    .collect(Collectors.groupingBy(Map.Entry::getValue, Collectors.mapping(Map.Entry::getKey, Collectors.toList())))
                    .entrySet().stream().max((o1, o2) -> o1.getKey().compareTo(o2.getKey())).map(Map.Entry::getValue)
                    .get().get(0);
            ploidyStates.add(ploidyStateMode);
        }
        return new PopulationMixture.VariantProfile(ploidyStates);
    }

    /**
     * Initialize variant profiles to normal.
     */
    private static PopulationMixture.VariantProfileCollection initializeNormalProfiles(final int numVariantPopulations,
                                                                                       final int numSegments,
                                                                                       final PloidyState normalPloidyState) {
        return new PopulationMixture.VariantProfileCollection(
                Collections.nCopies(numVariantPopulations, initializeNormalProfile(numSegments, normalPloidyState)));
    }

    /**
     * Initialize a variant profile to normal.
     */
    private static PopulationMixture.VariantProfile initializeNormalProfile(final int numSegments,
                                                                            final PloidyState normalPloidyState) {
        final PopulationMixture.VariantProfile ploidyStateIndicators =
                new PopulationMixture.VariantProfile(Collections.nCopies(numSegments, normalPloidyState));
        return new PopulationMixture.VariantProfile(ploidyStateIndicators);
    }
}
