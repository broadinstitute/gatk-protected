package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.*;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityModeller {
    private static final double EPSILON = 1E-12;

    protected static final double CONCENTRATION_MIN = EPSILON;
    protected static final double CONCENTRATION_MAX = 10.;

    private static final int NUM_SAMPLES_PER_LOG_ENTRY = 10;

    private final ParameterizedModel<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> model;
    private final TumorHeterogeneityData data;
    private final TumorHeterogeneityPriorCollection priors;

    private final List<Double> concentrationSamples = new ArrayList<>();
    private final List<TumorHeterogeneityState.PopulationIndicators> populationIndicatorsSamples = new ArrayList<>();
    private final List<TumorHeterogeneityState.PopulationFractions> populationFractionsSamples = new ArrayList<>();
    private final List<TumorHeterogeneityState.VariantProfileCollection> variantProfileCollectionSamples = new ArrayList<>();

    public TumorHeterogeneityModeller(final TumorHeterogeneityData data,
                                      final PloidyState normalPloidyState,
                                      final PloidyStatePrior variantPloidyStatePrior,
                                      final double concentrationPriorAlpha,
                                      final double concentrationPriorBeta,
                                      final double variantSegmentFractionPriorAlpha,
                                      final double variantSegmentFractionPriorBeta,
                                      final int numPopulations,
                                      final int numCells,
                                      final RandomGenerator rng) {
        Utils.nonNull(data);
        Utils.nonNull(normalPloidyState);
        Utils.nonNull(variantPloidyStatePrior);
        Utils.nonNull(rng);
        Utils.validateArg(Double.isNaN(variantPloidyStatePrior.logProbability(normalPloidyState)),
                "Variant-ploidy state prior should not be specified for normal ploidy state.");
        Utils.validateArg(concentrationPriorAlpha > 0, "Hyperparameter for concentration prior must be positive.");
        Utils.validateArg(concentrationPriorBeta > 0, "Hyperparameter for concentration prior must be positive.");
        Utils.validateArg(variantSegmentFractionPriorAlpha > 0, "Hyperparameter for variant-segment fraction must be positive.");
        Utils.validateArg(variantSegmentFractionPriorBeta > 0, "Hyperparameter for variant-segment fraction must be positive.");
        Utils.validateArg(numPopulations > 1, "Maximum number of populations must be strictly greater than 1.");
        Utils.validateArg(numCells > 0, "Number of auxiliary cells must be positive.");

        final int numVariantPopulations = numPopulations - 1;
        final int numVariantPloidyStates = variantPloidyStatePrior.numPloidyStates();

        //initialize population fractions to be evenly distributed
        final TumorHeterogeneityState.PopulationFractions initialPopulationFractions =
                new TumorHeterogeneityState.PopulationFractions(Collections.nCopies(numPopulations, 1. / numPopulations));
        //randomly initialize population indicators for each cell
        final TumorHeterogeneityState.PopulationIndicators initialPopulationIndicators =
                initializePopulationIndicators(numPopulations, numCells, rng);
        //initialize variant profiles
        final TumorHeterogeneityState.VariantProfileCollection initialVariantProfileCollection =
                initializeProfiles(numVariantPopulations, data.numSegments());

        //initialize TumorHeterogeneityState
        final double initialConcentration = concentrationPriorAlpha / concentrationPriorBeta;
        final double concentrationSliceSamplingWidth = initialConcentration;
        final TumorHeterogeneityPriorCollection priors = new TumorHeterogeneityPriorCollection(normalPloidyState, variantPloidyStatePrior,
                concentrationPriorAlpha, concentrationPriorBeta, variantSegmentFractionPriorAlpha, variantSegmentFractionPriorBeta);
        final TumorHeterogeneityState initialState = new TumorHeterogeneityState(
                initialConcentration, initialPopulationFractions, initialPopulationIndicators, initialVariantProfileCollection, priors);

        //define samplers
        final TumorHeterogeneitySamplers.ConcentrationSampler concentrationSampler =
                new TumorHeterogeneitySamplers.ConcentrationSampler(CONCENTRATION_MIN, CONCENTRATION_MAX, concentrationSliceSamplingWidth);
        final TumorHeterogeneitySamplers.PopulationFractionsSampler populationFractionsSampler =
                new TumorHeterogeneitySamplers.PopulationFractionsSampler();
        final TumorHeterogeneitySamplers.PopulationIndicatorsSampler populationIndicatorsSampler =
                new TumorHeterogeneitySamplers.PopulationIndicatorsSampler(numCells, numPopulations);
        final TumorHeterogeneitySamplers.VariantProfileCollectionSampler variantProfileCollectionSampler =
                new TumorHeterogeneitySamplers.VariantProfileCollectionSampler(numVariantPopulations, variantPloidyStatePrior);

        this.data = data;
        this.priors = priors;
        model = new ParameterizedModel.GibbsBuilder<>(initialState, data)
                .addParameterSampler(TumorHeterogeneityParameter.CONCENTRATION, concentrationSampler, Double.class)
                .addParameterSampler(TumorHeterogeneityParameter.POPULATION_INDICATORS, populationIndicatorsSampler, TumorHeterogeneityState.PopulationIndicators.class)
                .addParameterSampler(TumorHeterogeneityParameter.POPULATION_FRACTIONS, populationFractionsSampler, TumorHeterogeneityState.PopulationFractions.class)
                .addParameterSampler(TumorHeterogeneityParameter.VARIANT_PROFILES, variantProfileCollectionSampler, TumorHeterogeneityState.VariantProfileCollection.class)
                .build();
    }

    public TumorHeterogeneityModeller(final double initialConcentration,
                                      final TumorHeterogeneityState.PopulationFractions initialPopulationFractions,
                                      final TumorHeterogeneityState.PopulationIndicators initialPopulationIndicators,
                                      final TumorHeterogeneityState.VariantProfileCollection initialVariantProfileCollection,
                                      final TumorHeterogeneityData data,
                                      final PloidyState normalPloidyState,
                                      final PloidyStatePrior variantPloidyStatePrior,
                                      final double concentrationPriorAlpha,
                                      final double concentrationPriorBeta,
                                      final double variantSegmentFractionPriorAlpha,
                                      final double variantSegmentFractionPriorBeta,
                                      final int numPopulations,
                                      final int numCells,
                                      final RandomGenerator rng) {
        //TODO validate
        Utils.nonNull(normalPloidyState);
        Utils.nonNull(variantPloidyStatePrior);
        Utils.nonNull(rng);
        Utils.validateArg(Double.isNaN(variantPloidyStatePrior.logProbability(normalPloidyState)),
                "Variant-ploidy state prior should not be specified for normal ploidy state.");
        Utils.validateArg(concentrationPriorAlpha > 0, "Hyperparameter for concentration prior must be positive.");
        Utils.validateArg(concentrationPriorBeta > 0, "Hyperparameter for concentration prior must be positive.");
        Utils.validateArg(variantSegmentFractionPriorAlpha > 0, "Hyperparameter for variant-segment fraction must be positive.");
        Utils.validateArg(variantSegmentFractionPriorBeta > 0, "Hyperparameter for variant-segment fraction must be positive.");
        Utils.validateArg(numPopulations > 1, "Maximum number of populations must be strictly greater than 1.");
        Utils.validateArg(numCells > 0, "Number of auxiliary cells must be positive.");

        final int numVariantPopulations = numPopulations - 1;

        final double concentrationSliceSamplingWidth = initialConcentration;
        //define samplers
        final TumorHeterogeneitySamplers.ConcentrationSampler concentrationSampler =
                new TumorHeterogeneitySamplers.ConcentrationSampler(CONCENTRATION_MIN, CONCENTRATION_MAX, concentrationSliceSamplingWidth);
        final TumorHeterogeneitySamplers.PopulationFractionsSampler populationFractionsSampler =
                new TumorHeterogeneitySamplers.PopulationFractionsSampler();
        final TumorHeterogeneitySamplers.PopulationIndicatorsSampler populationIndicatorsSampler =
                new TumorHeterogeneitySamplers.PopulationIndicatorsSampler(numCells, numPopulations);
        final TumorHeterogeneitySamplers.VariantProfileCollectionSampler variantProfileCollectionSampler =
                new TumorHeterogeneitySamplers.VariantProfileCollectionSampler(numVariantPopulations, variantPloidyStatePrior);

        final TumorHeterogeneityPriorCollection priors = new TumorHeterogeneityPriorCollection(normalPloidyState, variantPloidyStatePrior,
                concentrationPriorAlpha, concentrationPriorBeta, variantSegmentFractionPriorAlpha, variantSegmentFractionPriorBeta);
        final TumorHeterogeneityState initialState = new TumorHeterogeneityState(
                initialConcentration, initialPopulationFractions, initialPopulationIndicators, initialVariantProfileCollection, priors);

        this.data = data;
        this.priors = priors;
        model = new ParameterizedModel.GibbsBuilder<>(initialState, data)
                .addParameterSampler(TumorHeterogeneityParameter.CONCENTRATION, concentrationSampler, Double.class)
                .addParameterSampler(TumorHeterogeneityParameter.POPULATION_INDICATORS, populationIndicatorsSampler, TumorHeterogeneityState.PopulationIndicators.class)
                .addParameterSampler(TumorHeterogeneityParameter.POPULATION_FRACTIONS, populationFractionsSampler, TumorHeterogeneityState.PopulationFractions.class)
                .addParameterSampler(TumorHeterogeneityParameter.VARIANT_PROFILES, variantProfileCollectionSampler, TumorHeterogeneityState.VariantProfileCollection.class)
                .build();
    }

    /**
     * Adds {@code numSamples - numBurnIn} Markov-Chain Monte-Carlo samples of the parameter posteriors (generated using
     * Gibbs sampling) to the collections held internally.  The current {@link TumorHeterogeneityState} held internally is used
     * to initialize the Markov Chain.
     * @param numSamples    total number of samples per posterior
     * @param numBurnIn     number of burn-in samples to discard
     */
    public void fitMCMC(final int numSamples, final int numBurnIn) {
        Utils.validateArg(numSamples > 0, "Total number of samples must be positive.");
        Utils.validateArg(0 <= numBurnIn && numBurnIn < numSamples,
                "Number of burn-in samples to discard must be non-negative and strictly less than total number of samples.");
        //run MCMC
        final GibbsSampler<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> gibbsSampler
                = new GibbsSampler<>(numSamples, model);
        gibbsSampler.setNumSamplesPerLogEntry(NUM_SAMPLES_PER_LOG_ENTRY);
        gibbsSampler.runMCMC();
        //update posterior samples
        concentrationSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.CONCENTRATION,
                Double.class, numBurnIn));
        populationIndicatorsSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.POPULATION_INDICATORS,
                TumorHeterogeneityState.PopulationIndicators.class, numBurnIn));
        populationFractionsSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.POPULATION_FRACTIONS,
                TumorHeterogeneityState.PopulationFractions.class, numBurnIn));
        variantProfileCollectionSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.VARIANT_PROFILES,
                TumorHeterogeneityState.VariantProfileCollection.class, numBurnIn));
    }

    /**
     * Returns an unmodifiable view of the list of samples of the concentration posterior.
     * @return  unmodifiable view of the list of samples of the concentration posterior
     */
    public List<Double> getConcentrationSamples() {
        return Collections.unmodifiableList(concentrationSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the population-indicators posterior, represented as a list of
     * {@link TumorHeterogeneityState.PopulationFractions} objects.
     * @return  unmodifiable view of the list of samples of the population-indicators posterior
     */
    public List<TumorHeterogeneityState.PopulationIndicators> getPopulationIndicatorsSamples() {
        return Collections.unmodifiableList(populationIndicatorsSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the population-fractions posterior, represented as a list of
     * {@link TumorHeterogeneityState.PopulationFractions} objects.
     * @return  unmodifiable view of the list of samples of the population-fractions posterior
     */
    public List<TumorHeterogeneityState.PopulationFractions> getPopulationFractionsSamples() {
        return Collections.unmodifiableList(populationFractionsSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the variant-profile-collection posterior, represented as a list of
     * {@link TumorHeterogeneityState.PopulationFractions} objects.
     * @return  unmodifiable view of the list of samples of the variant-profile-collection posterior
     */
    public List<TumorHeterogeneityState.VariantProfileCollection> getVariantProfileCollectionSamples() {
        return Collections.unmodifiableList(variantProfileCollectionSamples);
    }

    public List<Double> getPloidySamples() {
        final int numSamples = concentrationSamples.size();
        final List<Double> ploidySamples = new ArrayList<>(numSamples);
        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
            final double concentration = concentrationSamples.get(sampleIndex);
            final TumorHeterogeneityState.PopulationFractions populationFractions = populationFractionsSamples.get(sampleIndex);
            final TumorHeterogeneityState.PopulationIndicators populationIndicators = populationIndicatorsSamples.get(sampleIndex);
            final TumorHeterogeneityState.VariantProfileCollection variantProfileCollection = variantProfileCollectionSamples.get(sampleIndex);
            final TumorHeterogeneityState state = new TumorHeterogeneityState(concentration, populationFractions, populationIndicators, variantProfileCollection, priors);
            ploidySamples.add(state.calculatePopulationAndGenomicAveragedPloidy(data));
        }
        return ploidySamples;
    }

    /**
     * Returns a Map of {@link PosteriorSummary} elements summarizing the global parameters.
     * Should only be called after {@link TumorHeterogeneityModeller#fitMCMC(int, int)} has been called.
     * @param credibleIntervalAlpha credible-interval alpha, must be in (0, 1)
     * @param ctx                   {@link JavaSparkContext} used for mllib kernel density estimation
     * @return                      list of {@link PosteriorSummary} elements summarizing the global parameters
     */
    public Map<TumorHeterogeneityParameter, PosteriorSummary> getGlobalParameterPosteriorSummaries(final double credibleIntervalAlpha,
                                                                                                   final JavaSparkContext ctx) {
        Utils.validateArg(0. <= credibleIntervalAlpha && credibleIntervalAlpha <= 1., "Credible-interval alpha must be in [0, 1].");
        Utils.nonNull(ctx);
        final Map<TumorHeterogeneityParameter, PosteriorSummary> posteriorSummaries = new LinkedHashMap<>();
        posteriorSummaries.put(TumorHeterogeneityParameter.CONCENTRATION, PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(concentrationSamples, credibleIntervalAlpha, ctx));
        return posteriorSummaries;
    }

    static TumorHeterogeneityState.PopulationIndicators initializePopulationIndicators(final int numPopulations,
                                                                                       final int numCells,
                                                                                       final RandomGenerator rng) {
        final List<Integer> populationIndices = IntStream.range(0, numPopulations).boxed().collect(Collectors.toList());
        final Function<Integer, Double> probabilityFunction = i -> 1. / numPopulations;
        return new TumorHeterogeneityState.PopulationIndicators(IntStream.range(0, numCells).boxed()
                .map(p -> GATKProtectedMathUtils.randomSelect(populationIndices, probabilityFunction, rng))
                .collect(Collectors.toList()));
    }

    static TumorHeterogeneityState.VariantProfileCollection initializeProfiles(final int numVariantPopulations,
                                                                               final int numSegments) {
        return new TumorHeterogeneityState.VariantProfileCollection(Collections.nCopies(numVariantPopulations, initializeProfile(numSegments)));
    }

    static TumorHeterogeneityState.VariantProfile initializeProfile(final int numSegments) {
        final double variantSegmentFraction = 0.;
        final TumorHeterogeneityState.VariantProfile.VariantIndicators variantIndicators =
                new TumorHeterogeneityState.VariantProfile.VariantIndicators(Collections.nCopies(numSegments, false));
        final TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators variantPloidyStateIndicators =
                new TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators(Collections.nCopies(numSegments, 0));
        return new TumorHeterogeneityState.VariantProfile(variantSegmentFraction, variantIndicators, variantPloidyStateIndicators);
    }
}
