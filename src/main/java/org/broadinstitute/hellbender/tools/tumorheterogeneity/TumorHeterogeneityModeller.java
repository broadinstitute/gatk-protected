package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.mcmc.*;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityModeller {
    private static final double EPSILON = 1E-10;

    private static final double CONCENTRATION_INITIAL = 0.1;
    private static final double CONCENTRATION_MIN = EPSILON;
    private static final double CONCENTRATION_MAX = 10.;
    private static final double CONCENTRATION_SLICE_SAMPLING_WIDTH = 0.005;
    private static final double CONCENTRATION_PRIOR_ALPHA = 1.;
    private static final double CONCENTRATION_PRIOR_BETA = 100.;

    private static final double VARIANCE_INITIAL = 0.1;
    private static final double VARIANCE_MIN = EPSILON;
    private static final double VARIANCE_MAX = 10.;
    private static final double VARIANCE_SLICE_SAMPLING_WIDTH = 0.01;

    private static final double MEAN_MIN = -10.;
    private static final double MEAN_MAX = 10.;
    private static final double MEAN_SLICE_SAMPLING_WIDTH = 0.01;

    private final ParameterizedModel<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> model;

    private final List<Double> concentrationSamples = new ArrayList<>();
    private final List<Double> varianceSamples = new ArrayList<>();
    private final List<TumorHeterogeneityState.PopulationFractions> populationFractionsSamples = new ArrayList<>();
    private final List<TumorHeterogeneityState.Means> meansSamples = new ArrayList<>();
    private final List<TumorHeterogeneityState.PopulationIndicators> populationIndicatorsSamples = new ArrayList<>();

    /**
     */
    public TumorHeterogeneityModeller(final List<ACNVModeledSegment> segments, final int numPopulations, final RandomGenerator rng) {
        final TumorHeterogeneityData data = new TumorHeterogeneityData(segments);

        final TumorHeterogeneityState.PopulationFractions initialPopulationFractions =
                new TumorHeterogeneityState.PopulationFractions(Collections.nCopies(numPopulations, 1. / numPopulations));
        final TumorHeterogeneityState.Means initialMeans =
                new TumorHeterogeneityState.Means(IntStream.range(0, numPopulations).boxed().map(i -> (i + 0.5) * (MEAN_MAX - MEAN_MIN) / numPopulations + MEAN_MIN).collect(Collectors.toList()));
        final List<Integer> populationIndices = IntStream.range(0, numPopulations).boxed().collect(Collectors.toList());
        final Function<Integer, Double> probabilityFunction = j -> 1. / numPopulations;
        final TumorHeterogeneityState.PopulationIndicators initialPopulationIndicators =
                new TumorHeterogeneityState.PopulationIndicators(IntStream.range(0, numPoints).boxed()
                        .map(p -> GATKProtectedMathUtils.randomSelect(populationIndices, probabilityFunction, rng)).collect(Collectors.toList()));

        final TumorHeterogeneityState initialState = new TumorHeterogeneityState(
                CONCENTRATION_INITIAL, VARIANCE_INITIAL, initialPopulationFractions, initialMeans, initialPopulationIndicators);

        final TumorHeterogeneitySamplers.ConcentrationSampler concentrationSampler =
                new TumorHeterogeneitySamplers.ConcentrationSampler(CONCENTRATION_MIN, CONCENTRATION_MAX, CONCENTRATION_SLICE_SAMPLING_WIDTH,
                        CONCENTRATION_PRIOR_ALPHA, CONCENTRATION_PRIOR_BETA);
        final TumorHeterogeneitySamplers.VarianceSampler varianceSampler =
                new TumorHeterogeneitySamplers.VarianceSampler(VARIANCE_MIN, VARIANCE_MAX, VARIANCE_SLICE_SAMPLING_WIDTH);
        final TumorHeterogeneitySamplers.PopulationFractionsSampler populationFractionsSampler =
                new TumorHeterogeneitySamplers.PopulationFractionsSampler(rng, numPopulations);
        final TumorHeterogeneitySamplers.MeansSampler meansSampler =
                new TumorHeterogeneitySamplers.MeansSampler(MEAN_MIN, MEAN_MAX, MEAN_SLICE_SAMPLING_WIDTH);
        final TumorHeterogeneitySamplers.PopulationIndicatorsSampler populationIndicatorsSampler =
                new TumorHeterogeneitySamplers.PopulationIndicatorsSampler();


        model = new ParameterizedModel.GibbsBuilder<>(initialState, data)
                .addParameterSampler(TumorHeterogeneityParameter.CONCENTRATION, concentrationSampler, Double.class)
                .addParameterSampler(TumorHeterogeneityParameter.MEANS, meansSampler, TumorHeterogeneityState.Means.class)
                .addParameterSampler(TumorHeterogeneityParameter.VARIANCE, varianceSampler, Double.class)
                .addParameterSampler(TumorHeterogeneityParameter.POPULATION_INDICATORS, populationIndicatorsSampler, TumorHeterogeneityState.PopulationIndicators.class)
                .addParameterSampler(TumorHeterogeneityParameter.POPULATION_FRACTIONS, populationFractionsSampler, TumorHeterogeneityState.PopulationFractions.class)
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
        //run MCMC
        final GibbsSampler<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> gibbsSampler
                = new GibbsSampler<>(numSamples, model);
        gibbsSampler.runMCMC();
        //update posterior samples
        concentrationSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.CONCENTRATION,
                Double.class, numBurnIn));
        varianceSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.VARIANCE,
                Double.class, numBurnIn));
        populationFractionsSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.POPULATION_FRACTIONS,
                TumorHeterogeneityState.PopulationFractions.class, numBurnIn));
        meansSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.MEANS,
                TumorHeterogeneityState.Means.class, numBurnIn));
        populationIndicatorsSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.POPULATION_INDICATORS,
                TumorHeterogeneityState.PopulationIndicators.class, numBurnIn));
    }

    /**
     * Returns an unmodifiable view of the list of samples of the concentration posterior.
     * @return  unmodifiable view of the list of samples of the concentration posterior
     */
    public List<Double> getConcentrationSamples() {
        return Collections.unmodifiableList(concentrationSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the variance posterior.
     * @return  unmodifiable view of the list of samples of the variance posterior
     */
    public List<Double> getVarianceSamples() {
        return Collections.unmodifiableList(varianceSamples);
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
     * Returns an unmodifiable view of the list of samples of the means posterior, represented as a list of
     * {@link TumorHeterogeneityState.Means} objects.
     * @return  unmodifiable view of the list of samples of the means posterior
     */
    public List<TumorHeterogeneityState.Means> getMeansSamples() {
        return Collections.unmodifiableList(meansSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the population-indicators posterior, represented as a list of
     * {@link TumorHeterogeneityState.PopulationIndicators} objects.
     * @return  unmodifiable view of the list of samples of the population-indicators posterior
     */
    public List<TumorHeterogeneityState.PopulationIndicators> getPopulationIndicatorsSamples() {
        return Collections.unmodifiableList(populationIndicatorsSamples);
    }

    /**
     * Returns a Map of {@link PosteriorSummary} elements summarizing the global parameters.
     * Should only be called after {@link TumorHeterogeneityModeller#fitMCMC(int, int)} has been called.
     * @param credibleIntervalAlpha credible-interval alpha, must be in (0, 1)
     * @param ctx                   {@link JavaSparkContext} used for mllib kernel density estimation
     * @return                      list of {@link PosteriorSummary} elements summarizing the global parameters
     */
    public Map<TumorHeterogeneityParameter, PosteriorSummary> getGlobalParameterPosteriorSummaries(final double credibleIntervalAlpha, final JavaSparkContext ctx) {
        final Map<TumorHeterogeneityParameter, PosteriorSummary> posteriorSummaries = new LinkedHashMap<>();
        posteriorSummaries.put(TumorHeterogeneityParameter.CONCENTRATION, PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(concentrationSamples, credibleIntervalAlpha, ctx));
        posteriorSummaries.put(TumorHeterogeneityParameter.VARIANCE, PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(varianceSamples, credibleIntervalAlpha, ctx));
        return posteriorSummaries;
    }
}
