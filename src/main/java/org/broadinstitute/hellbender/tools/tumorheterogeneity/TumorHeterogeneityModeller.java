package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

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

    private final List<Boolean> doMetropolisStepSamples = new ArrayList<>();
    private final List<Double> concentrationSamples = new ArrayList<>();
    private final List<TumorHeterogeneityState.PopulationIndicators> populationIndicatorsSamples = new ArrayList<>();
    private final List<TumorHeterogeneityState.PopulationFractions> populationFractionsSamples = new ArrayList<>();
    private final List<TumorHeterogeneityState.VariantProfileCollection> variantProfileCollectionSamples = new ArrayList<>();

    public TumorHeterogeneityModeller(final TumorHeterogeneityData data,
                                      final TumorHeterogeneityPriorCollection priors,
                                      final int numPopulations,
                                      final int numCells,
                                      final RandomGenerator rng) {
        this(data, TumorHeterogeneityStateInitializationUtils.initializeState(priors.concentrationPriorAlpha() / priors.concentrationPriorBeta(), priors, data.numSegments(), numPopulations, numCells, rng), rng);
    }

    public TumorHeterogeneityModeller(final TumorHeterogeneityData data,
                                      final TumorHeterogeneityState initialState,
                                      final RandomGenerator rng) {
        Utils.nonNull(data);
        Utils.nonNull(initialState);
        Utils.nonNull(rng);

        this.data = data;
        this.priors = initialState.priors();

        final double concentrationSliceSamplingWidth = initialState.concentration();
        final int numPopulations = initialState.numPopulations();
        final int numCells = initialState.numCells();
        final int numVariantPopulations = numPopulations - 1;

        //define samplers
        final TumorHeterogeneitySamplers.DoMetropolisStepSampler metropolisStepSampler =
                new TumorHeterogeneitySamplers.DoMetropolisStepSampler();
        final TumorHeterogeneitySamplers.ConcentrationSampler concentrationSampler =
                new TumorHeterogeneitySamplers.ConcentrationSampler(CONCENTRATION_MIN, CONCENTRATION_MAX, concentrationSliceSamplingWidth);
        final TumorHeterogeneitySamplers.PopulationFractionsSampler populationFractionsSampler =
                new TumorHeterogeneitySamplers.PopulationFractionsSampler();
        final TumorHeterogeneitySamplers.PopulationIndicatorsSampler populationIndicatorsSampler =
                new TumorHeterogeneitySamplers.PopulationIndicatorsSampler(numCells, numPopulations);
        final TumorHeterogeneitySamplers.VariantProfileCollectionSampler variantProfileCollectionSampler =
                new TumorHeterogeneitySamplers.VariantProfileCollectionSampler(numVariantPopulations, priors.variantPloidyStatePrior());

        model = new ParameterizedModel.GibbsBuilder<>(initialState, data)
                .addParameterSampler(TumorHeterogeneityParameter.DO_METROPOLIS_STEP, metropolisStepSampler, Boolean.class)
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
        doMetropolisStepSamples.addAll(gibbsSampler.getSamples(TumorHeterogeneityParameter.DO_METROPOLIS_STEP,
                Boolean.class, numBurnIn));
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
     * Returns an unmodifiable view of the list of samples of the population-fractions posterior, represented as a list of
     * {@link TumorHeterogeneityState.PopulationFractions} objects.
     * @return  unmodifiable view of the list of samples of the population-fractions posterior
     */
    public List<TumorHeterogeneityState.PopulationFractions> getPopulationFractionsSamples() {
        return Collections.unmodifiableList(populationFractionsSamples);
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
            final boolean doMetropolisStep = doMetropolisStepSamples.get(sampleIndex);
            final double concentration = concentrationSamples.get(sampleIndex);
            final TumorHeterogeneityState.PopulationFractions populationFractions = populationFractionsSamples.get(sampleIndex);
            final TumorHeterogeneityState.PopulationIndicators populationIndicators = populationIndicatorsSamples.get(sampleIndex);
            final TumorHeterogeneityState.VariantProfileCollection variantProfileCollection = variantProfileCollectionSamples.get(sampleIndex);
            final TumorHeterogeneityState state =
                    new TumorHeterogeneityState(doMetropolisStep, concentration, populationFractions, populationIndicators, variantProfileCollection, priors);
            ploidySamples.add(state.calculatePopulationAndGenomicAveragedPloidy(data));
        }
        return ploidySamples;
    }

    public void output(final File outputFile,
                       final double credibleIntervalAlpha,
                       final JavaSparkContext ctx) {
        if (concentrationSamples.size() == 0) {
            throw new IllegalStateException("Cannot output modeller result before samples have been generated.");
        }
        try (final FileWriter writer = new FileWriter(outputFile)) {
            outputFile.createNewFile();

            //comments
            final PosteriorSummary ploidyPosteriorSummary =
                    PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(getPloidySamples(), credibleIntervalAlpha, ctx);
            final double ploidyPosteriorCenter = ploidyPosteriorSummary.getCenter();
            final double ploidyPosteriorWidth = calculatePosteriorWidth(ploidyPosteriorSummary);
            writer.write("#ploidy: " + ploidyPosteriorCenter + " " + ploidyPosteriorWidth);
            writer.write(System.getProperty("line.separator"));

            final PosteriorSummary concentrationPosteriorSummary =
                    PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(getConcentrationSamples(), credibleIntervalAlpha, ctx);
            final double concentrationPosteriorCenter = concentrationPosteriorSummary.getCenter();
            final double concentrationPosteriorWidth = calculatePosteriorWidth(concentrationPosteriorSummary);
            writer.write("#concentration: " + concentrationPosteriorCenter + " " + concentrationPosteriorWidth);
            writer.write(System.getProperty("line.separator"));

            final List<TumorHeterogeneityState.PopulationFractions> populationFractionsSamples = getPopulationFractionsSamples();
            final List<TumorHeterogeneityState.VariantProfileCollection> variantProfileCollectionSamples = getVariantProfileCollectionSamples();

            final int numPopulations = populationFractionsSamples.get(0).size();
            for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
                final int pi = populationIndex;
                final List<Double> populationFractionSamples = populationFractionsSamples.stream()
                        .map(s -> s.get(pi)).collect(Collectors.toList());
                final PosteriorSummary populationFractionPosteriorSummary =
                        PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(populationFractionSamples, credibleIntervalAlpha, ctx);
                final double populationFractionPosteriorCenter = populationFractionPosteriorSummary.getCenter();
                final double populationFractionPosteriorWidth = calculatePosteriorWidth(populationFractionPosteriorSummary);

                if (populationIndex != numPopulations - 1) {
                    final List<Double> variantSegmentFractionSamples = variantProfileCollectionSamples.stream()
                            .map(vpc -> vpc.get(pi).variantSegmentFraction()).collect(Collectors.toList());
                    final PosteriorSummary variantSegmentFractionPosteriorSummary =
                            PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(variantSegmentFractionSamples, credibleIntervalAlpha, ctx);
                    final double variantSegmentFractionPosteriorCenter = variantSegmentFractionPosteriorSummary.getCenter();
                    final double variantSegmentFractionPosteriorWidth = calculatePosteriorWidth(variantSegmentFractionPosteriorSummary);
                    writer.write("#population fraction " + populationIndex + ": " + populationFractionPosteriorCenter + " " + populationFractionPosteriorWidth + " " +
                            variantSegmentFractionPosteriorCenter + " " + variantSegmentFractionPosteriorWidth);
                    writer.write(System.getProperty("line.separator"));
                } else {
                    writer.write("#population fraction " + populationIndex + ": " + populationFractionPosteriorCenter + " " + populationFractionPosteriorWidth);
                    writer.write(System.getProperty("line.separator"));
                }
            }

            //column headers
            writer.write("POPULATION_INDEX\tSEGMENT_INDEX\tSEGMENT_INTERVAL\tIS_VARIANT_PROB\t");
            final List<PloidyState> variantPloidyStates = priors.variantPloidyStatePrior().ploidyStates();
            final int numVariantPloidyStates = variantPloidyStates.size();
            for (int variantPloidyStateIndex = 0; variantPloidyStateIndex < numVariantPloidyStates - 1; variantPloidyStateIndex++) {
                final PloidyState variantPloidyState = variantPloidyStates.get(variantPloidyStateIndex);
                writer.write(variantPloidyState.m() + "-" + variantPloidyState.n() + "\t");
            }
            final PloidyState variantPloidyState = variantPloidyStates.get(numVariantPloidyStates - 1);
            writer.write(variantPloidyState.m() + "-" + variantPloidyState.n());
            writer.write(System.getProperty("line.separator"));

            //rows
            for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
                final int pi = populationIndex;
                final double[] populationFractionSamples = populationFractionsSamples.stream().mapToDouble(s -> s.get(pi)).toArray();
                final double populationFractionPosteriorMean = new Mean().evaluate(populationFractionSamples);

                if (populationIndex != numPopulations - 1 && populationFractionPosteriorMean >= 0.01) {
                    for (int segmentIndex = 0; segmentIndex < data.numSegments(); segmentIndex++) {
                        final int si = segmentIndex;
                        final double[] isVariantSamples = variantProfileCollectionSamples.stream()
                                .mapToDouble(vpc -> vpc.get(pi).isVariant(si) ? 1. : 0.)
                                .toArray();
                        final double isVariantPosteriorMean = new Mean().evaluate(isVariantSamples);
                        writer.write(populationIndex + "\t" + segmentIndex + "\t" + data.segments().get(segmentIndex).getInterval() + "\t" + isVariantPosteriorMean + "\t");

                        for (int variantPloidyStateIndex = 0; variantPloidyStateIndex < numVariantPloidyStates; variantPloidyStateIndex++) {
                            final int vpsi = variantPloidyStateIndex;
                            final double[] isVariantPloidyStateSamples = variantProfileCollectionSamples.stream()
                                    .mapToDouble(vpc -> vpc.get(pi).variantPloidyStateIndex(si) == vpsi ? 1. : 0)
                                    .toArray();
                            final double variantPloidyStatePosteriorMean = new Mean().evaluate(isVariantPloidyStateSamples);
                            writer.write(String.format("%.3f", variantPloidyStatePosteriorMean));
                            if (variantPloidyStateIndex != numVariantPloidyStates - 1) {
                                writer.write("\t");
                            }
                        }
                        if (!(segmentIndex == data.numSegments() - 1 && populationIndex == numPopulations - 2)) {
                            writer.write(System.getProperty("line.separator"));
                        }
                    }
                }
            }
        } catch (final IOException e) {
            throw new GATKException("Error writing modeller result.");
        }
    }

    private static double calculatePosteriorWidth(final PosteriorSummary summary) {
        return (summary.getUpper() - summary.getLower()) / 2.;
    }
}
