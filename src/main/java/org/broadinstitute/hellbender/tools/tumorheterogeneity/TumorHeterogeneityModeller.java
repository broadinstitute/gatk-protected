package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Sets;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.ModelSampler;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedModel;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.WalkerPosition;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.PosteriorSummaryWriter;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.PosteriorSummary;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.PosteriorSummaryUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityModeller {
    private static final Logger logger = LogManager.getLogger(TumorHeterogeneityModeller.class);

    private static final int NUM_SAMPLES_PER_LOG_ENTRY = 100;
    private static final double INITIAL_BALL_SIZE = 1.;
    private static final int MAX_NUM_PROPOSALS_INITIAL_WALKER_BALL = 25;

    private final ParameterizedModel.EnsembleBuilder<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> builder;
    private final ParameterizedModel<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> model;
    private final TumorHeterogeneityData data;
    private final int numWalkers;

    private final List<Double> concentrationSamples = new ArrayList<>();
    private final List<Double> copyRatioNoiseConstantSamples = new ArrayList<>();
    private final List<Double> copyRatioNoiseFactorSamples = new ArrayList<>();
    private final List<Double> minorAlleleFractionNoiseFactorSamples = new ArrayList<>();
    private final List<Double> ploidySamples = new ArrayList<>();
    private final List<PopulationMixture> populationMixtureSamples = new ArrayList<>();

    public TumorHeterogeneityModeller(final TumorHeterogeneityData data,
                                      final int numPopulations,
                                      final int numWalkers,
                                      final RandomGenerator rng) {
        this(data, TumorHeterogeneityState.initializeState(data.priors(), data.numSegments(), numPopulations), numWalkers, rng);
    }

    public TumorHeterogeneityModeller(final TumorHeterogeneityData data,
                                      final TumorHeterogeneityState initialState,
                                      final int numWalkers,
                                      final RandomGenerator rng) {
        Utils.nonNull(data);
        Utils.nonNull(initialState);
        Utils.validateArg(numWalkers >= 2, "Number of walkers must be greater than or equal to two.");
        Utils.nonNull(rng);

        this.data = data;
        this.numWalkers = numWalkers;

        //define log-target function
        final Function<TumorHeterogeneityState, Double> logTargetTumorHeterogeneity = state ->
                TumorHeterogeneityUtils.calculateLogJacobianFactor(state, data) + TumorHeterogeneityUtils.calculateLogPosterior(state, data);

        //define walker transformation
        final int numPopulations = initialState.populationMixture().numPopulations();
        final List<PloidyState> ploidyStates = data.priors().ploidyStatePrior().ploidyStates();
        final Set<Integer> totalCopyNumberStates = ploidyStates.stream().map(PloidyState::total).collect(Collectors.toSet());
        final List<List<Integer>> totalCopyNumberProductStates =
                new ArrayList<>(Sets.cartesianProduct(Collections.nCopies(numPopulations, totalCopyNumberStates)));
        final Map<Integer, Set<PloidyState>> ploidyStateSetsMap = new HashMap<>();
        for (final int totalCopyNumber : totalCopyNumberStates) {
            final Set<PloidyState> ploidyStateSet = ploidyStates.stream().filter(ps -> ps.total() == totalCopyNumber).collect(Collectors.toSet());
            ploidyStateSetsMap.put(totalCopyNumber, ploidyStateSet);
        }
        final Function<WalkerPosition, TumorHeterogeneityState> transformWalkerPositionToState = walkerPosition ->
                TumorHeterogeneityUtils.transformWalkerPositionToState(walkerPosition, rng, data, totalCopyNumberProductStates, ploidyStateSetsMap);

        //initialize walker positions as a ball around initialState
        final List<WalkerPosition> initialWalkerPositions = initializeWalkerBall(rng, initialState, logTargetTumorHeterogeneity, transformWalkerPositionToState);

        builder = new ParameterizedModel.EnsembleBuilder<>(initialWalkerPositions, data, transformWalkerPositionToState, logTargetTumorHeterogeneity);
        model = builder.build();
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
        final ModelSampler<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> modelSampler
                = new ModelSampler<>(numWalkers * numSamples, model);
        modelSampler.setNumSamplesPerLogEntry(NUM_SAMPLES_PER_LOG_ENTRY);
        modelSampler.runMCMC();
        //update posterior samples
        concentrationSamples.addAll(modelSampler.getSamples(TumorHeterogeneityParameter.CONCENTRATION,
                Double.class, numWalkers * numBurnIn));
        copyRatioNoiseConstantSamples.addAll(modelSampler.getSamples(TumorHeterogeneityParameter.COPY_RATIO_NOISE_CONSTANT,
                Double.class, numWalkers * numBurnIn));
        copyRatioNoiseFactorSamples.addAll(modelSampler.getSamples(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FACTOR,
                Double.class, numWalkers * numBurnIn));
        minorAlleleFractionNoiseFactorSamples.addAll(modelSampler.getSamples(TumorHeterogeneityParameter.MINOR_ALLELE_FRACTION_NOISE_FACTOR,
                Double.class, numWalkers * numBurnIn));
        ploidySamples.addAll(modelSampler.getSamples(TumorHeterogeneityParameter.PLOIDY,
                Double.class, numWalkers * numBurnIn));
        populationMixtureSamples.addAll(modelSampler.getSamples(TumorHeterogeneityParameter.POPULATION_MIXTURE,
                PopulationMixture.class, numWalkers * numBurnIn));
    }

    public List<Double> getConcentrationSamples() {
        return Collections.unmodifiableList(concentrationSamples);
    }

    public List<Double> getCopyRatioNoiseConstantSamples() {
        return Collections.unmodifiableList(copyRatioNoiseConstantSamples);
    }

    public List<Double> getCopyRatioNoiseFactorSamples() {
        return Collections.unmodifiableList(copyRatioNoiseFactorSamples);
    }

    public List<Double> getMinorAlleleFractionNoiseFactorSamples() {
        return Collections.unmodifiableList(minorAlleleFractionNoiseFactorSamples);
    }

    public List<Double> getPloidySamples() {
        return Collections.unmodifiableList(ploidySamples);
    }

    public List<PopulationMixture.PopulationFractions> getPopulationFractionsSamples() {
        final PloidyState normalPloidyState = data.priors().normalPloidyState();
        return Collections.unmodifiableList(populationMixtureSamples.stream()
                .map(pm -> pm.collapseNormalPopulations(normalPloidyState))
                .map(PopulationMixture::populationFractions).collect(Collectors.toList()));
    }

    public List<PopulationMixture.VariantProfileCollection> getVariantProfileCollectionSamples() {
        final PloidyState normalPloidyState = data.priors().normalPloidyState();
        return Collections.unmodifiableList(populationMixtureSamples.stream()
                .map(pm -> pm.collapseNormalPopulations(normalPloidyState))
                .map(PopulationMixture::variantProfileCollection).collect(Collectors.toList()));
    }

    public TumorHeterogeneityData getData() {
        return data;
    }

    public TumorHeterogeneityState getMaxAPosterioriState() {
        return builder.getMaxLogTargetState();
    }

    public List<Integer> identifySamplesAtMode(final double purityModeBinSize,
                                               final double ploidyModeBinSize) {
        Utils.validateArg(0. <= purityModeBinSize && purityModeBinSize <= 1.,
                "Invalid purity bin size for determining mode.");
        Utils.validateArg(0. <= ploidyModeBinSize && ploidyModeBinSize <= data.priors().ploidyStatePrior().maxCopyNumber(),
                "Invalid ploidy bin size for determining mode.");
        if (getConcentrationSamples().size() == 0) {
            throw new IllegalStateException("Cannot output modeller result before samples have been generated.");
        }

        //get maximum a posteriori state and collapse normal populations
        final TumorHeterogeneityState maxLogPosteriorState = new TumorHeterogeneityState(
                getMaxAPosterioriState().concentration(),
                getMaxAPosterioriState().copyRatioNoiseConstant(),
                getMaxAPosterioriState().copyRatioNoiseFactor(),
                getMaxAPosterioriState().minorAlleleFractionNoiseFactor(),
                getMaxAPosterioriState().initialPloidy(),
                getMaxAPosterioriState().ploidy(),
                getMaxAPosterioriState().populationMixture().collapseNormalPopulations(data.priors().normalPloidyState()));

        final double purityMode = maxLogPosteriorState.populationMixture().populationFractions().tumorFraction();
        final double purityModeBinMin = Math.max(0., purityMode - purityModeBinSize / 2);
        final double purityModeBinMax = Math.min(1., purityMode + purityModeBinSize / 2);
        logger.info("Mode purity bin: [" + purityModeBinMin + ", " + purityModeBinMax + ")");

        final double ploidyMode = maxLogPosteriorState.ploidy();
        final double ploidyModeBinMin = Math.max(TumorHeterogeneityUtils.PLOIDY_MIN, ploidyMode - ploidyModeBinSize / 2);
        final double ploidyModeBinMax = Math.min(data.priors().ploidyStatePrior().maxCopyNumber(), ploidyMode + ploidyModeBinSize / 2);
        logger.info("Mode ploidy bin: [" + ploidyModeBinMin + ", " + ploidyModeBinMax + ")");

        final List<PopulationMixture.PopulationFractions> populationFractionsSamples = getPopulationFractionsSamples();
        final List<Integer> sampleIndices = IntStream.range(0, getPloidySamples().size()).boxed()
                .filter(i -> purityModeBinMin <= populationFractionsSamples.get(i).tumorFraction()
                        && populationFractionsSamples.get(i).tumorFraction() < purityModeBinMax
                        && ploidyModeBinMin <= getPloidySamples().get(i)
                        && getPloidySamples().get(i) < ploidyModeBinMax)
                .collect(Collectors.toList());
        logger.info("Number of samples in mode bin: " + sampleIndices.size());
        return sampleIndices;
    }

    private List<WalkerPosition> initializeWalkerBall(final RandomGenerator rng,
                                                      final TumorHeterogeneityState initialState,
                                                      final Function<TumorHeterogeneityState, Double> logTargetTumorHeterogeneity,
                                                      final Function<WalkerPosition, TumorHeterogeneityState> transformWalkerPositionToState) {
        //number of walker dimensions = 1 concentration parameter + 2 noise parameters + (numPopulations - 1) simplex parameters + 1 ploidy parameter
        final int numDimensions = TumorHeterogeneityUtils.NUM_GLOBAL_PARAMETERS + initialState.populationMixture().numPopulations() - 1;
        final NormalDistribution ballGaussian = new NormalDistribution(rng, 0., INITIAL_BALL_SIZE);
        final WalkerPosition walkerPositionOfInitialState = TumorHeterogeneityUtils.transformStateToWalkerPosition(initialState, data);
        final List<WalkerPosition> initialWalkerPositions = new ArrayList<>(numWalkers);
        for (int walkerIndex = 0; walkerIndex < numWalkers; walkerIndex++) {
            WalkerPosition initialWalkerPosition = walkerPositionOfInitialState;
            for (int proposalIndex = 0; proposalIndex < MAX_NUM_PROPOSALS_INITIAL_WALKER_BALL; proposalIndex++) {
                final WalkerPosition proposedWalkerPosition = new WalkerPosition(
                        IntStream.range(0, numDimensions).boxed()
                                .map(dimensionIndex -> walkerPositionOfInitialState.get(dimensionIndex) + ballGaussian.sample())
                                .collect(Collectors.toList()));
                if (Double.isFinite(logTargetTumorHeterogeneity.apply(transformWalkerPositionToState.apply(proposedWalkerPosition)))) {
                    initialWalkerPosition = proposedWalkerPosition;
                    break;
                }
            }
            initialWalkerPositions.add(initialWalkerPosition);
        }
        return initialWalkerPositions;
    }
}
