package org.broadinstitute.hellbender.utils.mcmc;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.WalkerPosition;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents a parameterized model.  The parameterized state of the model is represented by an
 * {@link ParameterizedState}, while the data is represented by an {@link DataCollection}.
 * Allows for sampling via Gibbs sampling (if {@link ParameterizedModel.GibbsBuilder} is used in construction)
 * or affine-invariant ensemble sampling (if {@link ParameterizedModel.EnsembleBuilder} is used in construction).
 * See GibbsSamplerSingleGaussianUnitTest and GibbsSamplerCopyRatioUnitTest for examples of Gibbs sampling.
 * @param <S1>  type of the ParameterizedState
 * @param <T1>  type of the DataCollection
 */
public final class ParameterizedModel<V1 extends Enum<V1> & ParameterEnum, S1 extends ParameterizedState<V1>, T1 extends DataCollection> {
    //enums for specifying method of updating ParameterizedState
    //updateMethod should be implemented accordingly within Builders and constructors corresponding to each update method
    private enum UpdateMethod {
        GIBBS, ENSEMBLE
    }

    private static final Logger logger = LogManager.getLogger(ParameterizedModel.class);

    private final S1 state;
    private final T1 dataCollection;
    private final Consumer<RandomGenerator> updateState;
    private final UpdateMethod updateMethod;

    /**
     * Builder for constructing a {@link ParameterizedModel} to be Gibbs sampled using {@link ModelSampler}.
     * Given an initial instance "initialState" of a ConcreteParameterizedState (which extends
     * {@link ParameterizedState}) and an instance "dataset" of a ConcreteDataCollection (which extends
     * {@link DataCollection}), as well as i = 1,...,N {@link ParameterSampler} objects SAMPLER_i that return samples
     * of type TYPE_i, a {@link ParameterizedModel} model can be constructed using the Builder pattern as:
     *
     *  ParameterizedModel<ConcreteParameterizedState, ConcreteDataCollection> model =
     *      new ParameterizedModel.GibbsBuilder<>(initialState, dataset, ConcreteParameterizedState.class)
     *                            .addParameterSampler(SAMPLER_1, TYPE_1.class)
     *                            .
     *                            .
     *                            .
     *                            .addParameterSampler(SAMPLER_N, TYPE_N.class)
     *                            .build()
     *
     * See GibbsSamplerSingleGaussianUnitTest and GibbsSamplerCopyRatioUnitTest for examples of use.
     * @param <V2>  type of the ParameterEnum
     * @param <S2>  type of the ParameterizedState
     * @param <T2>  type of the DataCollection
     */
    public static final class GibbsBuilder<V2 extends Enum<V2> & ParameterEnum, S2 extends ParameterizedState<V2>, T2 extends DataCollection> {
        private final S2 state;
        private final T2 dataCollection;
        private final Map<V2, ParameterSampler<?, V2, S2, T2>> samplerMap = new HashMap<>();

        /**
         * Constructor for {@link ParameterizedModel.GibbsBuilder}.
         * @param state             ParameterizedState held by the model
         * @param dataCollection    DataCollection used by the model
         */
        public GibbsBuilder(final S2 state, final T2 dataCollection) {
            Utils.nonNull(state);
            Utils.nonNull(dataCollection);
            this.state = state;
            this.dataCollection = dataCollection;
        }

        /**
         * Adds a {@link ParameterSampler} to the collection of parameter samplers using {@link ParameterizedModel.GibbsBuilder}.
         * @param parameterName         name of parameter to sample
         * @param parameterSampler      ParameterSampler that returns random samples of the parameter
         * @param parameterValueClass   class of the parameter value to sample
         * @param <U>                   type of the parameter value to sample
         */
        public <U> GibbsBuilder<V2, S2, T2> addParameterSampler(final V2 parameterName,
                                                                final ParameterSampler<U, V2, S2, T2> parameterSampler,
                                                                final Class<U> parameterValueClass) {
            Utils.nonNull(parameterName);
            Utils.nonNull(parameterSampler);
            Utils.nonNull(parameterValueClass);
            if (samplerMap.containsKey(parameterName)) {
                throw new UnsupportedOperationException("Cannot add more than one sampler per parameter.");
            }
            try {
                state.get(parameterName, parameterValueClass);
            } catch (final IllegalArgumentException e) {
                throw new IllegalArgumentException("Cannot add sampler for parameter that returns type different " +
                        "than that specified for parameter in initial state.");
            }
            samplerMap.put(parameterName, parameterSampler);
            return this;
        }

        /**
         * Builds the {@link ParameterizedModel} as specified via {@link ParameterizedModel.GibbsBuilder}.
         * @return {@link ParameterizedModel} as specified via GibbsBuilder
         * @throws UnsupportedOperationException if there is not a one-to-one mapping between Parameters in the
         *                                       {@link ParameterizedState} and the {@link ParameterSampler}s
         *                                       specified via GibbsBuilder
         */
        public ParameterizedModel<V2, S2, T2> build() {
            if (!samplerMap.keySet().equals(state.keySet())) {
                throw new UnsupportedOperationException("Each parameter must have a corresponding sampler specified.");
            }
            return new ParameterizedModel<>(this);
        }
    }

    /**
     * Constructor using {@link ParameterizedModel.GibbsBuilder}.
     */
    private ParameterizedModel(final GibbsBuilder<V1, S1, T1> builder) {
        state = builder.state;
        dataCollection = builder.dataCollection;
        updateMethod = UpdateMethod.GIBBS;
        updateState = rng -> doGibbsUpdate(builder, rng);
    }

    /**
     * Builder for constructing a {@link ParameterizedModel} to be ensemble sampled using {@link ModelSampler}.
     * Given a list of initial walker positions for all walkers in the ensemble, 
     * an instance "dataset" of a ConcreteDataCollection (which extends {@link DataCollection}), 
     * a function for transforming walker positions to ConcreteParameterizedStates, 
     * and a log-target function to be sampled, this builder holds the necessary state for
     * performing affine-invariant ensemble sampling according to Goodman & Weare 2010.
     * Walkers in the ensemble are iterated over in turn, with the ConcreteParameterizedState held internally 
     * by the {@link ParameterizedModel} corresponding to the currently selected walker; thus, a single iteration
     * over the entire ensemble corresponds to a number of sampled ConcreteParameterizedStates
     * held by the {@link ModelSampler} equal to the number of walkers.  The ConcreteParameterizedState
     * at which the maximum value of the log-target function is observed is also stored.
     * @param <V2>  type of the ParameterEnum
     * @param <S2>  type of the ParameterizedState
     * @param <T2>  type of the DataCollection
     */
    public static final class EnsembleBuilder<V2 extends Enum<V2> & ParameterEnum, S2 extends ParameterizedState<V2>, T2 extends DataCollection> {
        private static final double SCALE_PARAMETER = 2.;

        private final S2 state;
        private final T2 dataCollection;

        private int currentWalkerIndex;
        private int numAccepted;
        private int numSamples;
        private final int numWalkers;
        private final List<WalkerPosition> walkerPositions;
        private final Function<WalkerPosition, S2> transformWalkerPositionToState;
        private final Function<S2, Double> logTarget;
        private final List<S2> currentStates;
        private final List<Double> currentLogTargetValues;
        private double maxLogTarget = Double.NEGATIVE_INFINITY;
        private S2 maxLogTargetState;

        /**
         * Constructor for {@link ParameterizedModel.EnsembleBuilder}.
         * @param initialWalkerPositions            initial positions of the walkers in N-dimensional space
         * @param dataCollection                    DataCollection used by the model
         * @param transformWalkerPositionToState    function that transforms a walker position to a state
         * @param logTarget                         log of the target function to sample
         */
        public EnsembleBuilder(final List<WalkerPosition> initialWalkerPositions,
                               final T2 dataCollection,
                               final Function<WalkerPosition, S2> transformWalkerPositionToState,
                               final Function<S2, Double> logTarget) {
            Utils.nonNull(initialWalkerPositions);
            Utils.nonNull(dataCollection);
            Utils.nonNull(transformWalkerPositionToState);
            Utils.nonNull(logTarget);
            Utils.validateArg(initialWalkerPositions.size() >= 2,
                    "Number of walkers must be greater than or equal to two.");
            Utils.validateArg(initialWalkerPositions.stream()
                    .map(WalkerPosition::numDimensions).allMatch(n -> n.equals(initialWalkerPositions.get(0).numDimensions())),
                    "Dimension of walker space must be identical for all walkers.");

            state = transformWalkerPositionToState.apply(initialWalkerPositions.get(0));
            this.dataCollection = dataCollection;
            currentWalkerIndex = 0;
            numAccepted = 0;
            numSamples = 0;
            numWalkers = initialWalkerPositions.size();
            walkerPositions = new ArrayList<>(initialWalkerPositions);
            this.transformWalkerPositionToState = transformWalkerPositionToState;
            this.logTarget = logTarget;
            //initialize current states and log target values using initial walker positions
            currentStates = walkerPositions.stream().map(transformWalkerPositionToState).collect(Collectors.toList());
            currentLogTargetValues = currentStates.stream().map(logTarget).collect(Collectors.toList());
            maxLogTargetState = currentStates.get(0);
        }

        /**
         * Returns the {@link ParameterizedState} at which the maximum value of the log-target function was observed
         * over the entire sampling run.
         */
        public S2 getMaxLogTargetState() {
            return maxLogTargetState;
        }

        /**
         * Builds the {@link ParameterizedModel} as specified via {@link ParameterizedModel.EnsembleBuilder}.
         * @return {@link ParameterizedModel} as specified via EnsembleBuilder
         */
        public ParameterizedModel<V2, S2, T2> build() {
            return new ParameterizedModel<>(this);
        }
    }

    /**
     * Constructor using {@link ParameterizedModel.EnsembleBuilder}.
     */
    private ParameterizedModel(final EnsembleBuilder<V1, S1, T1> builder) {
        state = builder.state;
        dataCollection = builder.dataCollection;
        updateMethod = UpdateMethod.ENSEMBLE;
        updateState = rng -> doEnsembleUpdate(builder, rng);
    }

    /**
     * Returns a copy of the {@link ParameterizedState} held internally.
     * @return  copy of the {@link ParameterizedState} held internally
     */
    protected S1 state() {
        return state.copy();
    }

    /**
     * Updates the {@link ParameterizedState} held internally using the {@link ParameterSampler}s
     * and update method specified via the Builder pattern.
     * @param rng   {@link RandomGenerator} to pass to {@link ParameterSampler}s to generate samples
     */
    protected void update(final RandomGenerator rng) {
        updateState.accept(rng);
    }

    private void doGibbsUpdate(final GibbsBuilder<V1, S1, T1> builder, final RandomGenerator rng) {
        for (final V1 parameterName : state.keySet()) {
            state.update(parameterName, builder.samplerMap.get(parameterName).sample(rng, state, dataCollection));
        }
    }

    private void doEnsembleUpdate(final EnsembleBuilder<V1, S1, T1> builder, final RandomGenerator rng) {
        //select a walker other than the one currently under consideration for an update
        final int selectedWalkerIndex = IntStream.range(0, builder.numWalkers).boxed()
                .filter(i -> i != builder.currentWalkerIndex)
                .collect(Collectors.toList())
                .get(rng.nextInt(builder.numWalkers - 1));

        //get relevant walker positions
        final WalkerPosition currentWalkerPosition = builder.walkerPositions.get(builder.currentWalkerIndex);
        final WalkerPosition selectedWalkerPosition = builder.walkerPositions.get(selectedWalkerIndex);

        //propose a stretch move
        final int numDimensions = currentWalkerPosition.numDimensions();
        final double z = Math.pow((EnsembleBuilder.SCALE_PARAMETER - 1.) * rng.nextDouble() + 1, 2.) / EnsembleBuilder.SCALE_PARAMETER; //see Goodman & Weare 2010
        final WalkerPosition proposedWalkerPosition = new WalkerPosition(IntStream.range(0, numDimensions).boxed()
                .map(i -> selectedWalkerPosition.get(i) + z * (currentWalkerPosition.get(i) - selectedWalkerPosition.get(i)))
                .collect(Collectors.toList()));

        //transform to states
        final S1 currentState = builder.currentStates.get(builder.currentWalkerIndex);
        currentState.values().forEach(p -> logger.debug("Current " + p.getName().name() + ": " + p.getValue()));
        final S1 proposedState = builder.transformWalkerPositionToState.apply(proposedWalkerPosition);
        proposedState.values().forEach(p -> logger.debug("Proposed " + p.getName().name() + ": " + p.getValue()));

        //calculate log targets
        final double currentLogTarget = builder.currentLogTargetValues.get(builder.currentWalkerIndex);
        final double proposedLogTarget = builder.logTarget.apply(proposedState);

        //accept or reject
        final double acceptanceProbability = Math.min(1., Math.exp((numDimensions - 1.) * Math.log(z) + proposedLogTarget - currentLogTarget)); //see Goodman & Weare 2010
        logger.debug("Log target of current state: " + currentLogTarget);
        logger.debug("Log target of proposed state: " + proposedLogTarget);
        builder.numSamples++;
        if (rng.nextDouble() < acceptanceProbability) {
            builder.numAccepted++;
            logger.debug("Proposed state accepted.");
            //update the walker position for the current walker
            builder.walkerPositions.set(builder.currentWalkerIndex, proposedWalkerPosition);
            //update the current log target value for the current walker
            builder.currentLogTargetValues.set(builder.currentWalkerIndex, proposedLogTarget);
            //update the state held by the model using the accepted proposed state for the current walker
            proposedState.values().forEach(p -> state.update(p.getName(), p.getValue()));
            //update the maximum target state if appropriate
            if (proposedLogTarget > builder.maxLogTarget) {
                logger.debug("New maximum found.");
                builder.maxLogTarget = proposedLogTarget;
                builder.maxLogTargetState = proposedState;
            }
        } else {
            //update the state held by the model using the previous state for the current walker
            currentState.values().forEach(p -> state.update(p.getName(), p.getValue()));
        }
        logger.debug("Acceptance rate: " + (double) builder.numAccepted / builder.numSamples);
        state.values().forEach(p -> logger.debug("Sampled " + p.getName().name() + ": " + p.getValue()));

        //move to the next walker
        builder.currentWalkerIndex = (builder.currentWalkerIndex + 1) % builder.numWalkers;
    }
}
