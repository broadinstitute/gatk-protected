package org.broadinstitute.hellbender.utils.mcmc;


import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Implements Gibbs sampling of a multivariate probability density function.  Variables are separated into global
 * and local (i.e., segment-level or site-level) parameters.  The probability density function is assumed to be
 * independent across segments/sites for each local parameter.  Univariate conditional probability density functions are
 * sampled using {@link SliceSampler} and are assumed to be unimodal.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class GibbsSampler {
    private static final int RANDOM_SEED = 42;
    private static final RandomGenerator rng =
            RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

    private final Logger logger = LogManager.getLogger(GibbsSampler.class);
    private static final int NUMBER_OF_SAMPLES_PER_LOG_ENTRY = 100;

    private final int numSamples;

    private final List<List<Double>> data;                  //indexed by (dataset, datapoint)

    private final List<List<Double>> globalSamples;         //indexed by (sample number, parameter)
    private final List<List<List<Double>>> localSamples;    //indexed by (sample number, parameter, segment/site)

    private final List<GlobalSampler> globalSamplers;
    private final List<LocalSampler> localSamplers;

    private final int numGlobalParameters;
    private final int numLocalParameters;

    private boolean isMCMCRunComplete = false;

    /**
     * Class for constructing a GibbsSampler using the Builder pattern.
     */
    public static final class Builder {
        private final int numSamples;
        private final List<List<Double>> data;

        private final List<Double> initialGlobalParameters = new ArrayList<>();
        private final List<List<Double>> initialLocalParameters = new ArrayList<>();
        private final List<GlobalSampler> globalSamplers = new ArrayList<>();
        private final List<LocalSampler> localSamplers = new ArrayList<>();

        /**
         * Builder method for constructing a {@link GibbsSampler}. Example use:
         * <p>
         *     GibbsSampler gibbSampler = new GibbsSampler.Builder(numberOfSamplers, data)
         *                                                .addGlobalSampler(globalSampler1)
         *                                                ...
         *                                                .addGlobalSampler(globalSamplerN)
         *                                                .addLocalSampler(localSampler1)
         *                                                ...
         *                                                .addLocalSampler(localSamplerN)
         *                                                .initializeGlobalParameters(initialGlobalParameters)
         *                                                .initializeLocalParameters(initialLocalParameters)
         *                                                .build()
         * </p>
         * @param numSamples    total number of samples to generate (including burn-in)
         * @param data          list of datasets, which are themselves lists of double-type datapoints,
         *                      containing all information required to calculate the probability density function
         *                      to be sampled.
         *                      e.g., a list containing a list of coverages and a list of number of targets within
         *                      each segment suffices to calculate the likelihood in a simple copy-ratio model
         */
        public Builder(final int numSamples, List<List<Double>> data) {
            if (numSamples <= 0) {
                throw new IllegalArgumentException("Number of samples must be positive.");
            }
            Utils.nonNull("List of all data cannot be null.");
            for (final List<Double> datum : data) {
                Utils.nonNull(datum, "Data cannot be null.");
            }
            this.numSamples = numSamples;
            this.data = data;
        }

        /**
         * Initialize values for global parameters.
         * @param globalParameters  list of initial values for all global parameters
         * @return                  GibbsSampler Builder
         */
        public Builder initializeGlobalParameters(final List<Double> globalParameters) {
            Utils.nonNull(globalParameters, "List of initial global-parameter values cannot be null.");
            initialGlobalParameters.addAll(globalParameters);
            return this;
        }

        /**
         * Initialize values for local parameters.
         * @param localParameters   list of lists of initial values for all local parameters at all segments/sites
         * @return                  GibbsSampler Builder
         */
        public Builder initializeLocalParameters(final List<List<Double>> localParameters) {
            Utils.nonNull(localParameters, "List of lists of initial local-parameter values cannot be null.");
            for (final List<Double> localParameter : localParameters) {
                Utils.nonNull(localParameter, "List of initial local-parameter values cannot be null.");
            }
            initialLocalParameters.addAll(localParameters);
            return this;
        }

        /**
         * Add a GlobalSampler that samples the conditional probability density function for a global parameter.
         * @param globalSampler     GlobalSampler for conditional probability density function for global parameter
         * @return                  GibbsSampler Builder
         */
        public Builder addGlobalSampler(final GlobalSampler globalSampler) {
            Utils.nonNull(globalSampler, "Global sampler cannot be null.");
            globalSamplers.add(globalSampler);
            return this;
        }

        /**
         * Add a LocalSampler that samples the conditional probability density function for a local parameter at
         * all segments/sites.
         * @param localSampler      LocalSampler for conditional probability density function for local parameter
         * @return                  GibbsSampler Builder
         */
        public Builder addLocalSampler(final LocalSampler localSampler) {
            Utils.nonNull(localSampler, "Local sampler cannot be null.");
            localSamplers.add(localSampler);
            return this;
        }

        /**
         * Return the GibbsSampler constructed using the Builder pattern.
         * @return                  GibbsSampler
         * @see                     GibbsSampler.Builder#Builder(int, List)
         */
        public GibbsSampler build() {
            if (initialGlobalParameters.size() != globalSamplers.size()) {
                throw new IllegalArgumentException("Number of initial global-parameter values does not match " +
                        "number of global samplers.");
            }
            if (initialLocalParameters.size() != localSamplers.size()) {
                throw new IllegalArgumentException("Number of initial local-parameter values does not match " +
                        "number of local samplers.");
            }
            if (initialGlobalParameters.size() == 0 && initialLocalParameters.size() == 0) {
                throw new IllegalArgumentException("Total number of global and local parameters cannot be zero.");
            }
            return new GibbsSampler(this);
        }
    }

    private GibbsSampler(final Builder builder) {
        numSamples = builder.numSamples;

        data = builder.data;

        globalSamples = new ArrayList<>(numSamples);
        localSamples = new ArrayList<>(numSamples);
        globalSamples.add(builder.initialGlobalParameters);
        localSamples.add(builder.initialLocalParameters);

        globalSamplers = builder.globalSamplers;
        localSamplers = builder.localSamplers;

        numGlobalParameters = globalSamplers.size();
        numLocalParameters = localSamplers.size();
    }

    /**
     * Run Markov-Chain Monte-Carlo sampling of global and local parameters using the total number of samples,
     * initial parameter values, and samplers specified by {@link GibbsSampler.Builder}.
     */
    public void runMCMC() {
        logger.info("Initializing global and local parameters.");
        final List<Double> currentGlobalParameters = new ArrayList<>(globalSamples.get(0));
        final List<List<Double>> currentLocalParameters = new ArrayList<>(localSamples.get(0));

        //do one iteration to check that all LocalSamplers return lists of correct length
        //(i.e., consistent with initial values)
        logger.info("Checking multiplicity of local parameters.");
        for (int localParameter = 0; localParameter < numLocalParameters; localParameter++) {
            final List<Double> newLocalParameterSample =
                    localSamplers.get(localParameter)
                            .sample(rng, currentGlobalParameters, currentLocalParameters, data);
            if (newLocalParameterSample.size() != currentLocalParameters.get(localParameter).size()) {
                throw new IllegalStateException("Number of initial values for local parameter " + localParameter +
                        " does not match number of values returned by LocalSampler. " +
                        "Check that the order of LocalSamplers in GibbsSampler builder is the same as the " +
                        "order of parameters in initial values.");
            }
        }

        logger.info("Starting MCMC sampling.");
        for (int sample = 1; sample < numSamples; sample++) {
            if (sample % NUMBER_OF_SAMPLES_PER_LOG_ENTRY == 0) {
                logger.info(sample + " of " + numSamples + " samples generated.");
            }
            globalSamples.add(new ArrayList<>(numGlobalParameters));
            localSamples.add(new ArrayList<>(numLocalParameters));
            for (int globalParameter = 0; globalParameter < numGlobalParameters; globalParameter++) {
                final double newGlobalParameterSample =
                        globalSamplers.get(globalParameter)
                                .sample(rng, currentGlobalParameters, currentLocalParameters, data);
                globalSamples.get(sample).add(newGlobalParameterSample);
                currentGlobalParameters.set(globalParameter, newGlobalParameterSample);
            }
            for (int localParameter = 0; localParameter < numLocalParameters; localParameter++) {
                final List<Double> newLocalParameterSample =
                        localSamplers.get(localParameter)
                                .sample(rng, currentGlobalParameters, currentLocalParameters, data);
                localSamples.get(sample).add(newLocalParameterSample);
                currentLocalParameters.set(localParameter, newLocalParameterSample);
            }
        }
        logger.info(numSamples + " of " + numSamples + " samples generated.");
        logger.info("MCMC sampling complete.");
        isMCMCRunComplete = true;
    }

    /**
     * Return a deep copy of the global-parameter samples held internally, discarding a specified number of
     * burn-in samples.  If MCMC sampling was not run previously, {@link GibbsSampler#runMCMC} will be executed first.
     * @param numBurnIn     number of burn-in samples
     * @return              list of lists of global-parameter samples, indexed by (sample number, parameter)
     */
    public List<List<Double>> getGlobalSamples(final int numBurnIn) {
        if (numBurnIn < 0 || numBurnIn >= numSamples) {
            throw new IllegalArgumentException("Number of burn-in samples must be non-negative and less than " +
                    "total number of samples");
        }
        if (!isMCMCRunComplete) {
            runMCMC();
        }
        final List<List<Double>> globalSamplesCopy = new ArrayList<>();
        for (int sample = numBurnIn; sample < numSamples; sample++) {
            globalSamplesCopy.add(new ArrayList<>(globalSamples.get(sample)));
        }
        return globalSamplesCopy;
    }

    /**
     * Return a deep copy of the local-parameter samples held internally, discarding a specified number of
     * burn-in samples.  If MCMC sampling was not run previously, {@link GibbsSampler#runMCMC} will be executed first.
     * @param numBurnIn     number of burn-in samples
     * @return              list of list of lists of local-parameter samples at each segment/site,
     *                      indexed by (sample number, parameter, segment/site)
     */
    public List<List<List<Double>>> getLocalSamples(final int numBurnIn) {
        if (numBurnIn < 0 || numBurnIn >= numSamples) {
            throw new IllegalArgumentException("Number of burn-in samples must be non-negative and less than " +
                    "total number of samples");
        }
        if (!isMCMCRunComplete) {
            runMCMC();
        }
        final List<List<List<Double>>> localSamplesCopy = new ArrayList<>();
        for (int sample = numBurnIn; sample < numSamples; sample++) {
            final List<List<Double>> localSampleCopy = new ArrayList<>();
            for (int localParameter = 0; localParameter < numLocalParameters; localParameter++) {
                localSampleCopy.add(new ArrayList<>(localSamples.get(sample).get(localParameter)));
            }
            localSamplesCopy.add(localSampleCopy);
        }
        return localSamplesCopy;
    }
}