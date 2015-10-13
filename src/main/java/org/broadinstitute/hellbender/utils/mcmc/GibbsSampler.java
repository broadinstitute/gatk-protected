package org.broadinstitute.hellbender.utils.mcmc;


import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

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

    private static final Logger logger = LogManager.getLogger(GibbsSampler.class);
    private static final int NUMBER_OF_SAMPLES_PER_LOG_ENTRY = 100;

    private final int numSamples;

    private final ParameterizedModel model;

    private final List<ParameterizedState> samples;

    private boolean isMCMCRunComplete = false;

    public GibbsSampler(final int numSamples, final ParameterizedModel model) {
        if (numSamples <= 0) {
            throw new IllegalArgumentException("Number of samples must be positive.");
        }
        this.numSamples = numSamples;
        this.model = model;
        samples = new ArrayList<>(numSamples);
        samples.add(model.state());
    }

    public void runMCMC() {
        logger.info("Starting MCMC sampling.");
        for (int sample = 1; sample < numSamples; sample++) {
            if (sample % NUMBER_OF_SAMPLES_PER_LOG_ENTRY == 0) {
                logger.info(sample + " of " + numSamples + " samples generated.");
            }
            model.update(rng);
            samples.add(model.state());
        }
        logger.info(numSamples + " of " + numSamples + " samples generated.");
        logger.info("MCMC sampling complete.");
        isMCMCRunComplete = true;
    }

    public <T> List<T> getSamples(final String parameterName, final Class<T> type, final int numBurnIn) {
        if (numBurnIn < 0 || numBurnIn >= numSamples) {
            throw new IllegalArgumentException("Invalid number of burn-in samples.");
        }
        if (!isMCMCRunComplete) {
            runMCMC();
        }
        return samples.stream().map(s -> s.get(parameterName, type)).collect(Collectors.toList()).subList(numBurnIn, numSamples);
    }
}