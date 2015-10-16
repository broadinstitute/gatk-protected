package org.broadinstitute.hellbender.utils.mcmc;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

/**
 * Unit test for {@link GibbsSampler}.
 * <p>
 *     Test performs Bayesian inference of a Gaussian model with 2 global parameters specifying the variance and the mean.
 * </p>
 * <p>
 *     Data consists of a list of 10,000 datapoints drawn from a normal distribution with unity variance and mean.
 * </p>
 * <p>
 *     Success of the test is determined by recovery of the input variance and mean,
 *     as well as agreement of the standard deviations of the parameter posteriors with those given by both the
 *     python package emcee (see http://dan.iel.fm/emcee for details) and numerical evaluation in Mathematica
 *     of the analytic forms of the posteriors.
 * </p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class GibbsSamplerSingleGaussianUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/utils/mcmc/";

    private static final File COVERAGES_FILE = new File(TEST_SUB_DIR
            + "datapoints-for-gibbs-sampler-single-gaussian-test.txt");
    private static final String DATAPOINTS_NAME = "datapoints";

    private static final String VARIANCE_NAME = "variance";
    private static final double VARIANCE_MIN = 0.;
    private static final double VARIANCE_MAX = Double.POSITIVE_INFINITY;
    private static final double VARIANCE_WIDTH = 0.1;
    private static final double VARIANCE_INITIAL = 5.;
    private static final double VARIANCE_TRUTH = 1.;
    private static final double VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH = 0.014;

    private static final String MEAN_NAME = "mean";
    private static final double MEAN_WIDTH = 0.1;
    private static final double MEAN_INITIAL = 5.;
    private static final double MEAN_TRUTH = 1.;
    private static final double MEAN_POSTERIOR_STANDARD_DEVIATION_TRUTH = 0.01;

    private static final int NUM_SAMPLES = 500;
    private static final int NUM_BURN_IN = 250;

    /**
     * Calculates the exponent for a normal distribution; used in log-likelihood calculation below.
     */
    private static double normalTerm(final double quantity, final double mean, final double variance) {
        return (quantity - mean) * (quantity - mean) / (2 * variance);
    }

    private static double relativeError(final double x, final double xTrue) {
        return Math.abs((x - xTrue) / xTrue);
    }

    /**
     * Helper class to organize the model and its parameter samplers.
     */
    private final class GaussianModeller {
        //create fields in the Modeller for the model and samplers
        private final ParameterizedModel<ParameterizedState, DataCollection> model;
        private final Sampler<Double, ParameterizedState, DataCollection> varianceSampler;
        private final Sampler<Double, ParameterizedState, DataCollection> meanSampler;

        /**
         * Constructor for the Modeller takes as parameters all quantities needed to construct the ParameterizedState
         * (here, the initial variance and the initial mean) and the DataCollection (here, the datapoints file).
         */
        private GaussianModeller(final double varianceInitial, final double meanInitial, final File datapointsFile) {
            //construct the initial ParameterizedState by passing a list of Parameters of mixed type to the constructor
            //the names, initial values, and types of each of the parameters are set here
            final List<Parameter<?>> initialParameters =
                    Arrays.asList(new Parameter<>(VARIANCE_NAME, varianceInitial), new Parameter<>(MEAN_NAME, meanInitial));
            final ParameterizedState initialState = new ParameterizedState(initialParameters);

            //construct the DataCollection by passing a list of Data objects
            //the names, values, and types of each of the Data objects are set here
            final Data<Double> datapoints = new Data<>(DATAPOINTS_NAME, datapointsFile, Double::parseDouble);
            final DataCollection dataset = new DataCollection(Collections.singletonList(datapoints));

            //below, we override sample() for each of the Samplers; this can be done via a lambda that takes
            //(rng, state, data) and returns a new sample with type identical to that specified during initialization

            //Sampler for the variance global parameter.  The relevant logConditionalPDF is given by
            //the log of the product of Gaussian likelihoods for each datapoint c_t:
            //      log[product_t variance^(-1/2) * exp(-(c_t - mean)^2 / (2 * variance))] + constant
            //which reduces to the form in code below.  Slice sampling is used here to generate a new sample,
            //but, in general, any method can be used; e.g., if the conditional is from the exponential family,
            //one could simply use the corresponding Distribution from Apache Commons.
            varianceSampler = (rng, state, data) -> {
                final Function<Double, Double> logConditionalPDF =
                        newVariance -> -0.5 * Math.log(newVariance) * data.<Double>get(DATAPOINTS_NAME).size() +
                                data.<Double>get(DATAPOINTS_NAME).stream()
                                        .mapToDouble(c -> -normalTerm(c, state.get(MEAN_NAME, Double.class), newVariance))
                                        .sum();

                final SliceSampler sampler =
                        new SliceSampler(rng, logConditionalPDF, state.get(VARIANCE_NAME, Double.class),
                                VARIANCE_MIN, VARIANCE_MAX, VARIANCE_WIDTH);
                return sampler.sample();
            };

            //Sampler for the mean global parameter.  The relevant logConditionalPDF is given by
            //the log of the product of Gaussian likelihoods for each datapoint c_t:
            //     log[product_t exp(-(c_t - mean)^2 / (2 * variance))] + constant
            //which reduces to the form in code below.
            meanSampler = (rng, state, data) -> {
                final Function<Double, Double> logConditionalPDF =
                        newMean -> data.<Double>get(DATAPOINTS_NAME).stream()
                                .mapToDouble(c -> -normalTerm(c, newMean, state.get(VARIANCE_NAME, Double.class)))
                                .sum();

                final SliceSampler sampler =
                        new SliceSampler(rng, logConditionalPDF, state.get(MEAN_NAME, Double.class), MEAN_WIDTH);
                return sampler.sample();
            };

            //build the ParameterizedModel using the GibbsBuilder pattern
            //pass in the initial ParameterizedState, DataCollection, and specify the class of the ParameterizedState
            //then add samplers for each of the parameters, with names matching those used in initialization
            model = new ParameterizedModel.GibbsBuilder<>(initialState, dataset, ParameterizedState.class)
                    .addParameterSampler(VARIANCE_NAME, varianceSampler)
                    .addParameterSampler(MEAN_NAME, meanSampler)
                    .build();
        }
    }

    @Test
    public void testRunMCMCOnGaussianModel() {
        //create new instance of the Modeller helper class, passing all quantities needed to initialize state and data
        final GaussianModeller modeller = new GaussianModeller(VARIANCE_INITIAL, MEAN_INITIAL, COVERAGES_FILE);
        //create a GibbsSampler, passing the total number of samples (including burn-in samples) and the model held by the Modeller
        final GibbsSampler<ParameterizedState, DataCollection> gibbsSampler =
                new GibbsSampler<>(NUM_SAMPLES, modeller.model);
        //run the MCMC
        gibbsSampler.runMCMC();

        //get the samples of each of the parameter posteriors (discarding burn-in samples) by passing the
        //parameter name, type, and burn-in number to the getSamples method
        final double[] varianceSamples = Doubles.toArray(gibbsSampler.getSamples(VARIANCE_NAME, Double.class, NUM_BURN_IN));
        final double[] meanSamples = Doubles.toArray(gibbsSampler.getSamples(MEAN_NAME, Double.class, NUM_BURN_IN));

        //check that the statistics of the posteriors agree with those found by emcee/analytically
        final double variancePosteriorMean = new Mean().evaluate(varianceSamples);
        final double variancePosteriorStandardDeviation = new StandardDeviation().evaluate(varianceSamples);
        Assert.assertEquals(relativeError(variancePosteriorMean, VARIANCE_TRUTH), 0., 0.01);
        Assert.assertEquals(
                relativeError(variancePosteriorStandardDeviation, VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH),
                0., 0.1);
        final double meanPosteriorMean = new Mean().evaluate(meanSamples);
        final double meanPosteriorStandardDeviation = new StandardDeviation().evaluate(meanSamples);
        Assert.assertEquals(relativeError(meanPosteriorMean, MEAN_TRUTH), 0., 0.01);
        Assert.assertEquals(
                relativeError(meanPosteriorStandardDeviation, MEAN_POSTERIOR_STANDARD_DEVIATION_TRUTH),
                0., 0.1);
    }
}