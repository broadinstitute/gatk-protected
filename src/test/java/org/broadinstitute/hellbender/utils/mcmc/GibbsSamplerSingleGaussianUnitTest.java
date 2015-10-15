package org.broadinstitute.hellbender.utils.mcmc;

import com.google.common.primitives.Doubles;
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

public final class GibbsSamplerSingleGaussianUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/utils/mcmc/";

    private static final File COVERAGES_FILE = new File(TEST_SUB_DIR
            + "datapoints-for-gibbs-sampler-single-gaussian-test.txt");

    private static final double VARIANCE_MIN = 0.;
    private static final double VARIANCE_MAX = Double.POSITIVE_INFINITY;
    private static final double VARIANCE_WIDTH = 0.1;
    private static final double VARIANCE_INITIAL = 5.;
    private static final double VARIANCE_TRUTH = 1.;
    private static final double VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH = 0.014;

    private static final double MEAN_WIDTH = 0.1;
    private static final double MEAN_INITIAL = 5.;
    private static final double MEAN_TRUTH = 1.;
    private static final double MEAN_POSTERIOR_STANDARD_DEVIATION_TRUTH = 0.01;

    private static final int NUM_SAMPLES = 500;
    private static final int NUM_BURN_IN = 250;

    private static double normalTerm(final double quantity, final double mean, final double variance) {
        return (quantity - mean) * (quantity - mean) / (2 * variance);
    }

    private static double relativeError(final double x, final double xTrue) {
        return Math.abs((x - xTrue) / xTrue);
    }

    private final class GaussianModeller {
        private final ParameterizedModel<ParameterizedState, DataCollection> model;
        private final Sampler<Double, ParameterizedState, DataCollection> varianceSampler;
        private final Sampler<Double, ParameterizedState, DataCollection> meanSampler;

        private GaussianModeller(final double varianceInitial, final double meanInitial, final File datapointsFile) {
            final List<Parameter<?>> initialParameters =
                    Arrays.asList(new Parameter<>("variance", varianceInitial), new Parameter<>("mean", meanInitial));
            final ParameterizedState initialState = new ParameterizedState(initialParameters);

            final Data<Double> datapoints = new Data<>("datapoints", datapointsFile, Double::parseDouble);
            final DataCollection dataset = new DataCollection(Collections.singletonList(datapoints));

            varianceSampler = (rng, state, data) -> {
                final Function<Double, Double> logConditionalPDF =
                        newVariance -> -0.5 * Math.log(newVariance) * data.<Double>get("datapoints").size() +
                                data.<Double>get("datapoints").stream()
                                        .mapToDouble(c -> -normalTerm(c, state.get("mean", Double.class), newVariance))
                                        .sum();

                final SliceSampler sampler =
                        new SliceSampler(rng, logConditionalPDF, state.get("variance", Double.class),
                                VARIANCE_MIN, VARIANCE_MAX, VARIANCE_WIDTH);
                return sampler.sample();
            };

            meanSampler = (rng, state, data) -> {
                final Function<Double, Double> logConditionalPDF =
                        newMean -> data.<Double>get("datapoints").stream()
                                .mapToDouble(c -> -normalTerm(c, newMean, state.get("variance", Double.class)))
                                .sum();

                final SliceSampler sampler =
                        new SliceSampler(rng, logConditionalPDF, state.get("mean", Double.class), MEAN_WIDTH);
                return sampler.sample();
            };

            model = new ParameterizedModel.GibbsBuilder<>(initialState, dataset)
                    .addParameterSampler("variance", varianceSampler)
                    .addParameterSampler("mean", meanSampler)
                    .build();
        }
    }

    @Test
    public void testRunMCMCOnGaussianModel() {
        final GaussianModeller modeller = new GaussianModeller(VARIANCE_INITIAL, MEAN_INITIAL, COVERAGES_FILE);
        final GibbsSampler<ParameterizedState, DataCollection> gibbsSampler =
                new GibbsSampler<>(NUM_SAMPLES, modeller.model, ParameterizedState.class);
        gibbsSampler.runMCMC();

        final double[] varianceSamples = Doubles.toArray(gibbsSampler.getSamples("variance", Double.class, NUM_BURN_IN));
        final double variancePosteriorMean = new Mean().evaluate(varianceSamples);
        final double variancePosteriorStandardDeviation = new StandardDeviation().evaluate(varianceSamples);
        Assert.assertEquals(relativeError(variancePosteriorMean, VARIANCE_TRUTH), 0., 0.01);
        Assert.assertEquals(
                relativeError(variancePosteriorStandardDeviation, VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH),
                0., 0.1);

        final double[] meanSamples = Doubles.toArray(gibbsSampler.getSamples("mean", Double.class, NUM_BURN_IN));
        final double meanPosteriorMean = new Mean().evaluate(meanSamples);
        final double meanPosteriorStandardDeviation = new StandardDeviation().evaluate(meanSamples);
        Assert.assertEquals(relativeError(meanPosteriorMean, MEAN_TRUTH), 0., 0.01);
        Assert.assertEquals(
                relativeError(meanPosteriorStandardDeviation, MEAN_POSTERIOR_STANDARD_DEVIATION_TRUTH),
                0., 0.1);
    }
}