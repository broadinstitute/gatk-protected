package org.broadinstitute.hellbender.utils.mcmc;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;

public final class GibbsSamplerSimpleUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/utils/mcmc/";

    private static final File COVERAGES_FILE = new File(TEST_SUB_DIR
            + "coverages-for-gibbs-sampler-one-segment-copy-ratio-test.txt");

    private static final double VARIANCE_MIN = 0.;
    private static final double VARIANCE_MAX = 10.;
    private static final double VARIANCE_WIDTH = 0.1;
    private static final double VARIANCE_INITIAL = 5.;
    private static final double VARIANCE_TRUTH = 1.;
    private static final double VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH = 0.014;

    private static final double MEAN_MIN = 0.;
    private static final double MEAN_MAX = 10.;
    private static final double MEAN_WIDTH = 0.1;
    private static final double MEAN_INITIAL = 5.;
    private static final double MEAN_TRUTH = 1.;
    private static final double MEAN_POSTERIOR_STANDARD_DEVIATION_TRUTH = 0.01;

    private static final List<Parameter> INITIAL_PARAMETER_VALUES =
            Arrays.asList(new Parameter<>("variance", VARIANCE_INITIAL), new Parameter<>("mean", MEAN_INITIAL));
    
    private static final int NUM_SAMPLES = 500;
    private static final int NUM_BURN_IN = 250;

    private static double normalTerm(final double quantity, final double mean, final double variance) {
        return (quantity - mean) * (quantity - mean) / (2 * variance);
    }

    private static double relativeError(final double x, final double y) {
        return Math.abs(x - y) / (2 * (x + y));
    }

    private final class CopyRatioModel {
        private final ParameterizedModel<CopyRatioState, CopyRatioData> model;

        private CopyRatioModel(final CopyRatioState state, final CopyRatioData data) {
            model = new ParameterizedModel.GibbsBuilder<>(state, data)
                    .addParameterSampler("variance", varianceSampler)
                    .addParameterSampler("mean", meanSampler)
                    .build();
        }

        private final Sampler<Parameter, CopyRatioState, CopyRatioData> varianceSampler = (rng, state, data) -> {
            final Function<Double, Double> logConditionalPDF =
                    newVariance -> -0.5 * Math.log(newVariance) * data.coverages.size() +
                            data.coverages.values().stream()
                                    .mapToDouble(c -> -normalTerm(c, state.get("mean"), newVariance))
                                    .sum();

            final SliceSampler sampler =
                    new SliceSampler(rng, logConditionalPDF, state.get("variance"), VARIANCE_MIN, VARIANCE_MAX, VARIANCE_WIDTH);
            return new Parameter<>("variance", sampler.sample());
        };

        private final Sampler<Parameter, CopyRatioState, CopyRatioData> meanSampler = (rng, state, data) -> {
            final Function<Double, Double> logConditionalPDF =
                    newMean -> data.coverages.values().stream()
                            .mapToDouble(c -> -normalTerm(c, newMean, state.get("variance")))
                            .sum();

            final SliceSampler sampler =
                    new SliceSampler(rng, logConditionalPDF, state.get("mean"), MEAN_MIN, MEAN_MAX, MEAN_WIDTH);
            return new Parameter<>("mean", sampler.sample());
        };
    }

    private final class CopyRatioState extends ParameterizedState {
        public CopyRatioState(final List<Parameter> parameters) {
            super(parameters);
        }

        public <T> T get(final String parameterName) {
            return (T) getParameter(parameterName).value();
        }
    }

    private final class CopyRatioData extends DataCollection {
        private final Data<Double> coverages;

        private CopyRatioData(final Data<Double> coverages) {
            super(Arrays.asList(coverages));
            this.coverages = coverages;
        }

        private CopyRatioData(final File coveragesFile) {
            this(Data.loadData(coveragesFile, Double::parseDouble));
        }
    }

    @Test
    public void testRunMCMCOnCopyRatioModel() {
        final CopyRatioState state = new CopyRatioState(INITIAL_PARAMETER_VALUES);
        final CopyRatioData data = new CopyRatioData(COVERAGES_FILE);
        final CopyRatioModel model = new CopyRatioModel(state, data);
        final GibbsSampler<CopyRatioState, CopyRatioData> gibbsSampler = new GibbsSampler<>(NUM_SAMPLES, model.model);
        gibbsSampler.runMCMC();

        final double[] varianceSamples = Doubles.toArray(gibbsSampler.getSamples("variance", NUM_BURN_IN));
        final double variancePosteriorMean = new Mean().evaluate(varianceSamples);
        final double variancePosteriorStandardDeviation = new StandardDeviation().evaluate(varianceSamples);

        Assert.assertEquals(relativeError(variancePosteriorMean, VARIANCE_TRUTH), 0., 0.01);
        Assert.assertEquals(
                relativeError(variancePosteriorStandardDeviation, VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH),
                0., 0.05);

        final double[] meanSamples = Doubles.toArray(gibbsSampler.getSamples("mean", NUM_BURN_IN));
        final double meanPosteriorMean = new Mean().evaluate(meanSamples);
        final double meanPosteriorStandardDeviation = new StandardDeviation().evaluate(meanSamples);
        Assert.assertEquals(relativeError(meanPosteriorMean, MEAN_TRUTH), 0., 0.01);
        Assert.assertEquals(
                relativeError(meanPosteriorStandardDeviation, MEAN_POSTERIOR_STANDARD_DEVIATION_TRUTH),
                0., 0.05);
    }
}