package org.broadinstitute.hellbender.utils.mcmc;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Unit test for {@link GibbsSampler}.  Demonstrates application of {@link GibbsSampler} to a {@ParameterizedModel}
 * that is specified using helper classes extended from {@ParameterizedState} and {@DataCollection}.
 * <p>
 *     Test performs Bayesian inference of a simple copy-ratio model with 1 global and 100 segment-level parameters.
 *     We consider an exome with 100 segments, each of which has a number of targets that is drawn from a
 *     Poisson distribution with mean 100 in the generated data set.
 * </p>
 * <p>
 *     Data consists of a list of coverages at each target and a list of the number of targets in each segment.
 *     Coverages are assumed to be drawn from a normal distribution in each segment.
 * </p>
 * <p>
 *     The mean of each normal distribution is given by a segment-level parameter;
 *     the distribution of means across segments is taken to be uniform in (0, 10) in the generated data.
 * </p>
 * <p>
 *     The variance of each distribution is set by the global parameter, which is taken to be 1 in the generated data.
 * </p>
 * <p>
 *     Success of the test is determined by recovery of the input segment means and the global variance,
 *     as well as agreement of the standard deviations of the parameter posteriors with those given by both the
 *     python package emcee (see http://dan.iel.fm/emcee for details) and numerical evaluation in Mathematica
 *     of the analytic forms of the posteriors.
 * </p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class GibbsSamplerCopyRatioUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/utils/mcmc/";

    private static final File COVERAGES_FILE = new File(TEST_SUB_DIR
            + "coverages-for-gibbs-sampler-copy-ratio-test.txt");
    private static final File NUM_TARGETS_PER_SEGMENT_FILE =
            new File(TEST_SUB_DIR + "number-of-targets-per-segment-for-gibbs-sampler-copy-ratio-test.txt");
    private static final File MEANS_TRUTH_FILE = new File(TEST_SUB_DIR
            + "means-truth-for-gibbs-sampler-copy-ratio-test.txt");

    private static final double VARIANCE_MIN = 0.;
    private static final double VARIANCE_MAX = Double.POSITIVE_INFINITY;
    private static final double VARIANCE_WIDTH = 0.1;
    private static final double VARIANCE_INITIAL = 5.;
    private static final double VARIANCE_TRUTH = 1.;
    private static final double VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH = 0.014;    //checked with emcee & analytic result

    private static final double MEAN_WIDTH = 0.1;
    private static final double MEAN_INITIAL = 5.;
    private static final double MEAN_POSTERIOR_STANDARD_DEVIATION_MEAN_TRUTH = 0.1;     //checked with emcee & analytic result

    private static final int NUM_SAMPLES = 500;
    private static final int NUM_BURN_IN = 250;

    /**
     * Calculates the exponent for a normal distribution; used in log-likelihood calculation below.
     */
    private static double normalTerm(final double quantity, final double mean, final double variance) {
        return (quantity - mean) * (quantity - mean) / (2 * variance);
    }

    /**
     * Calculates relative error between x and xTrue, with respect to xTrue; used for checking statistics of
     * posterior samples below.
     */
    private static double relativeError(final double x, final double xTrue) {
        return Math.abs((x - xTrue) / xTrue);
    }

    //In contrast to

    private final class SegmentMeans extends AbstractParameterizedState {
        private static final String MEAN_IN_SEGMENT_PREFIX = "meanInSegment";

        public SegmentMeans(final List<Double> means) {
            super(MEAN_IN_SEGMENT_PREFIX, means);
        }

        public SegmentMeans(final SegmentMeans segmentMeans) {
            super(segmentMeans);
        }

        @Override
        protected <S extends AbstractParameterizedState> S copy(final Class<S> stateClass) {
            return stateClass.cast(new SegmentMeans(this));
        }

        public double getMeanInSegment(final int segment) {
            return get(MEAN_IN_SEGMENT_PREFIX + segment, Double.class);
        }
    }

    private final class CopyRatioState extends AbstractParameterizedState {
        private static final String VARIANCE_NAME = "variance";
        private static final String MEANS_NAME = "means";

        public CopyRatioState(final double variance, final double mean, final int numSegments) {
            super(Arrays.asList(
                    new Parameter<>(VARIANCE_NAME, variance),
                    new Parameter<>(MEANS_NAME, new SegmentMeans(Collections.nCopies(numSegments, mean)))));
        }

        public CopyRatioState(final CopyRatioState state) {
            super(state);
        }

        @Override
        protected <S extends AbstractParameterizedState> S copy(final Class<S> stateClass) {
            return stateClass.cast(new CopyRatioState(this));
        }

        public double variance() {
            return get(VARIANCE_NAME, Double.class);
        }

        public double meanInSegment(final int segment) {
            return get(MEANS_NAME, SegmentMeans.class).getMeanInSegment(segment);
        }
    }

    private final class CopyRatioDataCollection extends DataCollection {
        private static final String COVERAGES_NAME = "coverages";
        private static final String NUM_TARGETS_PER_SEGMENT_NAME = "numTargetsPerSegment";
        private static final String COVERAGES_IN_SEGMENT_PREFIX = "coveragesInSegment";

        private final int numSegments;

        public CopyRatioDataCollection(final File coveragesFile, final File numTargetsPerSegmentFile) {
            super();
            final Data<Double> coverages = new Data<>(COVERAGES_NAME, coveragesFile, Double::parseDouble);
            final Data<Integer> numTargetsPerSegment =
                    new Data<>(NUM_TARGETS_PER_SEGMENT_NAME, numTargetsPerSegmentFile, Integer::parseInt);
            numSegments = numTargetsPerSegment.size();
            final int numSegments = numTargetsPerSegment.size();
            int startTarget = 0;
            for (int segment = 0; segment < numSegments; segment++) {
                final int numTargetsInSegment = numTargetsPerSegment.values().get(segment);
                add(new Data<>(COVERAGES_IN_SEGMENT_PREFIX + segment,
                        coverages.values().subList(startTarget, startTarget + numTargetsInSegment)));
                startTarget += numTargetsInSegment;
            }
        }

        public List<Double> getCoveragesInSegment(final int segment) {
            return get(COVERAGES_IN_SEGMENT_PREFIX + segment);
        }

        public int getNumTargetsInSegment(final int segment) {
            return getCoveragesInSegment(segment).size();
        }
    }

    /**
     * We create a Modeller helper class to initialize the model state and specify the parameter samplers.
     */
    private final class CopyRatioModeller {
        //Create fields in the Modeller for the model and samplers.
        private final ParameterizedModel<CopyRatioState, CopyRatioDataCollection> model;
        private final Sampler<Double, CopyRatioState, CopyRatioDataCollection> varianceSampler;
        private final Sampler<SegmentMeans, CopyRatioState, CopyRatioDataCollection> meansSampler;  //note that this returns a SegmentMeans sample

        /**
         * Constructor for the Modeller takes as parameters all quantities needed to construct the ParameterizedState
         * (here, the initial variance and the initial mean) and the DataCollection
         * (here, the coverage and target-number files).
         */
        public CopyRatioModeller(final double initialVariance, final double initalMean,
                                 final File coveragesFile, final File numTargetsPerSegmentFile) {
            //TODO
            final CopyRatioDataCollection dataset = new CopyRatioDataCollection(coveragesFile, numTargetsPerSegmentFile);
            final CopyRatioState initialState = new CopyRatioState(initialVariance, initalMean, dataset.numSegments);

            //Implement Samplers for each parameter by overriding sample().  This can be done via a lambda that takes
            //(rng, state, dataCollection) and returns a new sample of the parameter with type identical to that
            //specified during initialization above.

            //Sampler for the variance global parameter.  The relevant logConditionalPDF is given by
            //the log of the product of Gaussian likelihoods for each target t:
            //      log[product_t variance^(-1/2) * exp(-(coverage_t - mean_t)^2 / (2 * variance))] + constant
            //which reduces to the form in code below.  Note that mean_t is identical for all targets in a segment.
            varianceSampler = (rng, state, data) -> {
                final Function<Double, Double> logConditionalPDF = newVariance -> {
                    double ll = 0.;
                    for (int segment = 0; segment < data.numSegments; segment++) {
                        final double meanInSegment = state.meanInSegment(segment);
                        ll += -0.5 * Math.log(newVariance) * data.getNumTargetsInSegment(segment) +
                                data.getCoveragesInSegment(segment).stream()
                                        .mapToDouble(c -> -normalTerm(c, meanInSegment, newVariance)).sum();
                    }
                    return ll;
                };

                final SliceSampler sampler =
                        new SliceSampler(rng, logConditionalPDF, state.variance(),
                                VARIANCE_MIN, VARIANCE_MAX, VARIANCE_WIDTH);
                return sampler.sample();
            };

            //Sampler for the segment-level mean parameters.  For each segment s, the relevant logConditionalPDF
            //is given by the log of the product of Gaussian likelihoods for all targets t in that segment:
            //     log[product_{t in s} exp(-(coverage_t - mean_s)^2 / (2 * variance))] + constant
            //which reduces to the form in code below.  Note that a lambda is specified for each segment,
            //so that meansSampler.sample() returns a SegmentMeans object that contains the means for all segments.
            //NOTE: It is possible to shoot yourself in the foot here if you incorrectly specify the sampler!
            //For example, note that there is no check that the number of segments sampled here matches the
            //number of segments given during initialization.  Such consistency checks can be added (e.g., in the
            //contructor of SegmentMeans), if desired, but may impact performance.
            meansSampler = (rng, state, data) -> {
                final List<Double> means = new ArrayList<>();
                for (int segment = 0; segment < data.numSegments; segment++) {
                    final List<Double> coveragesInSegment = data.getCoveragesInSegment(segment);
                    final Function<Double, Double> logConditionalPDF =
                            newMean -> coveragesInSegment.stream()
                                    .mapToDouble(c -> -normalTerm(c, newMean, state.variance()))
                                    .sum();

                    final SliceSampler sampler =
                            new SliceSampler(rng, logConditionalPDF, state.meanInSegment(segment), MEAN_WIDTH);
                    means.add(sampler.sample());
                }
                return new SegmentMeans(means);
            };

            //Build the ParameterizedModel using the GibbsBuilder pattern.
            //Pass in the initial CopyRatioState and CopyRatioDataCollection, and specify the class of the CopyRatioState.
            //Add samplers for each of the parameters, with names matching those used in initialization.
            model = new ParameterizedModel.GibbsBuilder<>(initialState, dataset, CopyRatioState.class)
                    .addParameterSampler(CopyRatioState.VARIANCE_NAME, varianceSampler)
                    .addParameterSampler(CopyRatioState.MEANS_NAME, meansSampler)
                    .build();
        }
    }

    /**
     * Tests Bayesian inference of a toy copy-ratio model via MCMC.
     * <p>
     *     Recovery of input values for the variance global parameter and the segment-level mean parameters is checked.
     *     In particular, the mean and standard deviation of the posterior for the variance must be recovered to within
     *     a relative error of 1% and 5%, respectively, in 250 samples (after 250 burn-in samples have been discarded).
     * </p>
     * <p>
     *     Furthermore, the number of truth values for the segment-level means falling outside confidence intervals of
     *     1-sigma, 2-sigma, and 3-sigma given by the posteriors in each segment should be roughly consistent with
     *     a normal distribution (i.e., ~32, ~5, and ~0, respectively; we allow for errors of 15, 6, and 2).
     *     Finally, the mean of the standard deviations of the posteriors for the segment-level means should be
     *     recovered to within a relative error of 5%.
     * </p>
     * <p>
     *     With these specifications, this unit test is not overly brittle (i.e., it should pass for a large majority
     *     of randomly generated data sets), but it is still brittle enough to check for correctness of the sampling
     *     (for example, specifying a sufficiently incorrect likelihood will cause the test to fail).
     * </p>
     */
    @Test
    public void testRunMCMCOnCopyRatioSegmentedModel() {
        //Create new instance of the Modeller helper class, passing all quantities needed to initialize state and data.
        final CopyRatioModeller modeller =
                new CopyRatioModeller(VARIANCE_INITIAL, MEAN_INITIAL, COVERAGES_FILE, NUM_TARGETS_PER_SEGMENT_FILE);
        //Create a GibbsSampler, passing the total number of samples (including burn-in samples)
        //and the model held by the Modeller.
        final GibbsSampler<CopyRatioState, CopyRatioDataCollection> gibbsSampler =
                new GibbsSampler<>(NUM_SAMPLES, modeller.model);
        //Run the MCMC.
        gibbsSampler.runMCMC();

        //Check that the statistics---i.e., the mean and standard deviation---of the variance posterior
        //agree with those found by emcee/analytically to a relative error of 1% and 10%, respectively.
        final double[] varianceSamples = Doubles.toArray(gibbsSampler.getSamples("variance", Double.class, NUM_BURN_IN));
        final double variancePosteriorMean = new Mean().evaluate(varianceSamples);
        final double variancePosteriorStandardDeviation = new StandardDeviation().evaluate(varianceSamples);
        Assert.assertEquals(relativeError(variancePosteriorMean, VARIANCE_TRUTH), 0., 0.01);
        Assert.assertEquals(
                relativeError(variancePosteriorStandardDeviation, VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH),
                0., 0.1);
        //Check statistics---i.e., the mean and standard deviation---of the segment-level mean posteriors.
        //In particular, check that the number of segments where the true mean falls outside confidence intervals
        //is roughly consistent with a normal distribution.
        final Data<Double> meansTruth = new Data<>("meansTruth", MEANS_TRUTH_FILE, Double::parseDouble);
        final int numSegments = meansTruth.size();
        final List<SegmentMeans> meansSamples = gibbsSampler.getSamples(CopyRatioState.MEANS_NAME, SegmentMeans.class, NUM_BURN_IN);
        int numMeansOutsideOneSigma = 0;
        int numMeansOutsideTwoSigma = 0;
        int numMeansOutsideThreeSigma = 0;
        final List<Double> meanPosteriorStandardDeviations = new ArrayList<>();
        for (int segment = 0; segment < numSegments; segment++) {
            final int j = segment;
            final double[] meanInSegmentSamples =
                    Doubles.toArray(meansSamples.stream()
                            .map(s -> s.get(SegmentMeans.MEAN_IN_SEGMENT_PREFIX + j, Double.class))
                            .collect(Collectors.toList()));
            final double meanPosteriorMean = new Mean().evaluate(meanInSegmentSamples);
            final double meanPosteriorStandardDeviation =
                    new StandardDeviation().evaluate(meanInSegmentSamples);
            meanPosteriorStandardDeviations.add(meanPosteriorStandardDeviation);
            final double absoluteDifferenceFromTruth = Math.abs(meanPosteriorMean - meansTruth.value(segment));
            if (absoluteDifferenceFromTruth > meanPosteriorStandardDeviation) {
                numMeansOutsideOneSigma++;
            }
            if (absoluteDifferenceFromTruth > 2 * meanPosteriorStandardDeviation) {
                numMeansOutsideTwoSigma++;
            }
            if (absoluteDifferenceFromTruth > 3 * meanPosteriorStandardDeviation) {
                numMeansOutsideThreeSigma++;
            }
        }
        final double meanPosteriorStandardDeviationsMean =
                new Mean().evaluate(Doubles.toArray(meanPosteriorStandardDeviations));
        Assert.assertEquals(numMeansOutsideOneSigma, 100 - 68, 15);
        Assert.assertEquals(numMeansOutsideTwoSigma, 100 - 95, 6);
        Assert.assertTrue(numMeansOutsideThreeSigma <= 2);
        Assert.assertEquals(
                relativeError(meanPosteriorStandardDeviationsMean, MEAN_POSTERIOR_STANDARD_DEVIATION_MEAN_TRUTH),
                0., 0.05);
    }

//    private static final CopyRatioData EMPTY_DATA =
//            new CopyRatioData(new Data<>(Collections.emptyList()), new Data<>(Collections.emptyList()));
//    private static final CopyRatioData SIMPLE_DATA =
//            new CopyRatioData(new Data<>(Collections.singletonList(1.)), new Data<>(Collections.singletonList(1)));
//
//    @Test(expectedExceptions = IllegalArgumentException.class)
//    public void testNegativeNumSamplesException() {
//        final GibbsSampler<CopyRatioData> gibbsSampler = new GibbsSampler.Builder<>(-1, EMPTY_DATA).build();
//    }
//
//    @Test(expectedExceptions = IllegalArgumentException.class)
//    public void testNoSamplersException() {
//        final GibbsSampler<CopyRatioData> gibbsSampler =
//                new GibbsSampler.Builder<>(1, EMPTY_DATA)
//                        .build();
//    }
//
//    @Test(expectedExceptions = IllegalArgumentException.class)
//    public void testNumberOfGlobalParametersAndSamplersMismatchException1() {
//        final GibbsSampler<CopyRatioData> gibbsSampler =
//                new GibbsSampler.Builder<>(1, EMPTY_DATA)
//                        .addGlobalSampler(varianceSampler)
//                        .build();
//    }
//
//    @Test(expectedExceptions = IllegalArgumentException.class)
//    public void testNumberOfGlobalParametersAndSamplersMismatchException2() {
//        final GibbsSampler<CopyRatioData> gibbsSampler =
//                new GibbsSampler.Builder<>(1, EMPTY_DATA)
//                        .initializeGlobalParameters(Arrays.asList(1., 1.))
//                        .addGlobalSampler(varianceSampler)
//                        .build();
//    }
//
//    @Test(expectedExceptions = IllegalArgumentException.class)
//    public void testNumberOfLocalParametersAndSamplersMismatchException1() {
//        final GibbsSampler<CopyRatioData> gibbsSampler =
//                new GibbsSampler.Builder<>(1, EMPTY_DATA)
//                        .addLocalSampler(meansSampler)
//                        .build();
//    }
//
//    @Test(expectedExceptions = IllegalArgumentException.class)
//    public void testNumberOfLocalParametersAndSamplersMismatchException2() {
//        final GibbsSampler<CopyRatioData> gibbsSampler =
//                new GibbsSampler.Builder<>(1, EMPTY_DATA)
//                        .initializeLocalParameters(Arrays.asList(Collections.singletonList(1.), Collections.singletonList(1.)))
//                        .addLocalSampler(meansSampler)
//                        .build();
//    }
//
//    @Test(expectedExceptions = IllegalArgumentException.class)
//    public void testNumBurnInException1() {
//        final GibbsSampler<CopyRatioData> gibbsSampler =
//                new GibbsSampler.Builder<>(10, SIMPLE_DATA)
//                        .initializeGlobalParameters(Collections.singletonList(1.))
//                        .initializeLocalParameters(Collections.singletonList(Collections.singletonList(1.)))
//                        .addGlobalSampler(varianceSampler)
//                        .addLocalSampler(meansSampler)
//                        .build();
//        gibbsSampler.runMCMC();
//        gibbsSampler.getGlobalSamples(10);
//    }
//
//    @Test(expectedExceptions = IllegalArgumentException.class)
//    public void testNumBurnInException2() {
//        final GibbsSampler<CopyRatioData> gibbsSampler =
//                new GibbsSampler.Builder<>(10, SIMPLE_DATA)
//                        .initializeGlobalParameters(Collections.singletonList(1.))
//                        .initializeLocalParameters(Collections.singletonList(Collections.singletonList(1.)))
//                        .addGlobalSampler(varianceSampler)
//                        .addLocalSampler(meansSampler)
//                        .build();
//        gibbsSampler.runMCMC();
//        gibbsSampler.getLocalSamples(-10);
//    }
//
//    @Test(expectedExceptions = IllegalStateException.class)
//    public void testLocalParameterAndSamplerMultiplicityMismatch() {
//        final GibbsSampler<CopyRatioData> gibbsSampler =
//                new GibbsSampler.Builder<>(10, SIMPLE_DATA)
//                        .initializeGlobalParameters(Collections.singletonList(1.))
//                        .initializeLocalParameters(Collections.singletonList(Arrays.asList(1., 1.)))
//                        .addGlobalSampler(varianceSampler)
//                        .addLocalSampler(meansSampler)
//                        .build();
//        gibbsSampler.runMCMC();
//    }
//
//    @Test
//    public void testGetGlobalSamplesWithoutRunMCMC() {
//        final GibbsSampler<CopyRatioData> gibbsSampler =
//                new GibbsSampler.Builder<>(10, SIMPLE_DATA)
//                        .initializeGlobalParameters(Collections.singletonList(1.))
//                        .initializeLocalParameters(Collections.singletonList(Collections.singletonList(1.)))
//                        .addGlobalSampler(varianceSampler)
//                        .addLocalSampler(meansSampler)
//                        .build();
//        final List<List<Double>> globalSamples = gibbsSampler.getGlobalSamples(0);
//        Assert.assertEquals(globalSamples.size(), 10);
//    }
//
//    @Test
//    public void testGetLocalSamplesWithoutRunMCMC() {
//        final GibbsSampler<CopyRatioData> gibbsSampler =
//                new GibbsSampler.Builder<>(10, SIMPLE_DATA)
//                        .initializeGlobalParameters(Collections.singletonList(1.))
//                        .initializeLocalParameters(Collections.singletonList(Collections.singletonList(1.)))
//                        .addGlobalSampler(varianceSampler)
//                        .addLocalSampler(meansSampler)
//                        .build();
//        final List<List<List<Double>>> localSamples = gibbsSampler.getLocalSamples(0);
//        Assert.assertEquals(localSamples.size(), 10);
//    }
}