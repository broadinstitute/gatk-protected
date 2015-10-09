package org.broadinstitute.hellbender.utils.mcmc;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link GibbsSampler}.
 * <p>
 *     Test case is a toy copy-ratio model with 1 global and 1 local (segment-level) parameter.
 *     We consider a case with 100 segments, each of which has a number of targets that is drawn from a
 *     Poisson distribution with mean 100 in the generated data set.
 *     Data consists of a list of coverages at each site and a list of the number of targets in each segment.
 *     Coverages are assumed to be drawn from a normal distribution in each segment.
 *     The mean of each normal distribution is given by a segment-level parameter;
 *     its distribution across segments is taken to be uniform in (0, 10) in the generated data.
 *     The variance of each distribution is set by the global parameter, which is taken to be 1 in the generated data.
 *     Success of the test is determined by recovery of the input segment means and the global variance,
 *     as well as agreement of the standard deviations of the parameter posteriors with those given by both the
 *     python package emcee (see http://dan.iel.fm/emcee for details) and numerical evaluation in Mathematica
 *     of the analytic forms of the posteriors.
 * </p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class GibbsSamplerUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/utils/mcmc/";

    private static final File COVERAGES_FILE = new File(TEST_SUB_DIR
            + "coverages-for-gibbs-sampler-copy-ratio-test.txt");
    private static final File NUM_TARGETS_PER_SEGMENT_FILE =
            new File(TEST_SUB_DIR + "number-of-targets-per-segment-for-gibbs-sampler-copy-ratio-test.txt");
    private static final File MEANS_TRUTH_FILE = new File(TEST_SUB_DIR
            + "means-truth-for-gibbs-sampler-copy-ratio-test.txt");

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
    private static final double MEAN_POSTERIOR_STANDARD_DEVIATION_MEAN_TRUTH = 0.1;
    
    private static final int NUM_SAMPLES = 500;
    private static final int NUM_BURN_IN = 250;

    private static List<Double> loadList(final File file) {
        final List<Double> list = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = br.readLine()) != null) {
                list.add(Double.parseDouble(line));
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(file, e);
        }
        return list;
    }

    private static List<List<Double>> loadData(final File coveragesFile, final File numTargetsPerSegmentFile) {
        final List<List<Double>> data = new ArrayList<>();
        data.add(loadList(coveragesFile));
        data.add(loadList(numTargetsPerSegmentFile));
        return data;
    }

    private static double normalTerm(final double quantity, final double mean, final double variance) {
        return (quantity - mean) * (quantity - mean) / (2 * variance);
    }

    private static double relativeError(final double x, final double y) {
        return Math.abs(x - y) / (2 * (x + y));
    }

    private final class CopyRatioModelState {
        //common variance of normal-distribution components in mixture model fit to log_2 coverages in each segment
        private final double variance;
        //contribution of uniform component in mixture model fit to log_2 coverages
        private final List<Double> means;

        private CopyRatioModelState(final double variance, final List<Double> means) {
            this.variance = variance;
            this.means = means;
        }
        private CopyRatioModelState(final List<Double> globalParameters, final List<List<Double>> localParameters) {
            this(globalParameters.get(0), localParameters.get(0));
        }
    }

    private final class CopyRatioData {
        private final List<Double> coverages;
        private final List<Integer> numTargetsPerSegment;

        private CopyRatioData(final List<List<Double>> dataAsDouble) {
            coverages = dataAsDouble.get(0);
            numTargetsPerSegment = dataAsDouble.get(1).stream().map(Double::intValue).collect(Collectors.toList());
        }
    }

    /**
     * Sampler for the variance global parameter.  The log of its conditional probability density function is
     * specified by the lambda logConditionalPDF and is given by the log of the product of Gaussian likelihoods
     * for each target t; up to an additive constant this is:
     * <p>
     *     log[product_t variance^(-1/2) * exp(-(coverage_t - mean_t)^2 / (2 * variance))]
     * </p>
     * which reduces to the form in code below.  Note that mean_t is identical for all targets in a segment.
     */
    private final GlobalSampler varianceSampler = (rng, globalParameters, localParameters, dataAsDouble) -> {
        final CopyRatioModelState model = new CopyRatioModelState(globalParameters, localParameters);
        final CopyRatioData data = new CopyRatioData(dataAsDouble);

        final Function<Double, Double> logConditionalPDF =
                newVariance -> {
                    double ll = 0.;
                    int targetStart = 0;
                    for (int segment = 0; segment < data.numTargetsPerSegment.size(); segment++) {
                        final double mean = model.means.get(segment);
                        final List<Double> coveragesInSegment =
                                data.coverages.subList(targetStart, targetStart + data.numTargetsPerSegment.get(segment));
                        final int numTargetsInSegment = data.numTargetsPerSegment.get(segment);
                        ll += -0.5 * Math.log(newVariance) * numTargetsInSegment +
                                coveragesInSegment.stream()
                                .mapToDouble(c -> -normalTerm(c, mean, newVariance))
                                .sum();
                        targetStart += numTargetsInSegment;
                    }
                    return ll;
                };

        final SliceSampler sampler =
                new SliceSampler(rng, logConditionalPDF, model.variance, VARIANCE_MIN, VARIANCE_MAX, VARIANCE_WIDTH);
        return sampler.sample();
    };

    /**
     * Sampler for the segment-level mean local parameter.  For a given segment s, the log of its
     * conditional probability density function is specified by the lambda logConditionalPDF and is given by the
     * log of the product of Gaussian likelihoods for all targets t in that segment; up to an additive constant, this is:
     * <p>
     *     log[product_{t in s} exp(-(coverage_t - mean_s)^2 / (2 * variance))]
     * </p>
     * which reduces to the form in code below.  Note that a lambda is specified for each segment,
     * so that meansSampler.sample() returns a sample that is a double[] of length given by the number of segments.
     */
    private final LocalSampler meansSampler = (rng, globalParameters, localParameters, dataAsDouble) -> {
        final CopyRatioModelState model = new CopyRatioModelState(globalParameters, localParameters);
        final CopyRatioData data = new CopyRatioData(dataAsDouble);

        final List<Double> sample = new ArrayList<>(model.means.size());
        int startTarget = 0;
        for (int segment = 0; segment < data.numTargetsPerSegment.size(); segment++) {
            final List<Double> coveragesInSegment =
                    data.coverages.subList(startTarget, startTarget + data.numTargetsPerSegment.get(segment));
            final int numTargetsInSegment = data.numTargetsPerSegment.get(segment);

            final Function<Double, Double> logConditionalPDF =
                    newMean -> coveragesInSegment.stream()
                            .mapToDouble(c -> -normalTerm(c, newMean, model.variance))
                            .sum();

            final SliceSampler sampler =
                    new SliceSampler(rng, logConditionalPDF, model.means.get(segment), MEAN_MIN, MEAN_MAX, MEAN_WIDTH);
            sample.add(sampler.sample());

            startTarget += numTargetsInSegment;
        }
        return sample;
    };

    /**
     * Tests the MCMC of a toy copy-ratio model.  Recovery of input values for the variance global parameter
     * and the segment-level mean local parameter are checked.  In particular, the mean and standard deviation
     * of the posterior for the variance must be recovered to within a relative error of 1% and 5%, respectively.
     * Furthermore, the number of truth values for the segment-level means falling outside confidence intervals of
     * 1-sigma, 2-sigma, and 3-sigma given by the posteriors in each segment should be roughly consistent with
     * a normal distribution (i.e., ~32, ~5, and ~0, respectively; we allow for errors of 15, 6, and 2).
     * Finally, the mean of the standard deviations of the posteriors for the segment-level means should be
     * recovered to within a relative error of 5%.  With these specifications, this unit test is not overly
     * brittle (i.e., it should pass for a large majority of randomly generated data sets), but it is still
     * brittle enough to check for correctness of the sampling (for example, specifying a sufficiently
     * incorrect likelihood will cause the test to fail).
     */
    @Test
    public void testRunMCMCOnCopyRatioSegmentedModel() {
        //load data (coverages and number of targets in each segment)
        final List<List<Double>> data = loadData(COVERAGES_FILE, NUM_TARGETS_PER_SEGMENT_FILE);
        final int numSegments = data.get(1).size();
        //initialize parameters
        final List<Double> initialGlobalParameters = Arrays.asList(VARIANCE_INITIAL);
        final List<List<Double>> initialLocalParameters =
                Arrays.asList(
                        IntStream.range(0, numSegments).mapToDouble(i -> MEAN_INITIAL).boxed()
                                .collect(Collectors.toList())
                );
        //build GibbsSampler and run MCMC
        final GibbsSampler gibbsSampler = new GibbsSampler.Builder(NUM_SAMPLES, data)
                .initializeGlobalParameters(initialGlobalParameters)
                .initializeLocalParameters(initialLocalParameters)
                .addGlobalSampler(varianceSampler)
                .addLocalSampler(meansSampler)
                .build();
        gibbsSampler.runMCMC();
        //get samples
        final List<List<Double>> globalSamples = gibbsSampler.getGlobalSamples(NUM_BURN_IN);
        final List<List<List<Double>>> localSamples = gibbsSampler.getLocalSamples(NUM_BURN_IN);
        //check statistics of global-parameter posterior samples (i.e., posterior mean and standard deviation)
        final List<Double> varianceSamples = IntStream.range(0, NUM_SAMPLES - NUM_BURN_IN).boxed()
                    .map(i -> globalSamples.get(i).get(0)).collect(Collectors.toList());
        final double variancePosteriorMean = new Mean().evaluate(Doubles.toArray(varianceSamples));
        final double variancePosteriorStandardDeviation =
                new StandardDeviation().evaluate(Doubles.toArray(varianceSamples));
        Assert.assertEquals(relativeError(variancePosteriorMean, VARIANCE_TRUTH), 0., 0.01);
        Assert.assertEquals(
                relativeError(variancePosteriorStandardDeviation, VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH),
                0., 0.05);
        //check statistics of local-parameter posterior samples (i.e., posterior means and standard deviations)
        final List<Double> meansTruth = loadList(MEANS_TRUTH_FILE);
        int numMeansOutsideOneSigma = 0;
        int numMeansOutsideTwoSigma = 0;
        int numMeansOutsideThreeSigma = 0;
        final List<Double> meanPosteriorStandardDeviations = new ArrayList<>();
        for (int segment = 0; segment < numSegments; segment++) {
            final int j = segment;
            final List<Double> meanSamples = IntStream.range(0, NUM_SAMPLES - NUM_BURN_IN)
                    .mapToDouble(i -> localSamples.get(i).get(0).get(j)).boxed().collect(Collectors.toList());
            final double meanPosteriorMean = new Mean().evaluate(Doubles.toArray(meanSamples));
            final double meanPosteriorStandardDeviation =
                    new StandardDeviation().evaluate(Doubles.toArray(meanSamples));
            meanPosteriorStandardDeviations.add(meanPosteriorStandardDeviation);
            final double absoluteDifferenceFromTruth = Math.abs(meanPosteriorMean - meansTruth.get(segment));
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

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNegativeNumSamplesException() {
        final GibbsSampler gibbsSampler = new GibbsSampler.Builder(-1, Arrays.asList()).build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNoSamplersException() {
        final GibbsSampler gibbsSampler =
                new GibbsSampler.Builder(1, Arrays.asList())
                        .build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNumberOfGlobalParametersAndSamplersMismatchException1() {
        final GibbsSampler gibbsSampler =
                new GibbsSampler.Builder(1, Arrays.asList())
                        .addGlobalSampler(varianceSampler)
                        .build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNumberOfGlobalParametersAndSamplersMismatchException2() {
        final GibbsSampler gibbsSampler =
                new GibbsSampler.Builder(1, Arrays.asList())
                        .initializeGlobalParameters(Arrays.asList(1., 1.))
                        .addGlobalSampler(varianceSampler)
                        .build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNumberOfLocalParametersAndSamplersMismatchException1() {
        final GibbsSampler gibbsSampler =
                new GibbsSampler.Builder(1, Arrays.asList())
                        .addLocalSampler(meansSampler)
                        .build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNumberOfLocalParametersAndSamplersMismatchException2() {
        final GibbsSampler gibbsSampler =
                new GibbsSampler.Builder(1, Arrays.asList())
                        .initializeLocalParameters(Arrays.asList(Arrays.asList(1.), Arrays.asList(1.)))
                        .addLocalSampler(meansSampler)
                        .build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNumBurnInException1() {
        final GibbsSampler gibbsSampler =
                new GibbsSampler.Builder(10, Arrays.asList(Arrays.asList(1.), Arrays.asList(1.)))
                        .initializeGlobalParameters(Arrays.asList(1.))
                        .initializeLocalParameters(Arrays.asList(Arrays.asList(1.)))
                        .addGlobalSampler(varianceSampler)
                        .addLocalSampler(meansSampler)
                        .build();
        gibbsSampler.runMCMC();
        gibbsSampler.getGlobalSamples(10);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNumBurnInException2() {
        final GibbsSampler gibbsSampler =
                new GibbsSampler.Builder(10, Arrays.asList(Arrays.asList(1.), Arrays.asList(1.)))
                        .initializeGlobalParameters(Arrays.asList(1.))
                        .initializeLocalParameters(Arrays.asList(Arrays.asList(1.)))
                        .addGlobalSampler(varianceSampler)
                        .addLocalSampler(meansSampler)
                        .build();
        gibbsSampler.runMCMC();
        gibbsSampler.getLocalSamples(-10);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testLocalParameterAndSamplerMultiplicityMismatch() {
        final GibbsSampler gibbsSampler =
                new GibbsSampler.Builder(10, Arrays.asList(Arrays.asList(1.), Arrays.asList(1.)))
                        .initializeGlobalParameters(Arrays.asList(1.))
                        .initializeLocalParameters(Arrays.asList(Arrays.asList(1., 1.)))
                        .addGlobalSampler(varianceSampler)
                        .addLocalSampler(meansSampler)
                        .build();
        gibbsSampler.runMCMC();
    }

    @Test
    public void testGetGlobalSamplesWithoutRunMCMC() {
        final GibbsSampler gibbsSampler =
                new GibbsSampler.Builder(10, Arrays.asList(Arrays.asList(1.), Arrays.asList(1.)))
                        .initializeGlobalParameters(Arrays.asList(1.))
                        .initializeLocalParameters(Arrays.asList(Arrays.asList(1.)))
                        .addGlobalSampler(varianceSampler)
                        .addLocalSampler(meansSampler)
                        .build();
        final List<List<Double>> globalSamples = gibbsSampler.getGlobalSamples(0);
        Assert.assertEquals(globalSamples.size(), 10);
    }

    @Test
    public void testGetLocalSamplesWithoutRunMCMC() {
        final GibbsSampler gibbsSampler =
                new GibbsSampler.Builder(10, Arrays.asList(Arrays.asList(1.), Arrays.asList(1.)))
                        .initializeGlobalParameters(Arrays.asList(1.))
                        .initializeLocalParameters(Arrays.asList(Arrays.asList(1.)))
                        .addGlobalSampler(varianceSampler)
                        .addLocalSampler(meansSampler)
                        .build();
        final List<List<List<Double>>> localSamples = gibbsSampler.getLocalSamples(0);
        Assert.assertEquals(localSamples.size(), 10);
    }
}