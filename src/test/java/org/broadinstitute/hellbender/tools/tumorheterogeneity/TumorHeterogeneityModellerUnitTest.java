package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Doubles;
import htsjdk.samtools.util.Log;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class TumorHeterogeneityModellerUnitTest extends BaseTest {
    private static final int NUM_POINTS = 500;
    private static final double VARIANCE_TRUTH = 0.04;
    private static final int NUM_POPULATIONS_TRUTH = 5;
    private static final List<Double> MEANS_TRUTH = Arrays.asList(-5., -3., 1., 3., 5.);
    private static final List<Double> POPULATION_FRACTIONS_TRUTH = Arrays.asList(0.1, 0.2, 0.5, 0.1, 0.1);

//    private static final int NUM_POPULATIONS_TRUTH = 2;
//    private static final List<Double> MEANS_TRUTH = Arrays.asList(-1., 1.);
//    private static final List<Double> POPULATION_FRACTIONS_TRUTH = Arrays.asList(0.7, 0.3);

    private static final int RANDOM_SEED = 13;
    private static final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));
    
    private static final double CREDIBLE_INTERVAL_ALPHA = 0.1;

    @Test
    public void testRunMCMC() throws IOException {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);

        rng.setSeed(RANDOM_SEED);

        final List<Integer> populationIndices = IntStream.range(0, NUM_POPULATIONS_TRUTH).boxed().collect(Collectors.toList());
        final List<Integer> randomPopulationIndices = IntStream.range(0, NUM_POINTS).boxed()
                .map(i -> GATKProtectedMathUtils.randomSelect(populationIndices, POPULATION_FRACTIONS_TRUTH::get, rng)).collect(Collectors.toList());
        final List<Double> points = randomPopulationIndices.stream().map(i -> new NormalDistribution(rng, MEANS_TRUTH.get(i), Math.sqrt(VARIANCE_TRUTH)).sample()).collect(Collectors.toList());
        ParamUtils.writeValuesToFile(Doubles.toArray(points), new File("points.txt"));

        final int numPopulations = 10;
        final int numSamples = 2000;
        final int numBurnIn = 1000;

        //run MCMC
        final TumorHeterogeneityModeller modeller = new TumorHeterogeneityModeller(points, numPopulations, rng);
        modeller.fitMCMC(numSamples, numBurnIn);
        
        //check statistics of global-parameter posterior samples (i.e., posterior mode and standard deviation)
        final Map<TumorHeterogeneityParameter, PosteriorSummary> globalParameterPosteriorSummaries =
                modeller.getGlobalParameterPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);

        final PosteriorSummary concentrationPosteriorSummary = globalParameterPosteriorSummaries.get(TumorHeterogeneityParameter.CONCENTRATION);
        final double concentrationPosteriorCenter = concentrationPosteriorSummary.getCenter();
        final double concentrationPosteriorStandardDeviation = (concentrationPosteriorSummary.getUpper() - concentrationPosteriorSummary.getLower()) / 2;
        System.out.println("concentration: " + concentrationPosteriorCenter + " " + concentrationPosteriorStandardDeviation);
        System.out.println();

        final PosteriorSummary variancePosteriorSummary = globalParameterPosteriorSummaries.get(TumorHeterogeneityParameter.VARIANCE);
        final double variancePosteriorCenter = variancePosteriorSummary.getCenter();
        final double variancePosteriorStandardDeviation = (variancePosteriorSummary.getUpper() - variancePosteriorSummary.getLower()) / 2;
        System.out.println("variance: " + variancePosteriorCenter + " " + variancePosteriorStandardDeviation);
        ParamUtils.writeValuesToFile(Doubles.toArray(modeller.getVarianceSamples()), new File("variance.txt"));
        System.out.println();
//        Assert.assertEquals(Math.abs(variancePosteriorCenter - VARIANCE_TRUTH),
//                0., MULTIPLES_OF_SD_THRESHOLD * VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH);
//        Assert.assertEquals(relativeError(variancePosteriorStandardDeviation, VARIANCE_POSTERIOR_STANDARD_DEVIATION_TRUTH),
//                0., RELATIVE_ERROR_THRESHOLD);

        for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
            final int j = populationIndex;
            final double[] meanSamples = Doubles.toArray(modeller.getMeansSamples().stream().map(s -> s.get(j)).collect(Collectors.toList()));
            final double meanMean = new Mean().evaluate(meanSamples);
            final double meanStandardDeviation = new StandardDeviation().evaluate(meanSamples);
            System.out.println("mean " + populationIndex + ": " + meanMean + " " + meanStandardDeviation);
            ParamUtils.writeValuesToFile(meanSamples, new File("mean-" + j + ".txt"));

            final double[] populationFractionSamples = Doubles.toArray(modeller.getPopulationFractionsSamples().stream().map(s -> s.get(j)).collect(Collectors.toList()));
            final double populationFractionMean = new Mean().evaluate(populationFractionSamples);
            final double populationFractionStandardDeviation = new StandardDeviation().evaluate(populationFractionSamples);
            System.out.println("fraction " + populationIndex + ": " + populationFractionMean + " " + populationFractionStandardDeviation);
            ParamUtils.writeValuesToFile(populationFractionSamples, new File("fraction-" + j + ".txt"));

            System.out.println();
        }

//        //check statistics of segment-mean posterior samples (i.e., posterior means and standard deviations)
//        int numMeansOutsideOneSigma = 0;
//        int numMeansOutsideTwoSigma = 0;
//        int numMeansOutsideThreeSigma = 0;
//        final List<PosteriorSummary> meanPosteriorSummaries =
//                modeller.getSegmentMeansPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);
//        final double[] meanPosteriorStandardDeviations = new double[numSegments];
//        for (int segment = 0; segment < numSegments; segment++) {
//            final double meanPosteriorCenter = meanPosteriorSummaries.get(segment).getCenter();
//            final double meanPosteriorStandardDeviation =
//                    (meanPosteriorSummaries.get(segment).getUpper() - meanPosteriorSummaries.get(segment).getLower()) / 2.;
//            meanPosteriorStandardDeviations[segment] = meanPosteriorStandardDeviation;
//            final double absoluteDifferenceFromTruth = Math.abs(meanPosteriorCenter - meansTruth.get(segment));
//            if (absoluteDifferenceFromTruth > meanPosteriorStandardDeviation) {
//                numMeansOutsideOneSigma++;
//            }
//            if (absoluteDifferenceFromTruth > 2 * meanPosteriorStandardDeviation) {
//                numMeansOutsideTwoSigma++;
//            }
//            if (absoluteDifferenceFromTruth > 3 * meanPosteriorStandardDeviation) {
//                numMeansOutsideThreeSigma++;
//            }
//        }
//        final double meanPosteriorStandardDeviationsMean = new Mean().evaluate(meanPosteriorStandardDeviations);
//        Assert.assertEquals(numMeansOutsideOneSigma, 100 - 68, DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_1_SIGMA);
//        Assert.assertEquals(numMeansOutsideTwoSigma, 100 - 95, DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_2_SIGMA);
//        Assert.assertTrue(numMeansOutsideThreeSigma <= DELTA_NUMBER_OF_MEANS_ALLOWED_OUTSIDE_3_SIGMA);
//        Assert.assertEquals(relativeError(meanPosteriorStandardDeviationsMean, MEAN_POSTERIOR_STANDARD_DEVIATION_MEAN_TRUTH),
//                0., RELATIVE_ERROR_THRESHOLD);
//
//        //check accuracy of latent outlier-indicator posterior samples
//        final List<TumorHeterogeneityState.OutlierIndicators> outlierIndicatorSamples =
//                modeller.getOutlierIndicatorsSamples();
//        int numIndicatorsCorrect = 0;
//        final int numIndicatorSamples = outlierIndicatorSamples.size();
//        final List<Integer> outlierIndicatorsTruthAsInt = loadList(OUTLIER_INDICATORS_TRUTH_FILE, Integer::parseInt);
//        final List<Boolean> outlierIndicatorsTruth =
//                outlierIndicatorsTruthAsInt.stream().map(i -> i == 1).collect(Collectors.toList());
//        for (int target = 0; target < coverage.targets().size(); target++) {
//            int numSamplesOutliers = 0;
//            for (final TumorHeterogeneityState.OutlierIndicators sample : outlierIndicatorSamples) {
//                if (sample.get(target)) {
//                    numSamplesOutliers++;
//                }
//            }
//            //take predicted state of indicator to be given by the majority of samples
//            if ((numSamplesOutliers >= numIndicatorSamples / 2.) == outlierIndicatorsTruth.get(target)) {
//                numIndicatorsCorrect++;
//            }
//        }
//        final double fractionOfOutlierIndicatorsCorrect = (double) numIndicatorsCorrect / coverage.targets().size();
//        Assert.assertTrue(fractionOfOutlierIndicatorsCorrect >= FRACTION_OF_OUTLIER_INDICATORS_CORRECT_THRESHOLD);
    }
}