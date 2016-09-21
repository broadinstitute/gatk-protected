package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import htsjdk.samtools.util.Log;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.function.Function;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class TumorHeterogeneityModellerUnitTest extends BaseTest {
    private static final int RANDOM_SEED = 13;
    private static final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

    private static final double CREDIBLE_INTERVAL_ALPHA = 0.95;

    private static final String TEST_DIR = "src/test/resources/org/broadinstitute/hellbender/tools/exome/cnlohcaller/";
    private static final File ACNV_SEG_FILE = new File(TEST_DIR, "cell_line-sim-final.seg");
    
    @Test
    public void testRunMCMC() throws IOException {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        LoggingUtils.setLoggingLevel(Log.LogLevel.DEBUG);

        rng.setSeed(RANDOM_SEED);

        final List<ACNVModeledSegment> segments = SegmentUtils.readACNVModeledSegmentFile(ACNV_SEG_FILE);

        final PloidyState normalPloidyState = new PloidyState(1, 1);
        final Function<PloidyState, Double> ploidyPDF = ps -> Math.log(Math.pow(0.75, (ps.m() == 0 ? 1 : 0) + (ps.n() == 0 ? 1 : 0)) / Math.pow(Math.abs(normalPloidyState.m() - ps.m()) + Math.abs(normalPloidyState.n() - ps.n()), 3));
        final Map<PloidyState, Double> unnormalizedLogProbabilityMassFunctionMap = new LinkedHashMap<>();
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 0), ploidyPDF.apply(new PloidyState(0, 0)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 1), ploidyPDF.apply(new PloidyState(0, 1)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 2), ploidyPDF.apply(new PloidyState(0, 2)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 3), ploidyPDF.apply(new PloidyState(0, 3)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(1, 2), ploidyPDF.apply(new PloidyState(1, 2)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 4), ploidyPDF.apply(new PloidyState(0, 4)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(1, 3), ploidyPDF.apply(new PloidyState(1, 3)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(2, 2), ploidyPDF.apply(new PloidyState(2, 2)));
        final PloidyStatePrior variantPloidyStatePrior = new PloidyStatePrior(unnormalizedLogProbabilityMassFunctionMap);

        final int numPopulations = 5;
        final int numCells = 10;

        final int numSamples = 100;
        final int numBurnIn = 20;

        final double concentrationPriorAlpha = 1.;
        final double concentrationPriorBeta = 100.;
        final double variantSegmentFractionPriorAlpha = 3.;
        final double variantSegmentFractionPriorBeta = 10.;

        //run MCMC
        final TumorHeterogeneityModeller modeller = new TumorHeterogeneityModeller(
                segments, normalPloidyState, variantPloidyStatePrior,
                concentrationPriorAlpha, concentrationPriorBeta, variantSegmentFractionPriorAlpha, variantSegmentFractionPriorBeta,
                numPopulations, numCells, rng);
        modeller.fitMCMC(numSamples, numBurnIn);

        //check statistics of global-parameter posterior samples (i.e., posterior mode and standard deviation)
        final Map<TumorHeterogeneityParameter, PosteriorSummary> globalParameterPosteriorSummaries =
                modeller.getGlobalParameterPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);

        final PosteriorSummary concentrationPosteriorSummary = globalParameterPosteriorSummaries.get(TumorHeterogeneityParameter.CONCENTRATION);
        final double concentrationPosteriorCenter = concentrationPosteriorSummary.getCenter();
        final double concentrationPosteriorStandardDeviation = (concentrationPosteriorSummary.getUpper() - concentrationPosteriorSummary.getLower()) / 2;
        System.out.println("concentration: " + concentrationPosteriorCenter + " " + concentrationPosteriorStandardDeviation);
        System.out.println();

        final List<TumorHeterogeneityState.PopulationFractions> populationFractionsSamples = modeller.getPopulationFractionsSamples();
        for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
            final int pi = populationIndex;
            final double[] populationFractionSamples = populationFractionsSamples.stream().mapToDouble(s -> s.get(pi)).toArray();
            final double populationFractionPosteriorMean = new Mean().evaluate(populationFractionSamples);
            final double populationFractionPosteriorStandardDeviation = new StandardDeviation().evaluate(populationFractionSamples);
            System.out.println("population fraction " + populationIndex + ": " + populationFractionPosteriorMean + " " + populationFractionPosteriorStandardDeviation);
        }
    }
}