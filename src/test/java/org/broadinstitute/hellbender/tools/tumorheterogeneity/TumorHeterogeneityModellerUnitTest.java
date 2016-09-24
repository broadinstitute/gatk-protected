package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import htsjdk.samtools.util.Log;
import org.apache.commons.collections4.MultiSet;
import org.apache.commons.collections4.multiset.HashMultiSet;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class TumorHeterogeneityModellerUnitTest extends BaseTest {
    private static final int RANDOM_SEED = 13;
    private static final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

    private static final double CREDIBLE_INTERVAL_ALPHA = 0.95;

//    private static final File ACNV_SEG_FILE = new File("/home/slee/working/ipython/purity-ploidy/clonal_test_data/seed-1_trunc-frac-1.0_segments-1000_length-20/purity-0.7_total_segments.acnv.seg");
    private static final File ACNV_SEG_FILE = new File("/home/slee/working/ipython/purity-ploidy/2_clones_test_data/seed-5_trunc-frac-0.5_segments-1000_length-20/purity-0.6_total_segments.acnv.seg");
//    private static final File ACNV_SEG_FILE = new File("/home/slee/working/ipython/purity-ploidy/purity-series/SM-74P4M-sim-final-edit.seg");
//    private static final File ACNV_SEG_FILE = new File("/home/slee/working/ipython/purity-ploidy/purity-series/SM-74P3M-sim-final-edit.seg");
    
    @Test
    public void testRunMCMC() throws IOException {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);

        rng.setSeed(RANDOM_SEED);

        final List<ACNVModeledSegment> segments = SegmentUtils.readACNVModeledSegmentFile(ACNV_SEG_FILE);//.stream().filter(s -> s.getInterval().size() > 100000).collect(Collectors.toList());

        final File mafCrFile = new File("/home/slee/working/ipython/purity-ploidy/maf-cr/maf-cr.tsv");
        mafCrFile.createNewFile();
        final TumorHeterogeneityData data = new TumorHeterogeneityData(segments);
        try (final FileWriter writer = new FileWriter(mafCrFile)) {
            for (double f = 0.; f <= 0.5; f += 0.01) {
                for (double c = 1E-10; c <= 5; c += 0.05) {
                    final double maf = f;
                    final double cr = c;
                    final double density = IntStream.range(0, data.numSegments()).mapToDouble(i -> Math.exp(data.logDensity(i, cr, maf))).sum();
                    writer.write(maf + "\t" + cr + "\t" + density + "\n");
                }
            }
        } catch (final IOException ioe) {
            throw new GATKException("Cannot write to file.", ioe);
        }

        System.out.println();

        final PloidyState normalPloidyState = new PloidyState(1, 1);
        final Function<PloidyState, Double> ploidyPDF = ps -> Math.log(Math.pow(0.75, (ps.m() == 0 ? 1 : 0) + (ps.n() == 0 ? 1 : 0)) / Math.pow(Math.abs(normalPloidyState.m() - ps.m()) + Math.abs(normalPloidyState.n() - ps.n()), 3));
//        final Function<PloidyState, Double> ploidyPDF = ps -> 0.;
        final Map<PloidyState, Double> unnormalizedLogProbabilityMassFunctionMap = new LinkedHashMap<>();
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 0), ploidyPDF.apply(new PloidyState(0, 0)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 1), ploidyPDF.apply(new PloidyState(0, 1)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 2), ploidyPDF.apply(new PloidyState(0, 2)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 3), ploidyPDF.apply(new PloidyState(0, 3)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(1, 2), ploidyPDF.apply(new PloidyState(1, 2)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 4), ploidyPDF.apply(new PloidyState(0, 4)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(1, 3), ploidyPDF.apply(new PloidyState(1, 3)));
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(2, 2), ploidyPDF.apply(new PloidyState(2, 2)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 5), ploidyPDF.apply(new PloidyState(0, 5)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(1, 4), ploidyPDF.apply(new PloidyState(1, 4)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(2, 3), ploidyPDF.apply(new PloidyState(2, 3)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 6), ploidyPDF.apply(new PloidyState(0, 6)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(1, 5), ploidyPDF.apply(new PloidyState(1, 5)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(2, 4), ploidyPDF.apply(new PloidyState(2, 4)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(3, 3), ploidyPDF.apply(new PloidyState(3, 3)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 7), ploidyPDF.apply(new PloidyState(0, 7)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(1, 6), ploidyPDF.apply(new PloidyState(1, 6)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(2, 5), ploidyPDF.apply(new PloidyState(2, 5)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(3, 4), ploidyPDF.apply(new PloidyState(3, 4)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 8), ploidyPDF.apply(new PloidyState(0, 8)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(1, 7), ploidyPDF.apply(new PloidyState(1, 7)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(2, 6), ploidyPDF.apply(new PloidyState(2, 6)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(3, 5), ploidyPDF.apply(new PloidyState(3, 5)));
//        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(4, 4), ploidyPDF.apply(new PloidyState(4, 4)));
        final PloidyStatePrior variantPloidyStatePrior = new PloidyStatePrior(unnormalizedLogProbabilityMassFunctionMap);

        final int numPopulations = 3;
        final int numCells = 200;

        final int numSamples = 100;
        final int numBurnIn = 50;

        final double concentrationPriorAlpha = 1.;
        final double concentrationPriorBeta = 10000.;
        final double variantSegmentFractionPriorAlpha = 4.;
        final double variantSegmentFractionPriorBeta = 4.;

        //run MCMC
        final TumorHeterogeneityModeller modeller = new TumorHeterogeneityModeller(
                segments, normalPloidyState, variantPloidyStatePrior,
                concentrationPriorAlpha, concentrationPriorBeta, variantSegmentFractionPriorAlpha, variantSegmentFractionPriorBeta,
                numPopulations, numCells, rng);
        modeller.fitMCMC(numSamples, numBurnIn);

        System.out.println();

        //check statistics of global-parameter posterior samples (i.e., posterior mode and standard deviation)
        final Map<TumorHeterogeneityParameter, PosteriorSummary> globalParameterPosteriorSummaries =
                modeller.getGlobalParameterPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);

        final PosteriorSummary concentrationPosteriorSummary = globalParameterPosteriorSummaries.get(TumorHeterogeneityParameter.CONCENTRATION);
        final double concentrationPosteriorCenter = concentrationPosteriorSummary.getCenter();
        final double concentrationPosteriorStandardDeviation = (concentrationPosteriorSummary.getUpper() - concentrationPosteriorSummary.getLower()) / 2;
        System.out.println("concentration: " + concentrationPosteriorCenter + " " + concentrationPosteriorStandardDeviation);
        System.out.println();

        final List<TumorHeterogeneityState.PopulationIndicators> populationIndicatorsSamples = modeller.getPopulationIndicatorsSamples();
        final List<TumorHeterogeneityState.PopulationFractions> populationFractionsSamples = modeller.getPopulationFractionsSamples();
        final List<TumorHeterogeneityState.VariantProfileCollection> variantProfileCollectionSamples = modeller.getVariantProfileCollectionSamples();

        for (int cellIndex = 0; cellIndex < numCells; cellIndex++) {
            final int ci = cellIndex;
            final MultiSet<Integer> populationCounts = new HashMultiSet<>(populationIndicatorsSamples.stream().map(s -> s.get(ci)).collect(Collectors.toList()));
            System.out.println("cell " + cellIndex + ": " + populationCounts);
        }
        System.out.println();

        for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
            final int pi = populationIndex;
            if (populationIndex != numPopulations - 1) {
                System.out.println("population " + populationIndex + " ==================================================");
                for (int segmentIndex = 0; segmentIndex < segments.size(); segmentIndex++) {
                    final int si = segmentIndex;
                    final double isVariantPosteriorProbability = (double) variantProfileCollectionSamples.stream()
                            .filter(vpc -> vpc.get(pi).isVariant(si))
                            .count() / (numSamples - numBurnIn);
                    System.out.println("segment " + segmentIndex + " (" + segments.get(segmentIndex).getInterval() + ") isVariant: " + isVariantPosteriorProbability);

                    for (int variantPloidyStateIndex = 0; variantPloidyStateIndex < variantPloidyStatePrior.numPloidyStates(); variantPloidyStateIndex++) {
                        final int vpsi = variantPloidyStateIndex;
                        final PloidyState variantPloidyState = variantPloidyStatePrior.ploidyStates().get(variantPloidyStateIndex);
                        final double variantPloidyStateProbability = (double) variantProfileCollectionSamples.stream()
                                .filter(vpc -> vpc.get(pi).variantPloidyStateIndex(si) == vpsi)
                                .count() / (numSamples - numBurnIn);
                        System.out.println(variantPloidyState + ": " + variantPloidyStateProbability);
                    }
                }
                System.out.println();
            }
        }
        System.out.println();

        for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
            final int pi = populationIndex;
            final double[] populationFractionSamples = populationFractionsSamples.stream().mapToDouble(s -> s.get(pi)).toArray();
            final double populationFractionPosteriorMean = new Mean().evaluate(populationFractionSamples);
            final double populationFractionPosteriorStandardDeviation = new StandardDeviation().evaluate(populationFractionSamples);

            if (populationIndex != numPopulations - 1) {
                final double variantSegmentFraction = (double) IntStream.range(0, segments.size()).filter(i ->
                        (double) variantProfileCollectionSamples.stream().filter(vpc -> vpc.get(pi).isVariant(i)).count() / (numSamples - numBurnIn) > 0.9)
                        .count() / segments.size();
                System.out.println("population fraction " + populationIndex + ": " + populationFractionPosteriorMean + " " + populationFractionPosteriorStandardDeviation + " " + variantSegmentFraction);
            } else {
                System.out.println("population fraction " + populationIndex + ": " + populationFractionPosteriorMean + " " + populationFractionPosteriorStandardDeviation);
            }
        }
    }
}