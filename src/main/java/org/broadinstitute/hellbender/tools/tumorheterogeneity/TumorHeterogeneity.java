package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.exome.AllelicCNV;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
 * Model tumor heterogeneity as a Dirichlet mixture of subclones with copy-number variation,
 * starting from {@link AllelicCNV} output.
 * (Alpha version: only a single tumor population is assumed.)
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Model tumor heterogeneity as a Dirichlet mixture of subclones with copy-number variation, " +
                "starting from AllelicCNV output. (Alpha version: only a single tumor population is assumed.)",
        oneLineSummary = "Model tumor heterogeneity as a Dirichlet mixture of subclones with copy-number variation, " +
                "starting from AllelicCNV output",
        programGroup = CopyNumberProgramGroup.class
)
public final class TumorHeterogeneity extends SparkCommandLineProgram {
    private static final long serialVersionUID = 19738246L;

    private static final long RANDOM_SEED = 13;
    private static final int NUM_POPULATIONS_CLONAL = 2;
    private static final PloidyState NORMAL_PLOIDY_STATE = new PloidyState(1, 1);
    private static final double EPSILON = TumorHeterogeneityUtils.EPSILON;

    //fixes concentration to practically unity for clonal-only version
    private static final double CONCENTRATION_PRIOR_ALPHA_CLONAL = 1E6;
    private static final double CONCENTRATION_PRIOR_BETA_CLONAL = 1E6;

    private static final double COPY_RATIO_NOISE_CONSTANT_MIN = TumorHeterogeneityUtils.COPY_RATIO_NOISE_CONSTANT_MIN;
    private static final double COPY_RATIO_NOISE_CONSTANT_MAX = TumorHeterogeneityUtils.COPY_RATIO_NOISE_CONSTANT_MAX;

    private static final double COPY_RATIO_NOISE_FACTOR_MIN = TumorHeterogeneityUtils.COPY_RATIO_NOISE_FACTOR_MIN;
    private static final double COPY_RATIO_NOISE_FACTOR_MAX = TumorHeterogeneityUtils.COPY_RATIO_NOISE_FACTOR_MAX;

    private static final double MINOR_ALLELE_FRACTION_NOISE_FACTOR_MIN = TumorHeterogeneityUtils.MINOR_ALLELE_FRACTION_NOISE_FACTOR_MIN;
    private static final double MINOR_ALLELE_FRACTION_NOISE_FACTOR_MAX = TumorHeterogeneityUtils.MINOR_ALLELE_FRACTION_NOISE_FACTOR_MAX;

    //filename tags for output
    protected static final String CLONAL_SAMPLES_FILE_SUFFIX = ".th.clonal.samples.tsv";
    protected static final String CLONAL_PROFILES_FILE_SUFFIX = ".th.clonal.profiles.tsv";
    protected static final String CLONAL_SUMMARY_FILE_SUFFIX = ".th.clonal.summary.tsv";

    //CLI arguments
    protected static final String OUTPUT_PREFIX_LONG_NAME = "outputPrefix";
    protected static final String OUTPUT_PREFIX_SHORT_NAME = "pre";

    protected static final String MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME = "maxAllelicCopyNumberClonal";
    protected static final String MAX_ALLELIC_COPY_NUMBER_CLONAL_SHORT_NAME = "maxACNClonal";

    protected static final String NUM_WALKERS_CLONAL_LONG_NAME = "numWalkersClonal";
    protected static final String NUM_WALKERS_CLONAL_SHORT_NAME = "numWalkClonal";

    protected static final String NUM_SAMPLES_CLONAL_LONG_NAME = "numSamplesClonal";
    protected static final String NUM_SAMPLES_CLONAL_SHORT_NAME = "numSampClonal";

    protected static final String NUM_BURN_IN_CLONAL_LONG_NAME = "numBurnInClonal";
    protected static final String NUM_BURN_IN_CLONAL_SHORT_NAME = "numBurnClonal";

    protected static final String COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_LONG_NAME = "copyRatioNoiseConstantPriorAlpha";
    protected static final String COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_SHORT_NAME = "crConstAlpha";

    protected static final String COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_LONG_NAME = "copyRatioNoiseConstantPriorBeta";
    protected static final String COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_SHORT_NAME = "crConstBeta";

    protected static final String COPY_RATIO_NOISE_FACTOR_PRIOR_ALPHA_LONG_NAME = "copyRatioNoiseFactorPriorAlpha";
    protected static final String COPY_RATIO_NOISE_FACTOR_PRIOR_ALPHA_SHORT_NAME = "crFactorAlpha";

    protected static final String COPY_RATIO_NOISE_FACTOR_PRIOR_BETA_LONG_NAME = "copyRatioNoiseFactorPriorBeta";
    protected static final String COPY_RATIO_NOISE_FACTOR_PRIOR_BETA_SHORT_NAME = "crFactorBeta";

    protected static final String MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_ALPHA_LONG_NAME = "minorAlleleFractionNoiseFactorPriorAlpha";
    protected static final String MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_ALPHA_SHORT_NAME = "mafFactorAlpha";

    protected static final String MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_BETA_LONG_NAME = "minorAlleleFractionNoiseFactorPriorBeta";
    protected static final String MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_BETA_SHORT_NAME = "mafFactorBeta";

    protected static final String PLOIDY_MISMATCH_PENALTY_LONG_NAME = "ploidyMismatchPenalty";
    protected static final String PLOIDY_MISMATCH_PENALTY_SHORT_NAME = "mismatchPen";

    protected static final String PLOIDY_STATE_PRIOR_HOMOZYGOUS_DELETION_PENALTY_LONG_NAME = "homozygousDeletionPenalty";
    protected static final String PLOIDY_STATE_PRIOR_HOMOZYGOUS_DELETION_PENALTY_SHORT_NAME = "homDelPen";

    protected static final String PLOIDY_STATE_PRIOR_CHANGE_PENALTY_LONG_NAME = "changePenalty";
    protected static final String PLOIDY_STATE_PRIOR_CHANGE_PENALTY_SHORT_NAME = "changePen";

    protected static final String MODE_PURITY_BIN_SIZE_LONG_NAME = "purityBinSize";
    protected static final String MODE_PURITY_BIN_SIZE_SHORT_NAME = "purityBin";

    protected static final String MODE_PLOIDY_BIN_SIZE_LONG_NAME = "ploidyBinSize";
    protected static final String MODE_PLOIDY_BIN_SIZE_SHORT_NAME = "ploidyBin";

    @Argument(
            doc = "Input file for AllelicCNV result.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME
    )
    protected File allelicCNVFile;

    @Argument(
            doc = "Prefix for output files.",
            fullName = OUTPUT_PREFIX_LONG_NAME,
            shortName = OUTPUT_PREFIX_SHORT_NAME
    )
    protected String outputPrefix;

    @Argument(
            doc = "Maximum allelic copy number for clonal model.",
            fullName = MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME,
            shortName = MAX_ALLELIC_COPY_NUMBER_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int maxAllelicCopyNumberClonal = 5;

    @Argument(
            doc = "Number of walkers for MCMC ensemble.",
            fullName = NUM_WALKERS_CLONAL_LONG_NAME,
            shortName = NUM_WALKERS_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int numWalkers = 50;

    @Argument(
            doc = "Total number of MCMC ensemble samples for clonal model. " +
                    "(Total number of samples will be number of walkers in ensemble  multiplied by this number.)",
            fullName = NUM_SAMPLES_CLONAL_LONG_NAME,
            shortName = NUM_SAMPLES_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int numSamplesClonal = 100;

    @Argument(
            doc = "Number of burn-in ensemble samples to discard for clonal model. " +
                    "(Total number of samples will be number of walkers in ensemble multiplied by this number.)",
            fullName = NUM_BURN_IN_CLONAL_LONG_NAME,
            shortName = NUM_BURN_IN_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int numBurnInClonal = 50;

    @Argument(
            doc = "Alpha hyperparameter for Gamma-distribution prior on copy-ratio noise-constant parameter.",
            fullName = COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_LONG_NAME,
            shortName = COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_SHORT_NAME,
            optional = true
    )
    protected double copyRatioNoiseConstantPriorAlpha = 1;

    @Argument(
            doc = "Beta hyperparameter for Gamma-distribution prior on copy-ratio noise-constant parameter.",
            fullName = COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_LONG_NAME,
            shortName = COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_SHORT_NAME,
            optional = true
    )
    protected double copyRatioNoiseConstantPriorBeta = 1E2;

    @Argument(
            doc = "Alpha hyperparameter for Beta-distribution prior on copy-ratio noise-factor parameter.",
            fullName = COPY_RATIO_NOISE_FACTOR_PRIOR_ALPHA_LONG_NAME,
            shortName = COPY_RATIO_NOISE_FACTOR_PRIOR_ALPHA_SHORT_NAME,
            optional = true
    )
    protected double copyRatioNoiseFactorPriorAlpha = 1E3;

    @Argument(
            doc = "Beta hyperparameter for Beta-distribution prior on copy-ratio noise-factor parameter.",
            fullName = COPY_RATIO_NOISE_FACTOR_PRIOR_BETA_LONG_NAME,
            shortName = COPY_RATIO_NOISE_FACTOR_PRIOR_BETA_SHORT_NAME,
            optional = true
    )
    protected double copyRatioNoiseFactorPriorBeta = 1;

    @Argument(
            doc = "Alpha hyperparameter for Beta-distribution prior on minor-allele-fraction noise-factor parameter.",
            fullName = MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_ALPHA_LONG_NAME,
            shortName = MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_ALPHA_SHORT_NAME,
            optional = true
    )
    protected double minorAlleleFractionNoiseFactorPriorAlpha = 1E3;

    @Argument(
            doc = "Beta hyperparameter for Beta-distribution prior on minor-allele-fraction noise-factor parameter.",
            fullName = MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_BETA_LONG_NAME,
            shortName = MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_BETA_SHORT_NAME,
            optional = true
    )
    protected double minorAlleleFractionNoiseFactorPriorBeta = 1;

    @Argument(
            doc = "Penalty factor for ploidy mismatch in proposal of variant profiles. " +
                    "(A strong (i.e., large) penalty factor will reduce the number of solutions found.)",
            fullName = PLOIDY_MISMATCH_PENALTY_LONG_NAME,
            shortName = PLOIDY_MISMATCH_PENALTY_SHORT_NAME,
            optional = true
    )
    protected double ploidyMismatchPenalty = 1E4;

    @Argument(
        doc = "Penalty factor for homozygous deletion per base in ploidy-state prior. " +
                "(A strong (i.e., large) penalty factor will decrease sensitivity at low purity.)",
        fullName = PLOIDY_STATE_PRIOR_HOMOZYGOUS_DELETION_PENALTY_LONG_NAME,
        shortName = PLOIDY_STATE_PRIOR_HOMOZYGOUS_DELETION_PENALTY_SHORT_NAME,
        optional = true
    )
    protected double ploidyStatePriorCompleteDeletionPenalty = 1E-6;

    @Argument(
            doc = "Penalty factor for copy change per base  in ploidy-state prior. " +
                    "(A strong (i.e., large) penalty factor will decrease sensitivity at low purity.)",
            fullName = PLOIDY_STATE_PRIOR_CHANGE_PENALTY_LONG_NAME,
            shortName = PLOIDY_STATE_PRIOR_CHANGE_PENALTY_SHORT_NAME,
            optional = true
    )
    protected double ploidyStatePriorChangePenalty = 1E-6;

    @Argument(
            doc = "Purity bin size for identifying samples at posterior mode.",
            fullName = MODE_PURITY_BIN_SIZE_LONG_NAME,
            shortName = MODE_PURITY_BIN_SIZE_SHORT_NAME,
            optional = true
    )
    protected double purityModeBinSize = 0.025;

    @Argument(
            doc = "Ploidy bin size for identifying samples at posterior mode.",
            fullName = MODE_PLOIDY_BIN_SIZE_LONG_NAME,
            shortName = MODE_PLOIDY_BIN_SIZE_SHORT_NAME,
            optional = true
    )
    protected double ploidyModeBinSize = 0.1;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        validateArguments();

        //initialize output files
        final File samplesFileClonal = new File(outputPrefix + CLONAL_SAMPLES_FILE_SUFFIX);
        final File profilesFileClonal = new File(outputPrefix + CLONAL_PROFILES_FILE_SUFFIX);
        final File summaryFileClonal = new File(outputPrefix + CLONAL_SUMMARY_FILE_SUFFIX);

        //load ACNV segments from input file
        final List<ACNVModeledSegment> segments = SegmentUtils.readACNVModeledSegmentFile(allelicCNVFile);

        //construct priors using input parameters
        final PloidyStatePrior ploidyStatePriorClonal = calculatePloidyStatePrior(ploidyStatePriorCompleteDeletionPenalty, ploidyStatePriorChangePenalty, maxAllelicCopyNumberClonal);
        final TumorHeterogeneityPriorCollection priorsClonal = new TumorHeterogeneityPriorCollection(
                NORMAL_PLOIDY_STATE, ploidyStatePriorClonal,
                CONCENTRATION_PRIOR_ALPHA_CLONAL, CONCENTRATION_PRIOR_BETA_CLONAL,
                copyRatioNoiseConstantPriorAlpha, copyRatioNoiseConstantPriorBeta,
                copyRatioNoiseFactorPriorAlpha, copyRatioNoiseFactorPriorBeta,
                minorAlleleFractionNoiseFactorPriorAlpha, minorAlleleFractionNoiseFactorPriorBeta,
                ploidyMismatchPenalty);
        //initialize data collection from ACNV input and priors
        final TumorHeterogeneityData data = new TumorHeterogeneityData(segments, priorsClonal);

        //initialize modeller and run MCMC
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));
        final TumorHeterogeneityModeller clonalModeller = new TumorHeterogeneityModeller(data, NUM_POPULATIONS_CLONAL, numWalkers, rng);
        clonalModeller.fitMCMC(numSamplesClonal, numBurnInClonal);

        //initialize writer
        final TumorHeterogeneityModellerWriter writer = new TumorHeterogeneityModellerWriter(clonalModeller);

        //write all MCMC samples to file
        writer.writePopulationFractionAndPloidySamples(samplesFileClonal);
        logger.info("Population-fraction--ploidy MCMC samples output to " + samplesFileClonal + ".");

        //identify samples in purity-ploidy bin centered on posterior mode
        final List<Integer> indicesOfSamplesAtMode = clonalModeller.identifySamplesAtMode(purityModeBinSize, ploidyModeBinSize);

        //average variant profiles of identified samples and write to file
        logger.info("Calculating averaged variant profiles at posterior mode...");
        writer.writeAveragedProfiles(profilesFileClonal, indicesOfSamplesAtMode);
        logger.info("Averaged variant profiles at posterior mode output to " + profilesFileClonal + ".");

        //calculate posterior summaries for global parameters from identified samples and write to file
        logger.info("Calculating summary for posterior mode...");
        writer.writePosteriorSummaries(summaryFileClonal, indicesOfSamplesAtMode, ctx);
        logger.info("Calculating summary for posterior-mode variant profiles output to " + summaryFileClonal + ".");

        logger.info("SUCCESS: TumorHeterogeneity run complete.");
    }

    //we implement a simple prior on ploidy states that penalizes copy-number changes, with an additional penalty for homozygous deletions
    private static PloidyStatePrior calculatePloidyStatePrior(final double ploidyStatePriorCompleteDeletionPenalty,
                                                              final double ploidyStatePriorChangePenalty,
                                                              final int maxAllelicCopyNumber) {
        final Function<PloidyState, Double> ploidyLogPDF = ps ->
                        Math.log(Math.max(EPSILON, 1. - ploidyStatePriorCompleteDeletionPenalty)) * (ps.m() == 0 && ps.n() == 0 ? 1 : 0)
                                + Math.log(Math.max(EPSILON, 1. - ploidyStatePriorChangePenalty)) * (Math.abs(NORMAL_PLOIDY_STATE.m() - ps.m()) + Math.abs(NORMAL_PLOIDY_STATE.n() - ps.n()));
        final Map<PloidyState, Double> unnormalizedLogProbabilityMassFunctionMap = new LinkedHashMap<>();
        for (int n = 0; n <= maxAllelicCopyNumber; n++) {
            for (int m = 0; m <= n; m++) {
                unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(m, n), ploidyLogPDF.apply(new PloidyState(m, n)));
            }
        }
        return new PloidyStatePrior(unnormalizedLogProbabilityMassFunctionMap);
    }

    //validate CLI arguments
    private void validateArguments() {
        Utils.regularReadableUserFile(allelicCNVFile);
        Utils.validateArg(maxAllelicCopyNumberClonal > 0, MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME + " must be positive.");
        Utils.validateArg(numSamplesClonal > 0, NUM_SAMPLES_CLONAL_LONG_NAME + " must be positive.");
        Utils.validateArg(numBurnInClonal >= 0 && numBurnInClonal < numSamplesClonal, NUM_BURN_IN_CLONAL_LONG_NAME + " must be non-negative and strictly less than " + NUM_SAMPLES_CLONAL_LONG_NAME);
        validatePriorHyperparameters(
                copyRatioNoiseConstantPriorAlpha, COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_LONG_NAME,
                copyRatioNoiseConstantPriorBeta, COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_LONG_NAME,
                COPY_RATIO_NOISE_CONSTANT_MIN, COPY_RATIO_NOISE_CONSTANT_MAX,
                (alpha, beta) -> alpha / beta);
        validatePriorHyperparameters(
                copyRatioNoiseFactorPriorAlpha, COPY_RATIO_NOISE_FACTOR_PRIOR_ALPHA_LONG_NAME,
                copyRatioNoiseFactorPriorBeta, COPY_RATIO_NOISE_FACTOR_PRIOR_BETA_LONG_NAME,
                COPY_RATIO_NOISE_FACTOR_MIN, COPY_RATIO_NOISE_FACTOR_MAX,
                (alpha, beta) -> alpha / (alpha + beta));
        validatePriorHyperparameters(
                minorAlleleFractionNoiseFactorPriorAlpha, MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_ALPHA_LONG_NAME,
                minorAlleleFractionNoiseFactorPriorBeta, MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_BETA_LONG_NAME,
                MINOR_ALLELE_FRACTION_NOISE_FACTOR_MIN, MINOR_ALLELE_FRACTION_NOISE_FACTOR_MAX,
                (alpha, beta) -> alpha / (alpha + beta));
        Utils.validateArg(0. <= ploidyMismatchPenalty, PLOIDY_MISMATCH_PENALTY_LONG_NAME + " must be non-negative.");
        Utils.validateArg(0. <= ploidyStatePriorCompleteDeletionPenalty && ploidyStatePriorCompleteDeletionPenalty <= 1.,
                PLOIDY_STATE_PRIOR_HOMOZYGOUS_DELETION_PENALTY_LONG_NAME + " must be in [0, 1].");
        Utils.validateArg(0. <= ploidyStatePriorChangePenalty && ploidyStatePriorChangePenalty <= 1.,
                PLOIDY_STATE_PRIOR_CHANGE_PENALTY_LONG_NAME + " must be in [0, 1].");
        Utils.validateArg(0. <= purityModeBinSize && purityModeBinSize <= 1.,
                "Invalid purity bin size for determining mode.");
        Utils.validateArg(0. <= ploidyModeBinSize && ploidyModeBinSize <= 2 * maxAllelicCopyNumberClonal,
                "Invalid ploidy bin size for determining mode.");
    }

    //validate hyperparameters for Beta or Gamma distributions
    private static void validatePriorHyperparameters(final double alpha, final String alphaName,
                                                     final double beta, final String betaName,
                                                     final double min, final double max,
                                                     final BiFunction<Double, Double, Double> calculateMeanFromHyperparameters) {
        Utils.validateArg(alpha > 0, alphaName + " must be positive.");
        Utils.validateArg(beta > 0, betaName + " must be positive.");
        final double boundedQuantity = calculateMeanFromHyperparameters.apply(alpha, beta);
        Utils.validateArg(min < boundedQuantity && boundedQuantity < max, "Mean calculated from hyperparameters alpha and beta must be in (" + min + ", " + max + ").");
    }
}