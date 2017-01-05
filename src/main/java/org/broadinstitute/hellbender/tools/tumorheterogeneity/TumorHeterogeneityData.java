package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.CauchyDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.DecileCollection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link DataCollection} for the tumor-heterogeneity model that allows the calculation of log posterior probabilities
 * for (copy ratio, minor-allele fraction) for each {@link ACNVModeledSegment}.  Fits a normal distribution to
 * the log_2 copy-ratio posterior deciles and a scaled beta distribution to the minor-allele fraction posterior deciles.
 * Allows for a copy-ratio noise constant (to model a constant contribution from stray reads).
 * Also stores any necessary prior information in {@link TumorHeterogeneityPriorCollection}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityData implements DataCollection {
    private static final Logger logger = LogManager.getLogger(TumorHeterogeneityData.class);

    //mathematical constants
    private static final double LN2 = GATKProtectedMathUtils.LN2;
    private static final double INV_LN2 = GATKProtectedMathUtils.INV_LN2;
    private static final double LN_LN2 = Math.log(LN2);

    //parameters for optimizer
    private static final double REL_TOLERANCE = 1E-5;
    private static final double ABS_TOLERANCE = 1E-10;
    private static final int NUM_MAX_EVAL = 1000;
    private static final double DEFAULT_SIMPLEX_STEP = 0.2;

    private static final double EPSILON = TumorHeterogeneityUtils.EPSILON;
    private static final double COPY_RATIO_EPSILON = TumorHeterogeneityUtils.COPY_RATIO_EPSILON; //below this, use mirrored minor-allele fraction posterior

    private static final MultivariateOptimizer optimizer = new SimplexOptimizer(REL_TOLERANCE, ABS_TOLERANCE);

    private final int numSegments;
    private final List<ACNVModeledSegment> segments;
    private final List<Integer> segmentLengths;
    private final List<Double> fractionalLengths;
    private final List<ACNVSegmentPosterior> segmentPosteriors;

    private final TumorHeterogeneityPriorCollection priors;

    public TumorHeterogeneityData(final List<ACNVModeledSegment> segments,
                                  final TumorHeterogeneityPriorCollection priors) {
        Utils.nonNull(segments);
        Utils.validateArg(segments.size() > 0, "Number of segments must be positive.");
        logger.info("Fitting copy-ratio and minor-allele-fraction posteriors to deciles...");
        numSegments = segments.size();
        this.segments = Collections.unmodifiableList(new ArrayList<>(segments));
        segmentLengths = Collections.unmodifiableList(segments.stream().map(s -> s.getInterval().size()).collect(Collectors.toList()));
        final double totalLength = segmentLengths.stream().mapToLong(Integer::longValue).sum();
        fractionalLengths = Collections.unmodifiableList(segmentLengths.stream().map(l -> (double) l / totalLength).collect(Collectors.toList()));
        segmentPosteriors = Collections.unmodifiableList(segments.stream().map(ACNVSegmentPosterior::new).collect(Collectors.toList()));
        this.priors = priors;
    }

    public TumorHeterogeneityData(final TumorHeterogeneityData data,
                                  final TumorHeterogeneityPriorCollection priors) {
        Utils.nonNull(data);
        numSegments = data.segments.size();
        this.segments = Collections.unmodifiableList(new ArrayList<>(data.segments));
        segmentLengths = Collections.unmodifiableList(new ArrayList<>(data.segmentLengths));
        fractionalLengths = Collections.unmodifiableList(new ArrayList<>(data.fractionalLengths));
        segmentPosteriors = Collections.unmodifiableList(new ArrayList<>(data.segmentPosteriors));
        this.priors = priors;
    }

    public int numSegments() {
        return numSegments;
    }

    public List<ACNVModeledSegment> segments() {
        return segments;
    }

    public double fractionalLength(final int segmentIndex) {
        Utils.validateArg(0 <= segmentIndex && segmentIndex < numSegments, "Segment index is not in valid range.");
        return fractionalLengths.get(segmentIndex);
    }

    public int length(final int segmentIndex) {
        Utils.validateArg(0 <= segmentIndex && segmentIndex < numSegments, "Segment index is not in valid range.");
        return segmentLengths.get(segmentIndex);
    }

    public double copyRatioLogDensity(final int segmentIndex, final double copyRatio, final double copyRatioNoiseConstant) {
        Utils.validateArg(0 <= segmentIndex && segmentIndex < numSegments, "Segment index is not in valid range.");
        Utils.validateArg(copyRatio >= 0., "Copy ratio must be non-negative.");
        Utils.validateArg(copyRatioNoiseConstant >= 0., "Copy-ratio noise constant must be non-negative.");
        return segmentPosteriors.get(segmentIndex).copyRatioLogDensity(copyRatio, copyRatioNoiseConstant);
    }

    public double logDensity(final int segmentIndex, final double copyRatio, final double minorAlleleFraction, final double copyRatioNoiseConstant) {
        Utils.validateArg(0 <= segmentIndex && segmentIndex < numSegments, "Segment index is not in valid range.");
        Utils.validateArg(copyRatio >= 0., "Copy ratio must be non-negative.");
        Utils.validateArg(0. <= minorAlleleFraction && minorAlleleFraction <= 0.5, "Minor-allele fraction must be in [0, 0.5].");
        Utils.validateArg(copyRatioNoiseConstant >= 0., "Copy-ratio noise constant must be non-negative.");
        if (copyRatio < COPY_RATIO_EPSILON) {
            //if copy ratio is below threshold, use mirrored minor-allele fraction posterior
            return segmentPosteriors.get(segmentIndex).copyRatioLogDensity(copyRatio, copyRatioNoiseConstant)
                    + FastMath.log(Math.max(EPSILON,
                    0.5 * (FastMath.exp(segmentPosteriors.get(segmentIndex).minorAlleleFractionLogDensity(minorAlleleFraction))
                            + FastMath.exp(segmentPosteriors.get(segmentIndex).minorAlleleFractionLogDensity(0.5 - minorAlleleFraction)))));
        }
        return segmentPosteriors.get(segmentIndex).copyRatioLogDensity(copyRatio, copyRatioNoiseConstant)
                + segmentPosteriors.get(segmentIndex).minorAlleleFractionLogDensity(minorAlleleFraction);
    }

    public TumorHeterogeneityPriorCollection priors() {
        return priors;
    }

    private final class ACNVSegmentPosterior {
        private final Function<Double, Double> log2CopyRatioPosteriorLogPDF;
        private final Function<Double, Double> minorAlleleFractionPosteriorLogPDF;

        ACNVSegmentPosterior(final ACNVModeledSegment segment) {
            logger.debug("Fitting segment: " + segment.getInterval());
            final List<Double> log2CopyRatioInnerDecilesList = segment.getSegmentMeanPosteriorSummary().getDeciles().getInner();
            final double[] log2CopyRatioInnerDeciles = Doubles.toArray(log2CopyRatioInnerDecilesList);
            logger.debug("Fitting normal distribution to inner deciles:\n" + log2CopyRatioInnerDecilesList);
            log2CopyRatioPosteriorLogPDF = fitCauchyLogPDFToInnerDeciles(log2CopyRatioInnerDeciles);

            final List<Double> minorAlleleFractionInnerDecilesList = segment.getMinorAlleleFractionPosteriorSummary().getDeciles().getInner();
            final double[] minorAlleleFractionInnerDeciles = Doubles.toArray(minorAlleleFractionInnerDecilesList);
            logger.debug("Fitting scaled beta distribution to inner deciles:\n" + minorAlleleFractionInnerDecilesList);
            final boolean isMinorAlleleFractionNaN = Double.isNaN(segment.getMinorAlleleFractionPosteriorSummary().getCenter());
            minorAlleleFractionPosteriorLogPDF = isMinorAlleleFractionNaN ?
                    point -> LN2 :       //flat over minor-allele fraction if NaN (i.e., no hets in segment) = log(1. / 0.5)
                    fitScaledBetaLogPDFToInnerDeciles(minorAlleleFractionInnerDeciles);
        }

        double copyRatioLogDensity(final double copyRatio, final double copyRatioNoiseConstant) {
            final double log2CopyRatio = FastMath.log(Math.max(EPSILON, copyRatio + copyRatioNoiseConstant)) * INV_LN2;
            return log2CopyRatioPosteriorLogPDF.apply(log2CopyRatio) - LN_LN2 - FastMath.log(Math.max(EPSILON, copyRatio));    //includes Jacobian: p(c) = p(log_2(c)) / (c * ln 2)
        }

        double minorAlleleFractionLogDensity(final double minorAlleleFraction) {
            final double minorAlleleFractionBounded = Math.max(EPSILON, Math.min(0.5 - EPSILON, minorAlleleFraction));
            return minorAlleleFractionPosteriorLogPDF.apply(minorAlleleFractionBounded);
        }

        //fit a Cauchy distribution to inner deciles (10th, 20th, ..., 90th percentiles) using least squares
        //and return the log PDF (for use as a posterior for log_2 copy ratio)
        private Function<Double, Double> fitCauchyLogPDFToInnerDeciles(final double[] innerDeciles) {
            final ObjectiveFunction innerDecilesL2LossFunction = new ObjectiveFunction(point -> {
                final double median = point[0];
                final double scale = FastMath.abs(point[1]);
                final CauchyDistribution cauchyDistribution = new CauchyDistribution(median, scale);
                final List<Double> cauchyInnerDeciles = IntStream.range(1, DecileCollection.NUM_DECILES - 1).boxed()
                        .map(i -> cauchyDistribution.inverseCumulativeProbability(i / 10.)).collect(Collectors.toList());
                return IntStream.range(0, DecileCollection.NUM_DECILES - 2)
                        .mapToDouble(i -> FastMath.pow(innerDeciles[i] - cauchyInnerDeciles.get(i), 2)).sum();
            });
            final double medianInitial = new Mean().evaluate(innerDeciles);
            final double scaleInitial = new StandardDeviation().evaluate(innerDeciles);
            logger.debug(String.format("Initial (median, scale) for Cauchy distribution: (%f, %f)", medianInitial, scaleInitial));
            final PointValuePair optimum = optimizer.optimize(
                            new MaxEval(NUM_MAX_EVAL),
                            innerDecilesL2LossFunction,
                            GoalType.MINIMIZE,
                            new InitialGuess(new double[]{medianInitial, scaleInitial}),
                            new NelderMeadSimplex(new double[]{DEFAULT_SIMPLEX_STEP, DEFAULT_SIMPLEX_STEP}));
            final double median = optimum.getPoint()[0];
            final double scale = FastMath.abs(optimum.getPoint()[1]);
            logger.debug(String.format("Final (median, scale) for Cauchy distribution: (%f, %f)", median, scale));
            return log2cr -> new CauchyDistribution(null, median, scale).logDensity(log2cr);
        }

        //fit a beta distribution to inner deciles (10th, 20th, ..., 90th percentiles) using least squares
        //and return the log PDF (scaled appropriately for use as a posterior for minor-allele fraction)
        private Function<Double, Double> fitScaledBetaLogPDFToInnerDeciles(final double[] innerDeciles) {
            final double[] scaledInnerDeciles = Arrays.stream(innerDeciles).map(d -> 2. * d).toArray();
            //scale minor-allele fraction deciles to [0, 1] and fit a beta distribution
            final ObjectiveFunction innerDecilesL2LossFunction = new ObjectiveFunction(point -> {
                final double alpha = FastMath.abs(point[0]);
                final double beta = FastMath.abs(point[1]);
                final BetaDistribution betaDistribution = new BetaDistribution(alpha, beta);
                final List<Double> betaInnerDeciles = IntStream.range(1, DecileCollection.NUM_DECILES - 1).boxed()
                        .map(i -> betaDistribution.inverseCumulativeProbability(i / 10.)).collect(Collectors.toList());
                return IntStream.range(0, DecileCollection.NUM_DECILES - 2)
                        .mapToDouble(i -> FastMath.pow(scaledInnerDeciles[i] - betaInnerDeciles.get(i), 2)).sum();
            });
            //use moment matching of deciles to initialize
            final double meanInitial = new Mean().evaluate(scaledInnerDeciles);
            final double varianceInitial = new Variance().evaluate(scaledInnerDeciles);
            final double commonFactor = FastMath.abs((meanInitial - meanInitial * meanInitial) / varianceInitial - 1.);
            final double alphaInitial = meanInitial * commonFactor;
            final double betaInitial = (1. - meanInitial) * commonFactor;
            logger.debug(String.format("Initial (alpha, beta) for scaled beta distribution: (%f, %f)", alphaInitial, betaInitial));
            final PointValuePair optimum = optimizer.optimize(
                    new MaxEval(NUM_MAX_EVAL),
                    innerDecilesL2LossFunction,
                    GoalType.MINIMIZE,
                    new InitialGuess(new double[]{alphaInitial, betaInitial}),
                    new NelderMeadSimplex(new double[]{DEFAULT_SIMPLEX_STEP, DEFAULT_SIMPLEX_STEP}));
            final double alpha = FastMath.abs(optimum.getPoint()[0]);
            final double beta = FastMath.abs(optimum.getPoint()[1]);
            logger.debug(String.format("Final (alpha, beta) for scaled beta distribution: (%f, %f)", alpha, beta));
            return maf -> new BetaDistribution(null, alpha, beta).logDensity(2. * maf) + LN2; //scale minor-allele fraction to [0, 1], including Jacobian factor
        }
    }
}