package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Doubles;
import org.apache.commons.math3.distribution.BetaDistribution;
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
import org.apache.commons.math3.stat.descriptive.moment.SecondMoment;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;

import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link DataCollection} for the tumor-heterogeneity model that allows the calculation of log posterior probabilities
 * for (copy ratio, minor-allele fraction) for each {@link ACNVModeledSegment}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityPosteriorData implements DataCollection {
    private static final double LN2 = GATKProtectedMathUtils.LN2;
    private static final double INV_LN2 = GATKProtectedMathUtils.INV_LN2;
    private static final double LN_LN2 = Math.log(LN2);
    private static final double REL_TOLERANCE = 1E-5;
    private static final double ABS_TOLERANCE = 1E-10;
    private static final int NUM_MAX_EVAL = 100;
    private static final double DEFAULT_SIMPLEX_STEP = 0.2;

    public static final Logger logger = LogManager.getLogger(TumorHeterogeneityPosteriorData.class);
    private static final MultivariateOptimizer optimizer = new SimplexOptimizer(REL_TOLERANCE, ABS_TOLERANCE);

    private final List<ACNVSegmentPosterior> segmentPosteriors;

    public TumorHeterogeneityPosteriorData(final List<ACNVModeledSegment> segments) {
        Utils.nonNull(segments);
        segmentPosteriors = segments.stream().map(ACNVSegmentPosterior::new).collect(Collectors.toList());
    }

    public double logDensity(final int segmentIndex, final double copyRatio, final double minorAlleleFraction) {
        Utils.validateArg(0 <= segmentIndex && segmentIndex < segmentPosteriors.size(), "Segment index is not in valid range.");
        Utils.validateArg(copyRatio >= 0, "Copy ratio must be non-negative.");
        Utils.validateArg(0 <= minorAlleleFraction && minorAlleleFraction <= 0.5, "Minor-allele fraction must be in [0, 0.5].");
        return segmentPosteriors.get(segmentIndex).logDensity(copyRatio, minorAlleleFraction);
    }
    
    private final class ACNVSegmentPosterior {
        private final boolean isMinorAlleleFractionNaN;
        private final Function<Double, Double> log2CopyRatioPosteriorLogPDF;
        private final Function<Double, Double> minorAlleleFractionPosteriorLogPDF;

        ACNVSegmentPosterior(final ACNVModeledSegment segment) {
            final double[] log2CopyRatioInnerDeciles = Doubles.toArray(segment.getSegmentMeanPosteriorSummary().getDeciles().getInner());
            final double[] minorAlleleFractionInnerDeciles = Doubles.toArray(segment.getMinorAlleleFractionPosteriorSummary().getDeciles().getInner());
            isMinorAlleleFractionNaN = Double.isNaN(segment.getMinorAlleleFractionPosteriorSummary().getCenter());
            log2CopyRatioPosteriorLogPDF = fitNormalLogPDFToInnerDeciles(log2CopyRatioInnerDeciles);
            minorAlleleFractionPosteriorLogPDF = isMinorAlleleFractionNaN ?
                    f -> LN2 :       //flat over minor-allele fraction if NaN (i.e., no hets in segment) = log(1. / 0.5)
                    fitBetaLogPDFToInnerDeciles(minorAlleleFractionInnerDeciles);
        }

        double logDensity(final double copyRatio, final double minorAlleleFraction) {
            final double log2CopyRatio = Math.log(copyRatio) * INV_LN2;
            final double copyRatioPosteriorLogDensity =
                    log2CopyRatioPosteriorLogPDF.apply(log2CopyRatio) - LN_LN2 - Math.log(copyRatio);    //includes Jacobian: p(c) = p(log_2(c)) / (c * ln 2)
            final double minorAlleleFractionPosteriorLogDensity = minorAlleleFractionPosteriorLogPDF.apply(minorAlleleFraction);
            return copyRatioPosteriorLogDensity + minorAlleleFractionPosteriorLogDensity;
        }

        private Function<Double, Double> fitNormalLogPDFToInnerDeciles(final double[] innerDeciles) {
            final ObjectiveFunction innerDecilesL2LossFunction = new ObjectiveFunction(point -> {
                final double mean = point[0];
                final double standardDeviation = Math.abs(point[1]);
                final NormalDistribution normalDistribution = new NormalDistribution(mean, standardDeviation);
                final List<Double> normalInnerDeciles = IntStream.range(1, DecileCollection.NUM_DECILES - 1).boxed()
                        .map(i -> normalDistribution.inverseCumulativeProbability(i / 10.)).collect(Collectors.toList());
                return IntStream.range(0, DecileCollection.NUM_DECILES - 2)
                        .mapToDouble(i -> FastMath.pow(innerDeciles[i] - normalInnerDeciles.get(i), 2)).sum();
            });
            final double meanInitial = new Mean().evaluate(innerDeciles);
            final double standardDeviationInitial = new StandardDeviation().evaluate(innerDeciles);
            final PointValuePair optimum = optimizer.optimize(
                            new MaxEval(NUM_MAX_EVAL),
                            innerDecilesL2LossFunction,
                            GoalType.MINIMIZE,
                            new InitialGuess(new double[]{meanInitial, standardDeviationInitial}),
                            new NelderMeadSimplex(new double[]{DEFAULT_SIMPLEX_STEP, DEFAULT_SIMPLEX_STEP}));
            final double mean = optimum.getPoint()[0];
            final double standardDeviation = Math.abs(optimum.getPoint()[1]);
            return log2cr -> new NormalDistribution(mean, standardDeviation).logDensity(log2cr);
        }

        private Function<Double, Double> fitBetaLogPDFToInnerDeciles(final double[] innerDeciles) {
            //scale minor-allele fraction deciles to [0, 1] and fit a beta distribution
            final ObjectiveFunction innerDecilesL2LossFunction = new ObjectiveFunction(point -> {
                final double alpha = Math.abs(point[0]);
                final double beta = Math.abs(point[1]);
                final BetaDistribution betaDistribution = new BetaDistribution(alpha, beta);
                final List<Double> betaInnerDeciles = IntStream.range(1, DecileCollection.NUM_DECILES - 1).boxed()
                        .map(i -> betaDistribution.inverseCumulativeProbability(i / 10.)).collect(Collectors.toList());
                return IntStream.range(0, DecileCollection.NUM_DECILES - 2)
                        .mapToDouble(i -> FastMath.pow(2 * innerDeciles[i] - betaInnerDeciles.get(i), 2)).sum();
            });
            //use moment matching of deciles to initialize
            final double meanInitial = new Mean().evaluate(innerDeciles);
            final double secondMomentInitial = new SecondMoment().evaluate(innerDeciles);
            final double denominator = secondMomentInitial - meanInitial * meanInitial;
            final double alphaInitial = meanInitial * (meanInitial - secondMomentInitial) / denominator;
            final double betaInitial = (1. - meanInitial) * (meanInitial - secondMomentInitial) / denominator;
            final PointValuePair optimum = optimizer.optimize(
                    new MaxEval(NUM_MAX_EVAL),
                    innerDecilesL2LossFunction,
                    GoalType.MINIMIZE,
                    new InitialGuess(new double[]{alphaInitial, betaInitial}),
                    new NelderMeadSimplex(new double[]{DEFAULT_SIMPLEX_STEP, DEFAULT_SIMPLEX_STEP}));
            final double alpha = Math.abs(optimum.getPoint()[0]);
            final double beta = Math.abs(optimum.getPoint()[1]);
            return maf -> new BetaDistribution(alpha, beta).logDensity(2. * maf) + LN2; //scale minor-allele fraction to [0, 1], including Jacobian factor
        }
    }
}