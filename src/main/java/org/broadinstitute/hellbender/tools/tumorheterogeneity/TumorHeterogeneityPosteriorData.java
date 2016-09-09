package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import akka.actor.FSM;
import breeze.stats.MeanAndVariance;
import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Doubles;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.analysis.MultivariateFunction;
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
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link DataCollection} for the tumor-heterogeneity model that allows the calculation of log posterior probabilities
 * for (copy ratio, minor-allele fraction) for each {@link ACNVModeledSegment}.  Given a list of
 * {@link ACNVModeledSegment}, a map from (copy ratio, minor-allele fraction) bins to
 * binned log posterior probability is constructed for each segment.  The binned log posterior probabilities
 * are determined by the copy-ratio and minor-allele-fraction posterior deciles for each segment and are uniform in
 * each bin.  These maps allow lookup via binary search of the log posterior probabilities
 * for arbitrary (copy ratio, minor allele fraction) points.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityPosteriorData implements DataCollection {
    private static final double INV_LN2 = GATKProtectedMathUtils.INV_LN2;
    private static final double REL_TOLERANCE = 1E-5;
    private static final double ABS_TOLERANCE = 1E-10;
    private static final int NUM_MAX_EVAL = 100;
    private static final double DEFAULT_SIMPLEX_STEP = 0.2;
    final MultivariateOptimizer optimizer = new SimplexOptimizer(REL_TOLERANCE, ABS_TOLERANCE);

    private final List<ACNVSegmentPosterior> segmentPosteriors;

    public TumorHeterogeneityPosteriorData(final List<ACNVModeledSegment> segments) {
        Utils.nonNull(segments);
        segmentPosteriors = segments.stream().map(ACNVSegmentPosterior::new).collect(Collectors.toList());
    }

    public double logProbability(final int segmentIndex, final double copyRatio, final double minorAlleleFraction) {
        Utils.validateArg(0 <= segmentIndex && segmentIndex < segmentPosteriors.size(), "Segment index is not in valid range.");
        Utils.validateArg(copyRatio >= 0, "Copy ratio must be non-negative.");
        Utils.validateArg(0 <= minorAlleleFraction && minorAlleleFraction <= 0.5, "Minor-allele fraction must be in [0, 0.5].");
        return segmentPosteriors.get(segmentIndex).logProbability(copyRatio, minorAlleleFraction);
    }
    
    private final class ACNVSegmentPosterior {
        private final NormalDistribution log2CopyRatioPosterior;
        private final BetaDistribution minorAlleleFractionPosterior;

        ACNVSegmentPosterior(final ACNVModeledSegment segment) {
            final double[] log2CopyRatioInnerDeciles = Doubles.toArray(segment.getSegmentMeanPosteriorSummary().getDeciles().getInner());
            final double[] minorAlleleFractionInnerDeciles = Doubles.toArray(segment.getMinorAlleleFractionPosteriorSummary().getDeciles().getInner());
            log2CopyRatioPosterior = fitNormalDistributionToInnerDeciles(log2CopyRatioInnerDeciles);
            minorAlleleFractionPosterior = fitBetaDistributionToInnerDeciles(minorAlleleFractionInnerDeciles);
        }

        double logProbability(final double copyRatio, final double minorAlleleFraction) {
            final double log2CopyRatio = Math.log(copyRatio) * INV_LN2;
            final double copyRatioPosteriorProbabilityDensity = log2CopyRatioPosterior.probability(log2CopyRatio) * INV_LN2 / copyRatio; //includes Jacobian: p(c) = p(log_2(c)) / (c * ln 2)
            final double minorAlleleFractionPosteriorProbabilityDensity = minorAlleleFractionPosterior.probability(minorAlleleFraction);
            return copyRatioPosteriorProbabilityDensity * minorAlleleFractionPosteriorProbabilityDensity;
        }

        private NormalDistribution fitNormalDistributionToInnerDeciles(final double[] innerDeciles) {
            final ObjectiveFunction innerDecilesL2LossFunction = new ObjectiveFunction(point -> {
                final double mean = point[0];
                final double standardDeviation = point[1];
                final NormalDistribution normalDistribution = new NormalDistribution(mean, standardDeviation);
                final List<Double> normalInnerDeciles = IntStream.range(1, DecileCollection.NUM_DECILES - 1).boxed()
                        .map(d -> normalDistribution.inverseCumulativeProbability(d / 10.)).collect(Collectors.toList());
                return IntStream.range(0, DecileCollection.NUM_DECILES - 2)
                        .mapToDouble(i -> FastMath.pow(innerDeciles[i] - normalInnerDeciles.get(i), 2)).sum();
            });
            final double meanInitial = new Mean().evaluate(innerDeciles);
            final double standarDeviationInitial = new StandardDeviation().evaluate(innerDeciles);
            final PointValuePair optimum = optimizer.optimize(
                            new MaxEval(NUM_MAX_EVAL),
                            innerDecilesL2LossFunction,
                            GoalType.MINIMIZE,
                            new InitialGuess(new double[]{meanInitial, standarDeviationInitial}),
                            new NelderMeadSimplex(new double[]{DEFAULT_SIMPLEX_STEP, DEFAULT_SIMPLEX_STEP}));
            final double mean = optimum.getPoint()[0];
            final double standardDeviation = optimum.getPoint()[1];
            return new NormalDistribution(mean, standardDeviation);
        }

        private BetaDistribution fitBetaDistributionToInnerDeciles(final double[] innerDeciles) {
            return new BetaDistribution(1., 1.);
        }
    }
}