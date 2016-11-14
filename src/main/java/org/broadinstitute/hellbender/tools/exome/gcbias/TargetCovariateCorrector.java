package org.broadinstitute.hellbender.tools.exome.gcbias;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Learn multiplicative correction factors as a function of some covariate content(such as GC or mappability) from coverage vs. covariate data.  Basically, learn a
 * regression curve of coverage vs. covariate, and divide by that curve to get covariate-corrected coverage.
 *
 * Our regression curve is obtained by filling covariate content bins of width 0.01 with the coverages of targets corresponding to each covariate
 * and then taking the median of coverages in each bin to get a robust estimate of the curve.  In order to smooth out bins with few data (i.e. extreme
 * covariate values that occur rarely) we then convolve these medians with an exponential kernel.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
class TargetCovariateCorrector {
    // Covariate bins are 0%, 1% . . . 100%
    private static final int NUMBER_OF_COVARIATE_BINS = 101;

    // scale (in units of CovariateContent from 0 to 1) over which covariate bias correlation decays
    // i.e. the bias at covariate content = 0.3 and at 0.2 are correlated ~exp(-0.1/correlationLength)
    // this is effectively a smoothness parameter used for regularizing the covariate bias estimate for
    // covariate content bins that have few targets)
    private static final double correlationLength = 0.02;
    private static final double correlationDecayRatePerBin = 1.0 / (correlationLength * NUMBER_OF_COVARIATE_BINS);

    // multiply by these to get a covariate correction as a function of the covariate
    private final double[] covariateCorrectionFactors;

    //Apache commons median doesn't work on empty arrays; this value is a placeholder to avoid exceptions
    private static final double DUMMY_VALUE_NEVER_USED = 1.0;

    /**
     * Learn multiplicative correction factors as a function of covariate from coverage vs. covariate data.  Basically, learn a
     * regression curve of coverage vs. covariate in order to divide by that curve later.
     *
     * @param covariateContents covariate content (from 0.0 to 1.0) of targets in {@code coverage}
     * @param coverage raw of proportional coverage
     */
    public TargetCovariateCorrector(final double[] covariateContents, final RealVector coverage) {
        Utils.nonNull(covariateContents);
        Utils.nonNull(coverage);
        Utils.validateArg(covariateContents.length > 0, "must have at least one datum");
        Utils.validateArg(covariateContents.length == coverage.getDimension(), "must have one gc value per coverage.");

        final List<List<Double>> coveragesByCovariate = new ArrayList<>(NUMBER_OF_COVARIATE_BINS);
        IntStream.range(0, NUMBER_OF_COVARIATE_BINS).forEach(n -> coveragesByCovariate.add(new ArrayList<>()));
        IntStream.range(0, covariateContents.length).forEach(n -> coveragesByCovariate.get(covariateContentToBinIndex(covariateContents[n])).add(coverage.getEntry(n)));
        covariateCorrectionFactors = calculateCorrectionFactors(coveragesByCovariate);
    }

    /**
     * As described above, calculate medians of each covariate bin and convolve with an exponential kernel.
     *
     * @param coveragesByCovariate list of coverages for each covariate bin
     * @return multiplicative correction factors for each covariate bin
     */
    private double[] calculateCorrectionFactors(final List<List<Double>> coveragesByCovariate) {
        final RealVector medians = new ArrayRealVector(coveragesByCovariate.stream().mapToDouble(TargetCovariateCorrector::medianOrDefault).toArray());
        return IntStream.range(0, NUMBER_OF_COVARIATE_BINS).mapToDouble(bin -> {
            final RealVector weights = new ArrayRealVector(IntStream.range(0, NUMBER_OF_COVARIATE_BINS)
                    .mapToDouble(n -> coveragesByCovariate.get(n).size() * Math.exp(-Math.abs(bin - n) * correlationDecayRatePerBin)).toArray());
            return weights.dotProduct(medians) / weights.getL1Norm();
        }).map(x -> 1/x).toArray();
    }

    /**
     *
     * @param inputCounts raw coverage before covariate correction
     * @param covariateContentByTarget      array of covariate contents, one per target of the input
     * @return covariate-corrected coverage
     */
    public static ReadCountCollection correctCoverage(final ReadCountCollection inputCounts, final double[] covariateContentByTarget) {
        // each column (sample) has its own covariate bias curve, hence its own covariateCorrector
        final List<TargetCovariateCorrector> covariateCorrectors = IntStream.range(0, inputCounts.columnNames().size())
                .mapToObj(n -> new TargetCovariateCorrector(covariateContentByTarget, inputCounts.counts().getColumnVector(n))).collect(Collectors.toList());

        // covariate-correct a copy of the input counts in-place
        final RealMatrix correctedCounts = inputCounts.counts().copy();
        correctedCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int target, int column, double coverage) {
                return covariateCorrectors.get(column).correctedCoverage(coverage, covariateContentByTarget[target]);
            }
        });

        // we would like the average correction factor to be 1.0 in the sense that average coverage before and after
        // correction should be equal
        final double[] columnNormalizationFactors = IntStream.range(0, inputCounts.columnNames().size())
                .mapToDouble(c -> inputCounts.counts().getColumnVector(c).getL1Norm() / correctedCounts.getColumnVector(c).getL1Norm()).toArray();
        correctedCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int target, int column, double coverage) {
                return coverage * columnNormalizationFactors[column];
            }
        });

        return new ReadCountCollection(inputCounts.targets(), inputCounts.columnNames(), correctedCounts);
    }

    private double correctedCoverage(final double coverage, final double covariateContent) {
        return covariateCorrectionFactors[covariateContentToBinIndex(covariateContent)] * coverage;
    }

    // return a median of coverages or dummy default value if no coverage exists at this covariate bin
    // this default is never used because empty bins get zero weight in {@code calculateCorrectionFactors}
    private static double medianOrDefault(final List<Double> list) {
        return list.size() > 0 ? new Median().evaluate(list.stream().mapToDouble(d->d).toArray()) : DUMMY_VALUE_NEVER_USED;
    }

    private static int covariateContentToBinIndex(final double covariateContent) {
        return (int) Math.round(covariateContent * (NUMBER_OF_COVARIATE_BINS -1));
    }
}
