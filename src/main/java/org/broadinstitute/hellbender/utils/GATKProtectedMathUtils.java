package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Extra MathUtils that should be moved to gatk-public
 * Created by davidben on 1/22/16.
 */
public class GATKProtectedMathUtils implements Serializable {

    private static final long serialVersionUID = -7587940242569731513L;

    private GATKProtectedMathUtils() {
    }

    public static final double INV_LN2 = 1.0 / Math.log(2.0);

    /**
     * Computes $\log(\sum_i e^{a_i})$ trying to avoid underflow issues by using the log-sum-exp trick.
     *
     * <p>
     * This trick consists on shift all the log values by the maximum so that exponent values are
     * much larger (close to 1) before they are summed. Then the result is shifted back down by
     * the same amount in order to obtain the correct value.
     * </p>
     * @return any double value.
     */
    public static double logSumExp(final double ... values) {
        double max = MathUtils.arrayMax(Utils.nonNull(values));
        double sum = 0.0;
        for (int i = 0; i < values.length; ++i) {
            if (values[i] != Double.NEGATIVE_INFINITY) {
                sum += Math.exp(values[i] - max);
            }
        }
        return max + Math.log(sum);
    }

    public static double logSumExp(final Collection<Double> values) {
        double max = Collections.max(values);
        double sum = 0.0;
        for (final double val : values) {
            if (val != Double.NEGATIVE_INFINITY) {
                sum += Math.exp(val - max);
            }
        }
        return max + Math.log(sum);
    }

    public static double mean(final double ... values) {
        Utils.nonNull(values);
        return MathUtils.mean(values, 0, values.length);
    }

    public static double[] rowMeans(final RealMatrix matrix) {
        Utils.nonNull(matrix);
        return IntStream.range(0, matrix.getRowDimension())
                .mapToDouble(r -> mean(matrix.getRow(r))).toArray();
    }

    public static double[] rowVariances(final RealMatrix matrix) {
        Utils.nonNull(matrix);
        final Variance varianceEvaluator = new Variance();
        return IntStream.range(0, matrix.getRowDimension())
                .mapToDouble(r -> varianceEvaluator.evaluate(matrix.getRow(r))).toArray();
    }

    /**
     * Calculates the standard deviation per row from a matrix.
     * @param matrix the input matrix.
     * @return never {@code null}, an array with as many positions as rows in {@code matrix}.
     * @throws IllegalArgumentException if {@code matrix} is {@code null}.
     */
    public static double[] rowStdDevs(final RealMatrix matrix) {
        Utils.nonNull(matrix);
        final Variance varianceEvaluator = new Variance();
        return IntStream.range(0, matrix.getRowDimension())
                .mapToDouble(r -> Math.sqrt(varianceEvaluator.evaluate(matrix.getRow(r)))).toArray();
    }

    /**
     * Calculates the mean per column from a matrix.
     * @param matrix the input matrix.
     * @return never {@code null}, an array with as many positions as columns in {@code matrix}.
     * @throws IllegalArgumentException if {@code matrix} is {@code null}.
     */
    public static double[] columnMeans(final RealMatrix matrix) {
        Utils.nonNull(matrix);
        return IntStream.range(0, matrix.getColumnDimension())
                .mapToDouble(c -> mean(matrix.getColumn(c))).toArray();
    }

    /**
     * Calculates the standard deviation per column from a matrix.
     * @param matrix the input matrix.
     * @return never {@code null}, an array with as many positions as columns in {@code matrix}.
     * @throws IllegalArgumentException if {@code matrix} is {@code null}.
     */
    public static double[] columnStdDevs(final RealMatrix matrix) {
        Utils.nonNull(matrix);
        final Variance varianceEvaluator = new Variance();
        return IntStream.range(0, matrix.getColumnDimension())
                .mapToDouble(c -> Math.sqrt(varianceEvaluator.evaluate(matrix.getColumn(c)))).toArray();
    }

    /**
     * Calculate the standard deviation of a collection of {@link Number} instances.
     * @param values the input values.
     * @return the standard deviation.
     * @throws IllegalArgumentException if {@code values} is {@code null} or it contains {@code null}.
     */
    public static double stdDev(final Collection<? extends Number> values) {
        Utils.nonNull(values);
        Utils.containsNoNull(values, "input values must not contain a null");
        final double[] doubleValues = values.stream()
                .mapToDouble(Number::doubleValue).toArray();
        return stdDev(doubleValues);
    }

    /**
     * Calculate the standard deviation of a double array.
     * @param values the input values.
     * @return the standard deviation.
     * @throws IllegalArgumentException if {@code values} is {@code null}.
     */
    public static double stdDev(final double ... values) {
        Utils.nonNull(values);
        return Math.sqrt(new Variance().evaluate(values));
    }

    /**
     * Create a set number of points, linearly spaced, between a minimum and maximum.
     *
     * Inspired by http://stackoverflow.com/questions/6208878/java-version-of-matlabs-linspace
     *
     * NaN are allowed, but will likely give useless results.
     *
     * @param min starting value
     * @param max ending value
     * @param points number of points, must be greater than -1
     * @return Never {@code null}
     */
    public static double[] createEvenlySpacedPoints(final double min, final double max, int points) {
        ParamUtils.isPositiveOrZero(points, "Number of points must be >= 0");
        if (points == 1) {
            return new double[] {max};
        }
        if (points == 0) {
            return new double[] {};
        }
        return IntStream.range(0, points).mapToDouble(n -> min + n * (max - min) / (points - 1)).toArray();
    }

    // given a list of options and a function for calculating their probabilities (these must sum to 1 over the whole list)
    // randomly choose one option from the implied categorical distribution
    public static <E> E randomSelect(final List<E> choices, final Function<E, Double> probabilityFunction, final RandomGenerator rng) {
        final List<Pair<E, Double>> pmf = choices.stream()
                .map(e -> new Pair<>(e, probabilityFunction.apply(e))).collect(Collectors.toList());
        return new EnumeratedDistribution<>(rng, pmf).sample();
    }

    /**
     *  Return an array with column sums in each entry
     *
     * @param matrix Never {@code null}
     * @return Never {@code null}
     */
    public static double[] columnSums(final RealMatrix matrix) {
        Utils.nonNull(matrix);

        return IntStream.range(0, matrix.getColumnDimension())
                .mapToDouble(c -> MathUtils.sum(matrix.getColumn(c))).toArray();
    }

    /**
     *  Return an array with row sums in each entry
     *
     * @param matrix Never {@code null}
     * @return Never {@code null}
     */
    public static double[] rowSums(final RealMatrix matrix) {
        Utils.nonNull(matrix);

        return IntStream.range(0, matrix.getRowDimension())
                .mapToDouble(r -> MathUtils.sum(matrix.getRow(r))).toArray();
    }

    /**
     *  Return sum of 3d array
     *
     * @param array Never {@code null}
     * @return sum of array
     */
    public static double sum(final double[][][] array) {
        Utils.nonNull(array);
        double result = 0;
        for (double[][] d: array) {
            for (double[] da: d){
                for (double daa: da) {
                    result += daa;
                }
            }

        }
        return result;
    }
    
    public static int minIndex(final int ... values) {
        Utils.nonNull(values);
        if (values.length == 0) {
            return -1;
        }
        int minValue = values[0];
        int minIndex = 0;
        for (int i = 0; i < values.length; i++) {
            final int nextValue = values[i];
            if (nextValue < minValue) {
                minValue = nextValue;
                minIndex = i;
            }
        }
        return minIndex;
    }

    /**
     * Given an int array returns the difference between the second smallest element
     * and the smallest element.
     * @param values the input values.
     * @param defaultValue value to return in case the input array has zero or only one element.
     * @return 0 if the input array has less than 2 elements, otherwise the second smallest
     * - the smallest.
     * @throws IllegalArgumentException if {@code values} is {@code null}.
     */
    public static int secondSmallestMinusSmallest(final int[] values, final int defaultValue) {
        Utils.nonNull(values);
        if (values.length <= 1) {
            return defaultValue;
        } else {
            int smallest = values[0];
            int secondSmallest = Integer.MAX_VALUE;
            for (int i = 1; i < values.length; i++) {
                if (values[i] < smallest) {
                    secondSmallest = smallest;
                    smallest = values[i];
                } else if (values[i] < secondSmallest) {
                    secondSmallest = values[i];
                }
            }
            return secondSmallest - smallest;
        }
    }

    /**
     * Returns the smallest power of 2 that exceeds or equals {@code val}
     * @param val input value
     * @return power of 2 integer
     * @throws IllegalArgumentException if an int overlow is encountered
     */
    public static int smallestPowerOfTwoGreaterThan(final int val)
            throws IllegalArgumentException {
        if (val > Integer.MAX_VALUE/2 || val < 0)
            throw new IllegalArgumentException("The mallest power of 2 greater than " + val + " is greater than Integer.MAX_VALUE" +
                    " or negative input.");
        return val > 1 ? Integer.highestOneBit(val - 1) << 1 : 1;
    }

    /**
     * Shrink or expand a double array uniformly using nearest neighbor interpolation
     * @param data original array
     * @param newLength length of the new array
     * @return interpolated array
     * @throws IllegalArgumentException if the input array is empty
     */
    public static double[] nearestNeighborUniform1DInterpolate(@Nonnull final double[] data, final int newLength)
            throws IllegalArgumentException {
        if (data.length == 0)
            throw new IllegalArgumentException("The input array is empty.");
        ParamUtils.isPositive(newLength - 1, "The new length of the array must be >= 2");
        final double fact = (double)(data.length - 1)/(newLength - 1);
        return IntStream.range(0, newLength).mapToDouble(i -> data[(int) FastMath.floor(i*fact + 0.5)]).toArray();
    }


    /**
     * Find the maximum difference between entries of two arrays.  This is useful for testing convergence, for example
     */
    public static double maxDifference(final List<Double> array1, final List<Double> array2) {
        Utils.validateArg(array1.size() == array2.size(), "arrays must have same length.");
        Utils.validateArg(array1.size() > 0, "arrays must be non-empty");
        return IntStream.range(0, array1.size()).mapToDouble(n -> Math.abs(array1.get(n) - array2.get(n))).max().getAsDouble();
    }

    public static int median(final int[] values) {
        return (int) FastMath.round(new Median().evaluate(Arrays.stream(values).mapToDouble(n -> n).toArray()));
    }
}
