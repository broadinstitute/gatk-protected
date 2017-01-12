package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.buffer.DataBuffer;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;
import org.nd4j.linalg.ops.transforms.Transforms;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

/**
 * This test is to ensure that the DType is correctly set to Double
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class Nd4jUnitTest extends BaseTest {
    private static final double EPS = 1e-10;

    @Test
    public void testDType() {
        Assert.assertTrue(Nd4j.dataType().equals(DataBuffer.Type.DOUBLE), "Data type for Nd4j must be set to" +
                " double; otherwise, coverage model EM algorithm will not function properly");
    }

    /**
     * A test for implementing a multiplication like:
     *
     *      X_{b c} = \sum_{a} W_{a b c} v_{a}
     *
     * using matrix products and successive reshapes.
     */
    @Test
    public void testMulTensorVector() {
        /* generate random data */
        final int A = 5;
        final int B = 6;
        final int C = 7;
        final INDArray W = Nd4j.rand(new int[] {A, B, C});
        final INDArray v = Nd4j.rand(new int[] {A, 1});

        /* result using reshapes and matrix products */
        final INDArray X = W.reshape(new int[] {A, B*C}).transpose().mmul(v).reshape(new int[] {B, C});

        /* check against brute force result */
        for (int b = 0; b < B; b++) {
            for (int c = 0; c < C; c++) {
                double prod = 0;
                for (int a = 0; a < A; a++) {
                    prod += W.get(NDArrayIndex.point(a), NDArrayIndex.point(b), NDArrayIndex.point(c)).getDouble(0) *
                            v.getDouble(a);
                }
                Assert.assertEquals(X.getScalar(b, c).getDouble(0), prod, EPS);
            }
        }
    }

    /**
     * A test for implementing a multiplication like:
     *
     *      X_{a} = \sum_{b, c} W_{a b c} V_{b c}
     */
    @Test
    public void testMulTensorMatrix() {
        /* generate random data */
        final int A = 5;
        final int B = 6;
        final int C = 7;
        final INDArray W = Nd4j.rand(new int[] {A, B, C});
        final INDArray V = Nd4j.rand(new int[] {B, C});

        /* result using reshapes and matrix products */
        final INDArray X = W.reshape(new int[] {A, B*C}).mmul(V.reshape(new int[] {B*C, 1}));

        /* check against brute force result */
        for (int a = 0; a < A; a++) {
            double prod = 0;
            for (int b = 0; b < B; b++) {
                for (int c = 0; c < C; c++) {
                    prod += W.get(NDArrayIndex.point(a), NDArrayIndex.point(b), NDArrayIndex.point(c)).getDouble(0) *
                            V.get(NDArrayIndex.point(b), NDArrayIndex.point(c)).getDouble(0);
                }
            }
            Assert.assertEquals(X.getScalar(a).getDouble(0), prod, EPS);
        }
    }

//    @Test
//    public void testStability() {
//        final int N = 10000;
//        final int m = 2;
//        final int seed = 1234;
//        final RandomGenerator rg = new MersenneTwister(seed);
//        final double[] rand = IntStream.range(0, N*N).mapToDouble(i -> 1 - 2*rg.nextDouble()).toArray();
//        final INDArray arr = Nd4j.create(rand, new int[] {N, N});
//        final long t0 = System.nanoTime();
//        for (int k = 0; k < m; k++) {
//            final RealMatrix mat = Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(arr);
//            final INDArray arr2 = Nd4jApacheAdapterUtils.convertApacheMatrixToINDArray(mat);
//            arr.get(NDArrayIndex.all()).assign(arr.mmul(arr2));
//        }
//        final long t1 = System.nanoTime();
//        System.out.println(arr.sumNumber());
//        System.out.println((double)(t1-t0)/1000000000);
//    }

    @Test
    public void testSpeed() {
        final int[] sizes = new int[]{1, 10, 100, 1000, 10000, 100000, 10000000};
        final int numTrials = 500;
        final RandomGenerator rng = new MersenneTwister();
        final List<DescriptiveStatistics> apacheStats = new ArrayList<>(sizes.length);
        final List<DescriptiveStatistics> nd4jStats = new ArrayList<>(sizes.length);
        final List<DescriptiveStatistics> nd4jCreationStats = new ArrayList<>(sizes.length);

        for (int idx = 0; idx < sizes.length; idx++) {
            final DescriptiveStatistics currentApacheStats = new DescriptiveStatistics();
            final DescriptiveStatistics currentNd4jStats = new DescriptiveStatistics();
            final DescriptiveStatistics currentNd4jCreationStats = new DescriptiveStatistics();
            apacheStats.add(currentApacheStats);
            nd4jStats.add(currentNd4jStats);
            nd4jCreationStats.add(currentNd4jCreationStats);

            for (int n = 0; n < numTrials; n++) {
                final double[] vals = IntStream.range(0, sizes[idx])
                        .mapToDouble(i -> 100 * FastMath.abs(rng.nextDouble())).toArray();

                final long t0 = System.nanoTime();
                apacheLog(vals);
                final long t1 = System.nanoTime();
                ndLog(vals);
                final long t2 = System.nanoTime();
                ndJustCreate(vals);
                final long t3 = System.nanoTime();

                currentApacheStats.addValue((t1 - t0) / 1000000.0);
                currentNd4jStats.addValue((t2 - t1) / 1000000.0);
                currentNd4jCreationStats.addValue((t3 - t2) / 1000000.0);
            }
        }

        for (int idx = 0; idx < sizes.length; idx++) {
            System.out.println(String.format("N = %d, ApacheFastMath = %f +/- %f ms, Nd4jLog = %f +/- %f ms, Nd4jOverhead = %f +/- %f ms",
                    sizes[idx],
                    apacheStats.get(idx).getMean(), apacheStats.get(idx).getStandardDeviation(),
                    nd4jStats.get(idx).getMean(), nd4jStats.get(idx).getStandardDeviation(),
                    nd4jCreationStats.get(idx).getMean(), nd4jCreationStats.get(idx).getStandardDeviation()));
        }

    }

    private double[] apacheLog(final double[] vals) {
        return Arrays.stream(vals).parallel().map(FastMath::log).toArray();
    }

    private double[] ndLog(final double[] vals) {
        return Transforms.log(Nd4j.create(vals), false).data().asDouble();
    }

    private double[] ndJustCreate(final double[] vals) {
        return Nd4j.create(vals).data().asDouble();
    }

}
