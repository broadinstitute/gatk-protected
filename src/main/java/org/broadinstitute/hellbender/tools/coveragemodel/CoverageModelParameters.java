package org.broadinstitute.hellbender.tools.coveragemodel;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.collect.Sets;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jIOUtils;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableReader;
import org.broadinstitute.hellbender.tools.exome.TargetWriter;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class reads, writes, and stores the coverage model parameters.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelParameters implements Serializable {

    private static final long serialVersionUID = -4350342293001054849L;

    private final List<Target> targetList;

    /* 1 x T */
    private final INDArray targetMeanLogBias;

    /* 1 x T */
    private final INDArray targetUnexplainedVariance;

    /* T x L */
    private final INDArray meanBiasCovariates;

    /* T x L x L */
    private final INDArray varBiasCovariates;

    /* 1 x L */
    private final INDArray biasCovariateARDCoefficients;

    private final int numTargets, numLatents;

    private final boolean ardEnabled, biasCovariatesEnabled;

    /**
     * Public constructor.
     *
     * Note:
     *
     * - If {@code meanBiasCovariates} and {@code varBiasCovariates} are both null, it is assumed that
     *   bias covariates are disabled.
     *
     * - If {@code biasCovariateARDCoefficients} is null, it is assumed that ARD for bias covariates is
     *   disabled.
     *
     * @param targetMeanLogBias target-specific mean log bias (m_t)
     * @param targetUnexplainedVariance target-specific unexplained bias variance (\Psi_t)
     * @param meanBiasCovariates mean bias covariates matrix (W_{t\mu})
     * @param varBiasCovariates bias covariates covariance tensor (cov[W_{t\mu}, W_{t\nu}])
     * @param biasCovariateARDCoefficients bias covariates ARD coefficients (\alpha_\mu)
     */
    public CoverageModelParameters(@Nonnull final List<Target> targetList,
                                   @Nonnull final INDArray targetMeanLogBias,
                                   @Nonnull final INDArray targetUnexplainedVariance,
                                   @Nullable final INDArray meanBiasCovariates,
                                   @Nullable final INDArray varBiasCovariates,
                                   @Nullable final INDArray biasCovariateARDCoefficients) {
        this.targetList = Utils.nonNull(targetList, "Target list must be non-null");
        this.targetMeanLogBias = Utils.nonNull(targetMeanLogBias, "Target-specific mean log bias must be non-null");
        this.targetUnexplainedVariance = Utils.nonNull(targetUnexplainedVariance, "Target-specific unexplained variance" +
                " must be non-null");
        Utils.validateArg((meanBiasCovariates == null && varBiasCovariates == null) ||
                (meanBiasCovariates != null && varBiasCovariates != null), "Either both mean and var of bias covariates" +
                " must be null, or both must be non-null");
        biasCovariatesEnabled = meanBiasCovariates != null && varBiasCovariates != null;
        ardEnabled = biasCovariateARDCoefficients != null;
        Utils.validateArg(biasCovariatesEnabled || !ardEnabled, "If bias covariates are disabled, ARD coefficients" +
                " must be null");
        this.meanBiasCovariates = meanBiasCovariates;
        this.varBiasCovariates = varBiasCovariates;
        this.biasCovariateARDCoefficients = biasCovariateARDCoefficients;

        this.numTargets = targetList.size();
        validateNDArrayShape(targetUnexplainedVariance, new int[] {1, numTargets}, "target-specific unexplained variance");
        validateNDArrayShape(targetMeanLogBias, new int[] {1, numTargets}, "target-specific mean log bias");
        if (biasCovariatesEnabled) {
            Utils.validateArg(meanBiasCovariates.rank() == 2, "The mean bias covariate NDArray must be rank 2");
            numLatents = meanBiasCovariates.size(1);
            validateNDArrayShape(meanBiasCovariates, new int[] {numTargets, numLatents}, "mean bias covariates");
            validateNDArrayShape(varBiasCovariates, new int[] {numTargets, numLatents, numLatents}, "var bias covariates");
            if (ardEnabled) {
                validateNDArrayShape(biasCovariateARDCoefficients, new int[] {1, numLatents}, "ARD coefficients");
            }
        } else {
            numLatents = 0;
        }
    }

    private void validateNDArrayShape(@Nonnull final INDArray arr, @Nonnull final int[] expected, @Nonnull final String name) {
        Utils.validateArg(Arrays.equals(expected, arr.shape()), "The shape of the provided NDArray for " + name +
                " is wrong; expected: " + Arrays.toString(expected) + ", provided: " + Arrays.toString(arr.shape()));
    }

    public boolean isARDEnabled() {
        return ardEnabled;
    }

    public boolean isBiasCovariatesEnabled() {
        return biasCovariatesEnabled;
    }

    public INDArray getTargetMeanLogBias() {
        return targetMeanLogBias;
    }

    public INDArray getTargetUnexplainedVariance() {
        return targetUnexplainedVariance;
    }

    public INDArray getMeanBiasCovariates() {
        return meanBiasCovariates;
    }

    public INDArray getVarBiasCovariates() {
        return varBiasCovariates;
    }

    public INDArray getBiasCovariateARDCoefficients() { return biasCovariateARDCoefficients; }

    public INDArray getTargetMeanBiasOnTargetBlock(@Nonnull final LinearSpaceBlock tb) {
        checkTargetBlock(tb);
        return targetMeanLogBias.get(NDArrayIndex.all(), NDArrayIndex.interval(tb.getBegIndex(), tb.getEndIndex()));
    }

    private void checkTargetBlock(@Nonnull LinearSpaceBlock tb) {
        ParamUtils.inRange(tb.getBegIndex(), 0, numTargets, "The begin index of target block is out of range");
        ParamUtils.inRange(tb.getEndIndex(), 0, numTargets, "The begin index of target block is out of range");
    }

    public INDArray getTargetUnexplainedVarianceOnTargetBlock(@Nonnull final LinearSpaceBlock tb) {
        checkTargetBlock(tb);
        return targetUnexplainedVariance.get(NDArrayIndex.all(), NDArrayIndex.interval(tb.getBegIndex(), tb.getEndIndex()));
    }

    public INDArray getMeanBiasCovariatesOnTargetBlock(@Nonnull final LinearSpaceBlock tb) {
        checkTargetBlock(tb);
        if (biasCovariatesEnabled) {
            return meanBiasCovariates.get(NDArrayIndex.interval(tb.getBegIndex(), tb.getEndIndex()),
                    NDArrayIndex.all());
        } else {
            return null;
        }
    }

    public INDArray getVarBiasCovariatesOnTargetBlock(@Nonnull final LinearSpaceBlock tb) {
        checkTargetBlock(tb);
        if (biasCovariatesEnabled) {
            return varBiasCovariates.get(NDArrayIndex.interval(tb.getBegIndex(), tb.getEndIndex()),
                    NDArrayIndex.all(), NDArrayIndex.all());
        } else {
            return null;
        }
    }

    private static void createOutputPath(final String outputPath) {
        final File outputPathFile = new File(outputPath);
        if (!outputPathFile.exists()) {
            if (!outputPathFile.mkdirs()) {
                throw new UserException.CouldNotCreateOutputFile(outputPathFile, "Could not create the output directory");
            }
        }
    }

    public int getNumTargets() {
        return numTargets;
    }

    public int getNumLatents() {
        return numLatents;
    }

    public List<Target> getTargetList() {
        return targetList;
    }

    private static double[] getNormalRandomNumbers(final int size, final double mean, final double std,
                                                   final RandomGenerator rng) {
        Utils.validateArg(std >= 0, "Standard deviation must be non-negative");
        return IntStream.range(0, size).mapToDouble(i -> mean + std * rng.nextGaussian()).toArray();
    }

    private static double[] getUniformRandomNumbers(final int size, final double min, final double max,
                                                    final RandomGenerator rng) {
        Utils.validateArg(max >= min, "Max value must be greater than min value");
        return IntStream.range(0, size).mapToDouble(i -> min + (max - min) * rng.nextDouble()).toArray();
    }

    /**
     * TODO
     *
     * @param targetList
     * @param numLatents
     * @param seed
     * @param randomMeanLogBiasStandardDeviation
     * @param randomBiasCovariatesStandardDeviation
     * @param randomMaxUnexplainedVariance
     * @param initialBiasCovariatesARDCoefficients
     * @return
     */
    public static CoverageModelParameters generateRandomModel(final List<Target> targetList,
                                                              final int numLatents,
                                                              final long seed,
                                                              final double randomMeanLogBiasStandardDeviation,
                                                              final double randomBiasCovariatesStandardDeviation,
                                                              final double randomMaxUnexplainedVariance,
                                                              final INDArray initialBiasCovariatesARDCoefficients) {
        Utils.validateArg(numLatents >= 0, "Dimension of the bias space must be non-negative");
        Utils.validateArg(randomBiasCovariatesStandardDeviation >= 0, "Standard deviation of random bias covariates" +
                " must be non-negative");
        Utils.validateArg(randomMeanLogBiasStandardDeviation >= 0, "Standard deviation of random mean log bias" +
                " must be non-negative");
        Utils.validateArg(randomMaxUnexplainedVariance >= 0, "Max random unexplained variance must be non-negative");

        final int numTargets = targetList.size();
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(seed));

        /* Gaussian random for mean log bias */
        final INDArray initialMeanLogBias = Nd4j.create(getNormalRandomNumbers(
                numTargets, 0, randomMeanLogBiasStandardDeviation, rng), new int[] {1, numTargets});

        /* Uniform random for unexplained variance */
        final INDArray initialUnexplainedVariance = Nd4j.create(getUniformRandomNumbers(
                numTargets, 0, randomMaxUnexplainedVariance, rng), new int[] {1, numTargets});

        final INDArray initialMeanBiasCovariates, initialVarBiasCovariates;

        if (numLatents > 0) {
            /* Gaussian random for bias covariates */
            initialMeanBiasCovariates = Nd4j.create(getNormalRandomNumbers(numTargets * numLatents, 0,
                    randomBiasCovariatesStandardDeviation, rng), new int[]{numTargets, numLatents});

            /* Zero variance of bias covariates */
            initialVarBiasCovariates = Nd4j.zeros(numTargets, numLatents, numLatents);
        } else {
            /* disable bias covariates: no mean and var bias covariates, no ARD */
            initialMeanBiasCovariates = null;
            initialVarBiasCovariates = null;
        }

        return new CoverageModelParameters(targetList, initialMeanLogBias, initialUnexplainedVariance,
                initialMeanBiasCovariates, initialVarBiasCovariates, initialBiasCovariatesARDCoefficients);
    }

    /**
     * Reads the model from disk
     *
     * TODO ARD, var log bias
     *
     * @param modelPath input model path
     * @return an instance of {@link CoverageModelParameters}
     */
    public static CoverageModelParameters read(@Nonnull final String modelPath) {
        final File modelPathFile = new File(Utils.nonNull(modelPath, "The input model path must be non-null"));
        Utils.validateArg(modelPathFile.exists(), "The model path does not exist: " + modelPathFile.getAbsolutePath());

        final File targetListFile = new File(modelPath, CoverageModelGlobalConstants.TARGET_LIST_OUTPUT_FILE);
        final List<Target> targetList;
        try (final Reader reader = new FileReader(targetListFile)) {
            targetList = TargetTableReader.readTargetFromReader(targetListFile.getAbsolutePath(), reader);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(targetListFile, "Could not read targets interval list");
        }

        final File targetMeanLogBiasFile = new File(modelPath, CoverageModelGlobalConstants.TARGET_MEAN_LOG_BIAS_OUTPUT_FILE);
        final INDArray targetMeanLogBias = Nd4jIOUtils.readNDArrayMatrixFromTextFile(targetMeanLogBiasFile);

        final File targetUnexplainedVarianceFile = new File(modelPath, CoverageModelGlobalConstants.TARGET_UNEXPLAINED_VARIANCE_OUTPUT_FILE);
        final INDArray targetUnexplainedVariance = Nd4jIOUtils.readNDArrayMatrixFromTextFile(targetUnexplainedVarianceFile);

        final File meanBiasCovariatesFile = new File(modelPath, CoverageModelGlobalConstants.MEAN_BIAS_COVARIATES_OUTPUT_FILE);
        final INDArray meanBiasCovariates = Nd4jIOUtils.readNDArrayMatrixFromTextFile(meanBiasCovariatesFile);

        final File varBiasCovariatesFile = new File(modelPath, CoverageModelGlobalConstants.VAR_BIAS_COVARIATES_OUTPUT_FILE);
        final INDArray varBiasCovariates = Nd4jIOUtils.readNDArrayTensorFromTextFile(varBiasCovariatesFile);

        final File biasCovariatesARDCoefficientsFile = new File(modelPath,
                CoverageModelGlobalConstants.BIAS_COVARIATES_ARD_COEFFICIENTS_OUTPUT_FILE);
        final INDArray biasCovariatesARDCoefficients = Nd4jIOUtils.readNDArrayMatrixFromTextFile(biasCovariatesARDCoefficientsFile);

        return new CoverageModelParameters(targetList, targetMeanLogBias, targetUnexplainedVariance,
                meanBiasCovariates, varBiasCovariates, biasCovariatesARDCoefficients);
    }

    /**
     * Writes the model to disk
     *
     * TODO ARD, var log bias
     *
     * @param outputPath model output path
     */
    public static void write(@Nonnull CoverageModelParameters model, @Nonnull final String outputPath) {
        /* create output directory if it doesn't exist */
        createOutputPath(Utils.nonNull(outputPath, "The output path string must be non-null"));

        /* write targets list */
        final File targetListFile = new File(outputPath, CoverageModelGlobalConstants.TARGET_LIST_OUTPUT_FILE);
        TargetWriter.writeTargetsToFile(targetListFile, model.getTargetList());

        final List<String> targetNames = model.getTargetList().stream()
                .map(Target::getName).collect(Collectors.toList());

        /* write target mean bias to file */
        final File targetMeanBiasFile = new File(outputPath, CoverageModelGlobalConstants.TARGET_MEAN_LOG_BIAS_OUTPUT_FILE);
        Nd4jIOUtils.writeNDArrayMatrixToTextFile(model.getTargetMeanLogBias(), targetMeanBiasFile, "MEAN_LOG_BIAS",
                null, targetNames);

        /* write target unexplained variance to file */
        final File targetUnexplainedVarianceFile = new File(outputPath, CoverageModelGlobalConstants.TARGET_UNEXPLAINED_VARIANCE_OUTPUT_FILE);
        Nd4jIOUtils.writeNDArrayMatrixToTextFile(model.getTargetUnexplainedVariance(), targetUnexplainedVarianceFile,
                "TARGET_UNEXPLAINED_VARIANCE", null, targetNames);

        if (model.isBiasCovariatesEnabled()) {
            /* write mean bias covariates to file */
            final List<String> meanBiasCovariatesNames = IntStream.range(0, model.getNumLatents())
                    .mapToObj(li -> String.format("BC_%d", li)).collect(Collectors.toList());
            final File meanBiasCovariatesFile = new File(outputPath, CoverageModelGlobalConstants.MEAN_BIAS_COVARIATES_OUTPUT_FILE);
            Nd4jIOUtils.writeNDArrayMatrixToTextFile(model.getMeanBiasCovariates(), meanBiasCovariatesFile,
                    "MEAN_BIAS_COVARIATES", targetNames, meanBiasCovariatesNames);

            /* write norm_2 of mean bias covariates to file */
            final double[] biasCovariatesNorm2 = new double[model.numLatents];
            final INDArray WTW = model.getMeanBiasCovariates().transpose().mmul(model.getMeanBiasCovariates());
            for (int li = 0; li < model.getNumLatents(); li++) {
                biasCovariatesNorm2[li] = WTW.getDouble(li, li);
            }
            final File biasCovariatesNorm2File = new File(outputPath, CoverageModelGlobalConstants.MEAN_BIAS_COVARIATES_NORM2_OUTPUT_FILE);
            Nd4jIOUtils.writeNDArrayMatrixToTextFile(Nd4j.create(biasCovariatesNorm2, new int[]{1, model.getNumLatents()}),
                    biasCovariatesNorm2File, "MEAN_BIAS_COVARIATES_NORM_2", null, meanBiasCovariatesNames);

            /* write var bias covariates to file */
            final File varBiasCovariatesFile = new File(outputPath, CoverageModelGlobalConstants.VAR_BIAS_COVARIATES_OUTPUT_FILE);
            Nd4jIOUtils.writeNDArrayTensorToTextFile(model.getVarBiasCovariates(), varBiasCovariatesFile,
                    "VAR_BIAS_COVARIATES", null);

            /* if ARD is enabled, write the ARD coefficients as well */
            if (model.isARDEnabled()) {
                final File biasCovariatesARDCoefficientsFile = new File(outputPath, CoverageModelGlobalConstants.BIAS_COVARIATES_ARD_COEFFICIENTS_OUTPUT_FILE);
                Nd4jIOUtils.writeNDArrayMatrixToTextFile(model.getBiasCovariateARDCoefficients(),
                        biasCovariatesARDCoefficientsFile, "BIAS_COVARIATES_ARD_COEFFICIENTS", null, meanBiasCovariatesNames);
            }
        }
    }

    /**
     * This method "adapts" a model to a read count collection in the following sense:
     *
     *     - removes targets that are not included in the model from the read counts collection
     *     - removes targets that are in the read count collection from the model
     *     - rearranges model targets in the same order as read count collection targets
     *
     * The modifications are not done in-plane and the original input parameters remain intact.
     *
     * @param model a model
     * @param readCounts a read count collection
     * @return a pair of model and read count collection
     */
    public static ImmutablePair<CoverageModelParameters, ReadCountCollection> adaptModelToReadCountCollection(
            @Nonnull final CoverageModelParameters model, @Nonnull final ReadCountCollection readCounts,
            @Nonnull final Logger logger) {
        logger.info("Adapting model to read counts...");
        Utils.nonNull(model, "The model parameters must be non-null");
        Utils.nonNull(readCounts, "The read count collection must be non-null");
        Utils.nonNull(logger, "The logger must be non-null");

        final List<Target> modelTargetList = model.getTargetList();
        final List<Target> readCountsTargetList = readCounts.targets();
        final Set<Target> mutualTargetList = Sets.intersection(new HashSet<>(modelTargetList),
                new HashSet<>(readCountsTargetList));
        final List<Target> finalTargetList = readCountsTargetList.stream()
                .filter(mutualTargetList::contains)
                .collect(Collectors.toList());
        final Set<Target> finalTargetsSet = new LinkedHashSet<>(finalTargetList);

        logger.info("Number of mutual targets: " + finalTargetList.size());
        Utils.validateArg(finalTargetList.size() > 0, "The intersection between model targets and targets from read count" +
                    " collection is empty. Please check there the model is compatible with the given read count" +
                    " collection.");

        if (modelTargetList.size() > finalTargetList.size()) {
            logger.info("The following targets dropped from the model: " + Sets.difference(new HashSet<>(modelTargetList),
                    finalTargetsSet).stream().map(Target::getName).collect(Collectors.joining(", ", "[", "]")));
        }

        if (readCountsTargetList.size() > finalTargetList.size()) {
            logger.info("The following targets dropped from read counts: " + Sets.difference(new HashSet<>(readCountsTargetList),
                    finalTargetsSet).stream().map(Target::getName).collect(Collectors.joining(", ", "[", "]")));
        }

        /* the targets in {@code subsetReadCounts} follow the original order of targets in {@code readCounts} */
        final ReadCountCollection subsetReadCounts = readCounts.subsetTargets(finalTargetsSet);

        /* fetch original model parameters */
        final INDArray originalModelTargetMeanBias = model.getTargetMeanLogBias();
        final INDArray originalModelTargetUnexplainedVariance = model.getTargetUnexplainedVariance();
        final INDArray originalModelMeanBiasCovariates = model.getMeanBiasCovariates();
        final INDArray originalModelVarBiasCovariates = model.getVarBiasCovariates();

        /* re-arrange targets */
        final int[] newTargetIndicesInOriginalModel = finalTargetList.stream()
                .mapToInt(modelTargetList::indexOf)
                .toArray();
        final INDArray newModelTargetMeanBias = Nd4j.create(new int[] {1, finalTargetList.size()});
        final INDArray newModelTargetUnexplainedVariance = Nd4j.create(new int[] {1, finalTargetList.size()});
        final INDArray newModelMeanBiasCovariates = Nd4j.create(new int[] {finalTargetList.size(),
                model.getNumLatents()});
        final INDArray newModelVarBiasCovariates = Nd4j.create(new int[] {finalTargetList.size(),
                model.getNumLatents(), model.getNumLatents()});
        IntStream.range(0, finalTargetList.size())
                .forEach(ti -> {
                    newModelTargetMeanBias.put(0, ti,
                            originalModelTargetMeanBias.getDouble(0, newTargetIndicesInOriginalModel[ti]));
                    newModelTargetUnexplainedVariance.put(0, ti,
                            originalModelTargetUnexplainedVariance.getDouble(0, newTargetIndicesInOriginalModel[ti]));
                    newModelMeanBiasCovariates.get(NDArrayIndex.point(ti), NDArrayIndex.all())
                            .assign(originalModelMeanBiasCovariates.get(NDArrayIndex.point(newTargetIndicesInOriginalModel[ti]),
                                    NDArrayIndex.all()));
                    newModelVarBiasCovariates.get(NDArrayIndex.point(ti), NDArrayIndex.all(), NDArrayIndex.all())
                            .assign(originalModelVarBiasCovariates.get(NDArrayIndex.point(newTargetIndicesInOriginalModel[ti]),
                                    NDArrayIndex.all(), NDArrayIndex.all()));
                });

        return ImmutablePair.of(new CoverageModelParameters(finalTargetList, newModelTargetMeanBias,
                newModelTargetUnexplainedVariance, newModelMeanBiasCovariates, newModelVarBiasCovariates,
                model.getBiasCovariateARDCoefficients()), subsetReadCounts);
    }

}
