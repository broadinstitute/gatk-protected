package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents the panel of normals used for allele-bias correction.  See docs/CNVs/CNV-methods.pdf.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicPanelOfNormals {
    private static final Logger logger = LogManager.getLogger(AllelicPanelOfNormals.class);

    public static final AllelicPanelOfNormals EMPTY_PON = new AllelicPanelOfNormals();

    //local (per site) hyperparameter values
    private final Map<SimpleInterval, HyperparameterValues> siteToHyperparameterPairMap = new HashMap<>();

    //global hyperparameter values
    private final HyperparameterValues mleHyperparameterValues;
    private final double mleMeanBias;
    private final double mleBiasVariance;

    private static final String MLE_GLOBAL_ALPHA_COMMENT_STRING = "MLE_GLOBAL_ALPHA=";
    private static final String MLE_GLOBAL_BETA_COMMENT_STRING = "MLE_GLOBAL_BETA=";

    private AllelicPanelOfNormals() {
        mleHyperparameterValues = new HyperparameterValues(Double.NaN, Double.NaN);
        mleMeanBias = Double.NaN;
        mleBiasVariance = Double.NaN;
    }

    /**
     * Constructs an allelic panel of normals from a tab-separated file.
     * @param inputFile    contains lines specifying hyperparameters at each site, as well as global hyperparameters in comment lines
     */
    public AllelicPanelOfNormals(final File inputFile) {
        Utils.regularReadableUserFile(inputFile);

        try (final AllelicPanelOfNormalsReader reader = new AllelicPanelOfNormalsReader(inputFile)) {
            //parse comment lines for global hyperparameter values
            final List<String> commentLines = FileUtils.readLines(inputFile).stream().filter(l -> l.contains(TableUtils.COMMENT_PREFIX)).collect(Collectors.toList());
            final List<String> mleGlobalAlphaCommentLines = commentLines.stream().filter(l -> l.contains(MLE_GLOBAL_ALPHA_COMMENT_STRING)).collect(Collectors.toList());
            final List<String> mleGlobalBetaCommentLines = commentLines.stream().filter(l -> l.contains(MLE_GLOBAL_BETA_COMMENT_STRING)).collect(Collectors.toList());
            if (mleGlobalAlphaCommentLines.size() != 1 || mleGlobalBetaCommentLines.size() != 1) {
                throw new UserException.BadInput("MLE global hyperparameter values for allelic panel of normals were not specified correctly in comment lines of file " + inputFile.getAbsolutePath() + ".");
            }
            final double mleGlobalAlpha = Double.valueOf(mleGlobalAlphaCommentLines.get(0).replace(MLE_GLOBAL_ALPHA_COMMENT_STRING, ""));
            final double mleGlobalBeta = Double.valueOf(mleGlobalBetaCommentLines.get(0).replace(MLE_GLOBAL_BETA_COMMENT_STRING, ""));
            ParamUtils.isPositive(mleGlobalAlpha, "MLE global hyperparameter alpha must be positive.");
            ParamUtils.isPositive(mleGlobalBeta, "MLE global hyperparameter alpha must be positive.");
            mleHyperparameterValues = new HyperparameterValues(mleGlobalAlpha, mleGlobalBeta);
            mleMeanBias = meanBias(mleGlobalAlpha, mleGlobalBeta);
            mleBiasVariance = biasVariance(mleGlobalAlpha, mleGlobalBeta);
            //parse data lines for local hyperparameter values
            reader.stream().forEach(s -> siteToHyperparameterPairMap.put(s.getKey(), s.getValue()));
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        }
    }

    /**
     * Constructs an allelic panel of normals from an {@link AllelicCountCollection} that contains
     * total alt and ref counts observed across all normals at each site.
     * @param counts    total alt and ref counts observed across all normals at each site
     */
    public AllelicPanelOfNormals(final AllelicCountCollection counts) {
        mleHyperparameterValues = AlleleFractionInitializer.calculateMLEHyperparameterValues(counts);   //use AlleleFractionInitializer to fit MLE global hyperparameter values
        mleMeanBias = meanBias(mleHyperparameterValues.alpha, mleHyperparameterValues.beta);
        mleBiasVariance = biasVariance(mleHyperparameterValues.alpha, mleHyperparameterValues.beta);
        initializeSiteToHyperparameterPairMap(counts);
    }

    /**
     * Gets the reference-bias alpha hyperparameter at a given SNP site if it is in the panel of normals
     * and the MLE alpha hyperparameter across all sites if it is not.
     * @param site  SNP site
     * @return      reference-bias alpha hyperparameter if site is in panel of normals,
     *              MLE alpha hyperparameter across all sites if it is not
     */
    public double getAlpha(final SimpleInterval site) {
        throwExceptionIfPoNIsEmpty();
        Utils.nonNull(site);
        return siteToHyperparameterPairMap.getOrDefault(site, mleHyperparameterValues).alpha;
    }

    /**
     * Gets the reference-bias beta hyperparameter at a given SNP site if it is in the panel of normals
     * and the MLE beta hyperparameter across all sites if it is not.
     * @param site  SNP site
     * @return      reference-bias beta hyperparameter if site is in panel of normals,
     *              MLE beta hyperparameter across all sites if it is not
     */
    public double getBeta(final SimpleInterval site) {
        throwExceptionIfPoNIsEmpty();
        Utils.nonNull(site);
        return siteToHyperparameterPairMap.getOrDefault(site, mleHyperparameterValues).beta;
    }

    /**
     * Gets the MLE mean-bias hyperparameter across all sites.
     * @return  MLE mean-bias hyperparameter across all sites
     */
    public double getMLEMeanBias() {
        throwExceptionIfPoNIsEmpty();
        return mleMeanBias;
    }

    /**
     * Writes out the {@link AllelicPanelOfNormals} to the specified file.
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     */
    public void write(final File outputFile) {
        final List<Map.Entry<SimpleInterval, HyperparameterValues>> sortedMapEntries = collectSortedMapEntries(siteToHyperparameterPairMap);
        try (final AllelicPanelOfNormalsWriter writer = new AllelicPanelOfNormalsWriter(outputFile)) {
            writer.writeComment(MLE_GLOBAL_ALPHA_COMMENT_STRING + AllelicPanelOfNormalsWriter.formatDouble(mleHyperparameterValues.alpha));
            writer.writeComment(MLE_GLOBAL_BETA_COMMENT_STRING + AllelicPanelOfNormalsWriter.formatDouble(mleHyperparameterValues.beta));
            writer.writeAllRecords(sortedMapEntries);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    /**
     * Gets the MLE bias-variance hyperparameter across all sites.
     * @return  MLE bias-variance hyperparameter across all sites
     */
    public double getMLEBiasVariance() {
        throwExceptionIfPoNIsEmpty();
        return mleBiasVariance;
    }

    static class HyperparameterValues {
        private final double alpha;
        private final double beta;

        HyperparameterValues(final double alpha, final double beta) {
            this.alpha = alpha;
            this.beta = beta;
        }

        /**
         * Initializes the hyperparameter values at a site given the observed counts in all normals.
         * @param alpha global hyperparameter MLE value for alpha
         * @param beta  global hyperparameter MLE value for beta
         * @param a     total alt counts observed across all normals at site
         * @param r     total ref counts observed across all normals at site
         */
        private HyperparameterValues(final double alpha, final double beta, final int a, final int r) {
            final double f = 0.5;
            final int n = a + r;
            final double lambda0 = AlleleFractionLikelihoods.biasPosteriorMode(alpha, beta, f, a, r);
            final double kappa = AlleleFractionLikelihoods.biasPosteriorCurvature(alpha, f, r, n, lambda0);
            this.alpha = AlleleFractionLikelihoods.biasPosteriorEffectiveAlpha(lambda0, kappa);
            this.beta = AlleleFractionLikelihoods.biasPosteriorEffectiveBeta(lambda0, kappa);
        }

        double getAlpha() {
            return alpha;
        }

        double getBeta() {
            return beta;
        }
    }

    private void initializeSiteToHyperparameterPairMap(final AllelicCountCollection counts) {
        logger.info("Initializing allelic panel of normals...");
        for (final AllelicCount count : counts.getCounts()) {
            final SimpleInterval site = count.getInterval();
            final HyperparameterValues hyperparameterValues = new HyperparameterValues(
                    mleHyperparameterValues.alpha, mleHyperparameterValues.beta,
                    count.getAltReadCount(), count.getRefReadCount());
            if (siteToHyperparameterPairMap.containsKey(site)) {
                throw new UserException.BadInput("Input AllelicCountCollection for allelic panel of normals contains duplicate sites.");
            } else {
                siteToHyperparameterPairMap.put(site, hyperparameterValues);
            }
        }
        logger.info("Allelic panel of normals initialized.");
    }

    private void throwExceptionIfPoNIsEmpty() {
        if (this.equals(EMPTY_PON)) {
            throw new UnsupportedOperationException("Cannot get MLE hyperparameters for empty panel of normals.");
        }
    }

    private static double meanBias(final double alpha, final double beta) {
        return alpha / beta;
    }

    private static double biasVariance(final double alpha, final double beta) {
        return alpha / (beta * beta);
    }

    //returns a list of the sites in siteToHyperarameterPairMap that are sorted by SimpleInterval (contigs are sorted by lexicographical order)
    private static List<Map.Entry<SimpleInterval, HyperparameterValues>> collectSortedMapEntries(final Map<SimpleInterval, HyperparameterValues> siteToHyperparameterPairMap) {
        final List<SimpleInterval> sortedMapKeys = new ArrayList<>(siteToHyperparameterPairMap.keySet());
        Collections.sort(sortedMapKeys, IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
        return sortedMapKeys.stream().map(si -> new AbstractMap.SimpleEntry<>(si, siteToHyperparameterPairMap.get(si))).collect(Collectors.toList());
    }
}
