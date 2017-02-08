package org.broadinstitute.hellbender.tools.exome.allelefractionrevision;

import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collection;
import java.util.stream.IntStream;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;
import static org.broadinstitute.hellbender.utils.MathUtils.log10Factorial;
import static org.broadinstitute.hellbender.utils.MathUtils.log10ToLog;

/**
 * Contains likelihood methods for the allele-fraction model.
 * See docs/CNVs/CNV-methods.pdf for a thorough description of the model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionLikelihoods {
    private static final double LN_2 = Math.log(2.);

    private AlleleFractionLikelihoods() {}

    /**
     * Compute the log-likelihood of a alt reads and r ref reads given minor fraction f and gamma hyperparameters
     * (specifying the distribution on allelic biases) mu (mean) and beta (rate = mean/variance) and given
     * an alt minor, ref minor, or outlier indicator state.  Note that this is a partially collapsed log-likelihood in that the
     * latent variable corresponding to the allelic bias at this site has been marginalized out but the indicator
     * variable has not been marginalized out.
     * <p>
     * See docs/CNVs/CNV-methods.pdf for derivation.
     * <p>
     * Finally, note that this is a static method and does not get mu, beta, and minorFraction from an AlleleFractionState object
     * We need such functionality because MCMC evaluates the likelihood under proposed parameter changes.
     *
     * @param parameters global parameters mean, variance, and outlier probability of allele fraction model
     * @param minorFraction minor allele fraction of segment containing this het site
     * @param count AllelicCount of alt and ref reads
     * @param indicator the hidden state (alt minor / ref minor / outlier)
     *
     * @return if indicator == ALT_MINOR:
     * <p>
     * log { [beta^alpha / Gamma(alpha)][(1-pi)/2] * int_{0 to infty} f^a * (1-f)^r * lambda^(alpha + r - 1) * exp(-beta*lambda)/(f + (1-f)*lambda)^n d lambda }
     * <p>
     * if indicator == REF_MINOR, same as ALT_MINOR but with f <--> 1 - f
     * <p>
     * if indicator == OUTLIER log {pi * a!r!/(n+1)!}
     * <p>
     * where alpha = mu*beta and n = a + r
     */
    public static double hetLogLikelihood(final AlleleFractionGlobalParameters parameters, final double minorFraction, final AllelicCount count, final AlleleFractionIndicator indicator) {
        return hetLogLikelihood(parameters, minorFraction, count, indicator, AllelicPanelOfNormals.EMPTY_PON);
    }

    /**
     * Computes {@link AlleleFractionLikelihoods#hetLogLikelihood} using allelic-bias priors derived from an
     * {@link AllelicPanelOfNormals}.  See docs/CNVs/CNV-methods.pdf for details.
     */
    public static double hetLogLikelihood(final AlleleFractionGlobalParameters parameters, final double minorFraction,
                                          final AllelicCount count, final AlleleFractionIndicator indicator,
                                          final AllelicPanelOfNormals allelicPoN) {
        final SimpleInterval site = count.getInterval();
        final double alpha = allelicPoN.equals(AllelicPanelOfNormals.EMPTY_PON) ? parameters.getAlpha() : allelicPoN.getAlpha(site);
        final double beta = allelicPoN.equals(AllelicPanelOfNormals.EMPTY_PON) ? parameters.getBeta() : allelicPoN.getBeta(site);
        final double pi = parameters.getOutlierProbability();
        final int a = count.getAltReadCount();
        final int r = count.getRefReadCount();
        return hetLogLikelihood(alpha, beta, minorFraction, pi, indicator, a, r);
    }

    private static double hetLogLikelihood(final double alpha, final double beta, final double minorFraction,
                                           final double outlierProbability, final AlleleFractionIndicator indicator,
                                           final int a, final int r) {
        if (indicator == AlleleFractionIndicator.OUTLIER) {
            return log(outlierProbability) + LN_2;  //outlier distribution is flat in [0, 0.5]
        } else {
            final double fAlt = indicator == AlleleFractionIndicator.ALT_MINOR ? minorFraction : 1 - minorFraction;
            return log((1 - outlierProbability) / 2) + logIntegralOverAllelicBias(alpha, beta, fAlt, a, r);
        }
    }

    /**
     * the log likelihood summed (marginalized) over indicator states, which we use in the fully collapsed model
     * in which latent variables (bias and indicator) are marginalized out
     *
     * @param parameters global parameters mean, variance, and outlier probability of allele fraction model
     * @param minorFraction minor allele fraction of segment containing this het site
     * @param count AllelicCount of alt and ref reads
     * @param allelicPoN allelic panel of normals constructed from total alt and ref counts observed across all normals at each site
     * @return the log of the likelihood at this het site, marginalized over indicator states.
     */
    public static double collapsedHetLogLikelihood(final AlleleFractionGlobalParameters parameters, final double minorFraction,
                                                   final AllelicCount count, final AllelicPanelOfNormals allelicPoN) {
        return GATKProtectedMathUtils.logSumExp(
                hetLogLikelihood(parameters, minorFraction, count, AlleleFractionIndicator.ALT_MINOR, allelicPoN),
                hetLogLikelihood(parameters, minorFraction, count, AlleleFractionIndicator.REF_MINOR, allelicPoN),
                hetLogLikelihood(parameters, minorFraction, count, AlleleFractionIndicator.OUTLIER, allelicPoN));
    }

    /**
     * the log-likelihood of all the hets in a segment
     *
     * @param parameters global parameters mean, variance, and outlier probability of allele fraction model
     * @param minorFraction minor allele fraction of segment containing this het site
     * @param counts AllelicCount of alt and ref reads in this segment
     * @param allelicPoN allelic panel of normals constructed from total alt and ref counts observed across all normals at each site
     * @return the sum of log-likelihoods over all het sites in a segment
     */
    public static double segmentLogLikelihood(final AlleleFractionGlobalParameters parameters, final double minorFraction,
                                              final Collection<AllelicCount> counts, final AllelicPanelOfNormals allelicPoN) {
        return counts.stream().mapToDouble(c -> collapsedHetLogLikelihood(parameters, minorFraction, c, allelicPoN)).sum();
    }

    /**
     * the total log likelihood of all segments
     * @param parameters parameters
     * @param data data
     * @return sum of log likelihoods of all segments
     */
    public static double logLikelihood(final AlleleFractionGlobalParameters parameters, final AlleleFractionState.MinorFractions minorFractions,
                                       final AlleleFractionData data) {
        return IntStream.range(0, data.getNumSegments())
                .mapToDouble(s -> segmentLogLikelihood(parameters, minorFractions.get(s), data.getCountsInSegment(s), data.getPoN())).sum();
    }

    protected static double logIntegralOverAllelicBias(final double sampleDepth, final double sampleBias, final double fAlt, final int a, final int r) {

    }
}
