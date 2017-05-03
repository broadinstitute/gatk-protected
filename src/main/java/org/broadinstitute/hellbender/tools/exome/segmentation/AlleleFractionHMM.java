package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionGlobalParameters;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionLikelihoods;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;

/**
 * The Allele Fraction Hidden Markov Model is a generative model describing latent CNV states
 * with different minor allele fractions and observed ref and alt allele counts at a set of
 * het sites.
 *
 * Hidden states are essentially minor allele fractions f, but for convenience we represent them
 * as integers 0, 1, . . . K - 1 and store minor allele fractions in a corresponding array.
 *
 * The model contains a memory length parameter representing the prior probability for the minor
 * allele fraction state to be "forgotten" as a function of distance d between consecutive SNPs --
 * the probability to remember a state is exp(-d/memoryLength).  If the state is forgotten, a new state
 * is chosen with probabilities given by an array of weights.
 *
 * Thus our transition probabilities are P(i -> j) = exp(-d/D) delta_{ij} + (1 - exp(-d/D)) weights[j]
 * where delta is the Kronecker delta and D is the memory length.
 *
 * Emission likelihoods as a function of minor allele fraction are defined and described in
 * {@link AlleleFractionLikelihoods}.  These likelihoods require this class to contain a set of
 * {@link AlleleFractionGlobalParameters} and an {@link AllelicPanelOfNormals}.
 *
 * This model, along with the likelihoods model, is described thoroughly in docs/CNVs/CNV-methods.pdf
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionHMM extends ClusteringGenomicHMM<AllelicCount, Double> {
    private final AlleleFractionGlobalParameters parameters;
    private final AllelicPanelOfNormals allelicPoN;

    /**
     * @param minorAlleleFractions array of minor allele fractions corresponding to the hidden states
     * @param weights array of (real-space, not log) prior probabilities of each hidden state
     *                when memory of the previous state is lost.  These may be unnormalized relative
     *                probabilities, which is useful when using variational Bayes.
     * @param memoryLength when consecutive SNPs are a distance d bases apart, the prior probability
     *                     for memory of the CNV state to be kept is exp(-d/memoryLength)
     * @param allelicPoN allelic panel of normals containing prior knowledge of allelic biases at common SNPs
     * @param parameters the global parameters of the allelic bias model: mean bias, bias variance, and
     *                   outlier probability
     */
    public AlleleFractionHMM(final List<Double> minorAlleleFractions, final List<Double> weights,
                             final double memoryLength, final AllelicPanelOfNormals allelicPoN,
                             final AlleleFractionGlobalParameters parameters) {
        super(minorAlleleFractions, weights, memoryLength);
        minorAlleleFractions.forEach(f -> ParamUtils.inRange(f, 0, 0.5, "minor fractions must be between 0 and 1/2, found " + f));
        this.allelicPoN = Utils.nonNull(allelicPoN);
        this.parameters = Utils.nonNull(parameters);
    }

    @Override
    public double logEmissionProbability(final AllelicCount data, final Integer state, final SimpleInterval position) {
        return logEmissionProbability(data, getHiddenStateValue(state));
    }

    /**
     * Visible for {@link AlleleFractionSegmenter}
     */
    @Override
    public double logEmissionProbability(final AllelicCount data, final Double minorFraction) {
        ParamUtils.inRange(minorFraction, 0.0, 0.5, "Minor fraction must be in [0, 1/2].");
        return logEmissionProbability(data, minorFraction, parameters, allelicPoN);
    }

    /**
     * Visible for {@link AlleleFractionSegmenter}
     */
    protected static double logEmissionProbability(final AllelicCount data, final double minorFraction,
                                                   final AlleleFractionGlobalParameters parameters, final AllelicPanelOfNormals allelicPoN) {
        return AlleleFractionLikelihoods.collapsedHetLogLikelihood(parameters, minorFraction, data, allelicPoN);
    }

    public double getMinorAlleleFraction(final int state) { return getHiddenStateValue(state); }
    public AllelicPanelOfNormals getAllelicPoN() { return allelicPoN; }
    public AlleleFractionGlobalParameters getParameters() { return parameters; }
}
