package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;

import java.util.List;

/**
 * Model prior for heterozygous pileups
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */

public abstract class HeterozygousPileupPriorModel {

    /**
     * Calculate the log likelihood for a pileup being heterzygous given a fixed allele fraction.
     *
     * Entry k in the pileup gives an additive contribution l_k,
     *
     *    l_k = log(alpha_k * alleleFraction + beta_k),
     *
     * to the log likelihood. The list of (alpha_k, beta_k) tuples are provided by {@code coeffs}. These
     * coefficients are calculated according to the read and mapping qualtities (see CNV-methods.pdf for details).
     *
     * Note that the full heterozygosity likelihood must be further adjusted according to the non-ref and non-alt
     * reads. For example, see {@link BayesianHetPulldownCalculator::getHetLogLikelihood}.
     *
     * @param alleleFraction the ref to alt allele fraction
     * @param coeffs list of (alpha, beta) tuples
     * @param minErrorProbability the theoretical minimum error probability for each nucleotide in the pileup
     * @return any double value.
     */
    protected double getHetLogLikelihoodFixedAlleleFraction(final double alleleFraction,
                                                            final List<ImmutablePair<Double, Double>> coeffs,
                                                            final double minErrorProbability) {
        return coeffs.stream()
                .mapToDouble(dat -> dat.getLeft() + alleleFraction * dat.getRight())
                .map(x -> GATKProtectedMathUtils.safeLog(x, minErrorProbability / 3))
                .sum();
    }

    /**
     * To be implemented by concrete models that extend this class. Concrete models weight
     * the fixed allele fraction likelihoods, as given by {@code getHetLogLikelihoodFixedAlleleFraction}, with
     * the implemented allele fraction prior.
     *
     * @param coeffs list of (alpha, beta) tuples
     * @param minErrorProbability the theoretical minimum error probability for each nucleotide in the pileup
     * @return
     */
    public abstract double getHetLogLikelihood(final List<ImmutablePair<Double, Double>> coeffs,
                                               final double minErrorProbability);
}
