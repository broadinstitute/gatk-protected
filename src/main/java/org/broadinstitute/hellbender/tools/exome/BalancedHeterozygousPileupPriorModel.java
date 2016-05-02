package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.lang3.tuple.ImmutablePair;

import java.util.List;

/**
 * Balanced model prior for heterozygous pileups, i.e. allele fraction = 1/2.
 *
 * This prior is suitable for detecting heterozygous sites when the ref to alt allele ratio is expected to be 1/2,
 * which is the case for sequenced reads from pure germline samples.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class BalancedHeterozygousPileupPriorModel extends HeterozygousPileupPriorModel {

    /**
     * Calculates the log likelihood of heterzygosity assuming allele fraction = 1/2
     * @param coeffs list of (alpha, beta) tuples
     * @param minErrorProbability the theoretical minimum error probability for each nucleotide in the pileup
     * @return
     */
    public double getHetLogLikelihood(final List<ImmutablePair<Double, Double>> coeffs,
                                      final double minErrorProbability) {
        return getHetLogLikelihoodFixedAlleleFraction(0.5, coeffs, minErrorProbability);
    }

}
