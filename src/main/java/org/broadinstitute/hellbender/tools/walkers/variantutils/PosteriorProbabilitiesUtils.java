package org.broadinstitute.hellbender.tools.walkers.variantutils;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.stream.Collectors;

public final class PosteriorProbabilitiesUtils {

    private PosteriorProbabilitiesUtils(){}

    /**
     * Calculates phred-scaled posterior probabilities for genotypes given the data and allele frequency priors.
     */
    public static VariantContext calculatePosteriorProbs(final VariantContext vc1,
                                                         final Collection<VariantContext> resources,
                                                         final int numRefSamplesFromMissingResources,
                                                         final double globalFrequencyPriorDirichlet,
                                                         final boolean useInputSamples,
                                                         final boolean useAC,
                                                         final boolean useACoff) {
        Utils.nonNull(vc1, "VariantContext vc1 is null");
        final Map<Allele,Integer> totalAlleleCounts = new HashMap<>();
        //only use discovered allele count if there are at least 10 samples
        final boolean useDiscoveredAC = !useACoff && vc1.getNSamples() >= 10;

        if(vc1.isSNP()) {
            //store the allele counts for each allele in the variant priors
            resources.forEach(r -> addAlleleCounts(totalAlleleCounts, r, useAC));

            //add the allele counts from the input samples (if applicable)
            if ( useInputSamples ) {
                addAlleleCounts(totalAlleleCounts,vc1,useAC);
            }

            //add zero allele counts for any reference alleles not seen in priors (if applicable)
            final int existingRefCounts = totalAlleleCounts.getOrDefault(vc1.getReference(), 0);
            totalAlleleCounts.put(vc1.getReference(), existingRefCounts + numRefSamplesFromMissingResources);
        }

        // now extract the counts of the alleles in vc1 in order
        final double[] alleleCounts = vc1.getAlleles().stream()
                .mapToDouble(a -> globalFrequencyPriorDirichlet + totalAlleleCounts.getOrDefault(a, 0)).toArray();

        final List<double[]> likelihoods = vc1.getGenotypes().stream().map(g -> parseLikelihoods(g)).collect(Collectors.toList());

        //TODO: for now just use priors that are SNPs because indel priors will bias SNP calls
        final boolean useFlatPriors = !vc1.isSNP() || (resources.isEmpty() && !useDiscoveredAC) || resources.stream().anyMatch(r -> !r.isSNP()) ;

        final List<double[]> posteriors = calculatePosteriorProbs(likelihoods,alleleCounts,vc1.getMaxPloidy(2), useFlatPriors);

        final GenotypesContext newContext = GenotypesContext.create();
        for ( int genoIdx = 0; genoIdx < vc1.getNSamples(); genoIdx ++ ) {
            final GenotypeBuilder builder = new GenotypeBuilder(vc1.getGenotype(genoIdx));
            builder.phased(vc1.getGenotype(genoIdx).isPhased());
            if ( posteriors.get(genoIdx) != null ) {
                GATKVariantContextUtils.makeGenotypeCall(vc1.getMaxPloidy(2), builder,
                        GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, posteriors.get(genoIdx), vc1.getAlleles());
                builder.attribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY,
                        Utils.listFromPrimitives(GenotypeLikelihoods.fromLog10Likelihoods(posteriors.get(genoIdx)).getAsPLs()));
            }
            newContext.add(builder.make());
        }

        final List<Integer> priors = Utils.listFromPrimitives(
                GenotypeLikelihoods.fromLog10Likelihoods(getDirichletPrior(alleleCounts, vc1.getMaxPloidy(2),useFlatPriors)).getAsPLs());

        final VariantContextBuilder builder = new VariantContextBuilder(vc1).genotypes(newContext).attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, priors);
        // add in the AC, AF, and AN attributes
        VariantContextUtils.calculateChromosomeCounts(builder, true);
        return builder.make();
    }

    private static double[] parseLikelihoods(Genotype genotype) {
        final Object PPfromVCF = genotype.getExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY);

        if (PPfromVCF == null){
            return getLikelihoodsVector(genotype);
        } else if (PPfromVCF instanceof String) {
            final String PPstring = (String) PPfromVCF;
            //samples not in trios will have PP tag like ".,.,." if family priors are applied
            return PPstring.charAt(0)=='.' ? getLikelihoodsVector(genotype) :
                    Arrays.stream(PPstring.split(",")).mapToDouble(s -> Double.parseDouble(s)/-10.0).toArray();
        } else {
            return Arrays.stream(extractInts(PPfromVCF)).mapToDouble(i -> i/-10.0).toArray();
        }
    }

    // return the double[] of likelihoods if available, otherwise null
    private static double[] getLikelihoodsVector(Genotype genotype) {
        return genotype.hasLikelihoods() ? genotype.getLikelihoods().getAsVector() : null;
    }

    /**
     * Given genotype likelihoods and known allele counts, calculate the posterior probabilities
     * over the genotype states
     * @param genotypeLikelihoods - the genotype likelihoods for the individual
     * @param knownAlleleCountsByAllele - the known allele counts in the population. For AC=2 AN=12 site, this is {10,2}
     * @param ploidy - the ploidy to assume
     * @param useFlatPriors - if true, apply flat priors to likelihoods in order to calculate posterior probabilities
     * @return - the posterior genotype likelihoods
     */
    protected static List<double[]> calculatePosteriorProbs(final List<double[]> genotypeLikelihoods,
                                                            final double[] knownAlleleCountsByAllele,
                                                            final int ploidy,
                                                            final boolean useFlatPriors) {
        Utils.validate( ploidy == 2, "Genotype posteriors not yet implemented for ploidy != 2");

        final double[] genotypePriorByAllele = getDirichletPrior(knownAlleleCountsByAllele,ploidy, useFlatPriors);
        final List<double[]> posteriors = new ArrayList<>(genotypeLikelihoods.size());
        for ( final double[] likelihoods : genotypeLikelihoods ) {
            double[] posteriorProbabilities = null;

            if ( likelihoods != null ) {
                Utils.validate( likelihoods.length == genotypePriorByAllele.length,
                        () -> String.format("Likelihoods not of correct size: expected %d, observed %d",
                            knownAlleleCountsByAllele.length*(knownAlleleCountsByAllele.length+1)/2,likelihoods.length));

                posteriorProbabilities = new double[genotypePriorByAllele.length];
                for ( int genoIdx = 0; genoIdx < likelihoods.length; genoIdx ++ ) {
                    posteriorProbabilities[genoIdx] = likelihoods[genoIdx] + genotypePriorByAllele[genoIdx];
                }

                posteriorProbabilities = MathUtils.normalizeLog10(posteriorProbabilities);

            }

            posteriors.add(posteriorProbabilities);
        }

        return posteriors;
    }

    // convenience function for a single genotypelikelihoods array. Just wraps.
    @VisibleForTesting
    static double[] calculatePosteriorProbs(final double[] genotypeLikelihoods,
                                            final double[] knownAlleleCountsByAllele,
                                            final int ploidy,
                                            final boolean useFlatPriors) {
        return calculatePosteriorProbs(Arrays.asList(genotypeLikelihoods),knownAlleleCountsByAllele,ploidy, useFlatPriors).get(0);
    }


    /**
     * Given known allele counts (whether external, from the sample, or both), calculate the prior distribution
     * over genotype states. This assumes
     *   1) Random sampling of alleles (known counts are unbiased, and frequency estimate is Dirichlet)
     *   2) Genotype states are independent (Hardy-Weinberg)
     * These assumptions give rise to a Dirichlet-Multinomial distribution of genotype states as a prior
     * (the "number of trials" for the multinomial is simply the ploidy)
     * @param knownCountsByAllele - the known counts per allele. For an AC=2, AN=12 site this is {10,2}
     * @param ploidy - the number of chromosomes in the sample. For now restricted to 2.
     * @return - the Dirichlet-Multinomial distribution over genotype states
     */
    @VisibleForTesting
    static double[] getDirichletPrior(final double[] knownCountsByAllele, final int ploidy, final boolean useFlatPrior) {
        Utils.validate( ploidy == 2, "Genotype priors not yet implemented for ploidy != 2");

        // multi-allelic format is
        // AA AB BB AC BC CC AD BD CD DD ...
        final double sumOfKnownCounts = MathUtils.sum(knownCountsByAllele);
        final double[] priors = new double[knownCountsByAllele.length*(knownCountsByAllele.length+1)/2];
        int priorIndex = 0;
        for ( int allele2 = 0; allele2 < knownCountsByAllele.length; allele2++ ) {
            for ( int allele1 = 0; allele1 <= allele2; allele1++) {
                if (useFlatPrior) {
                    priors[priorIndex++] = 1.0;
                } else {
                    final int[] counts = new int[knownCountsByAllele.length];
                    counts[allele1] += 1;
                    counts[allele2] += 1;
                    priors[priorIndex++] = MathUtils.dirichletMultinomial(knownCountsByAllele,counts);
                }
            }
        }

        return priors;
    }

    /**
     * Parse counts for each allele
     * @param counts - Map to store and return data
     * @param context - line to be parsed from the input VCF file
     * @param useAC - use allele count annotation value from VariantContext (vs. MLEAC)
     */
    private static void addAlleleCounts(final Map<Allele,Integer> counts, final VariantContext context, final boolean useAC) {
        final int[] ac;
        //use MLEAC value...
        if ( context.hasAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY) && ! useAC ) {
            ac = getAlleleCounts(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, context);
        }
        //...unless specified by the user in useAC or unless MLEAC is absent
        else if ( context.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) ) {
            ac = getAlleleCounts(VCFConstants.ALLELE_COUNT_KEY, context);
        }
        //if VariantContext annotation doesn't contain AC or MLEAC then get the data from direct evaluation
        else {
            ac = new int[context.getAlternateAlleles().size()];
            int idx = 0;
            for ( final Allele allele : context.getAlternateAlleles() ) {
                ac[idx++] = context.getCalledChrCount(allele);
            }
        }

        //since the allele count for the reference allele is not given in the VCF format,
        //calculate it from the allele number minus the total counts for alternate alleles
        for ( final Allele allele : context.getAlleles() ) {
            final int count;
            if ( allele.isReference() ) {
                if ( context.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY) ) {
                    count = Math.max(context.getAttributeAsInt(VCFConstants.ALLELE_NUMBER_KEY,-1) - (int) MathUtils.sum(ac),0); //occasionally an MLEAC value will sneak in that's greater than the AN
                } else {
                    count = Math.max(context.getCalledChrCount() - (int) MathUtils.sum(ac),0);
                }
            } else {
                count = ac[context.getAlternateAlleles().indexOf(allele)];
            }
            //if this allele isn't in the map yet, add it
            if ( ! counts.containsKey(allele) ) {
                counts.put(allele,0);
            }
            //add the count for the current allele to the existing value in the map
            counts.put(allele,count + counts.get(allele));
        }
    }

    /**
     * Retrieve allele count data from VariantContext using VCFkey, checks for correct number of values in VCF
     * @param VCFkey VariantContext annotation tag of interest (should be AC or MLEAC)
     * @param context VariantContext from which to extract the data
     * @return int[] with allele count data
     */
    private static int[] getAlleleCounts(final String VCFkey, final VariantContext context) {
        final Object alleleCountsFromVCF = context.getAttribute(VCFkey);
        if ( alleleCountsFromVCF instanceof List) {
            if ( ((List) alleleCountsFromVCF).size() != context.getAlternateAlleles().size() ) {
                throw new UserException(String.format("Variant does not contain the same number of MLE allele counts as alternate alleles for record at %s:%d", context.getContig(), context.getStart()));
            }
        }
        else if ( alleleCountsFromVCF instanceof String || alleleCountsFromVCF instanceof Integer) {//here length is 1
            if (context.getAlternateAlleles().size() != 1) {
                throw new UserException(String.format("Variant does not contain the same number of MLE allele counts as alternate alleles for record at %s:%d", context.getContig(), context.getStart()));
            }
        }
        return extractInts(alleleCountsFromVCF);
    }

    /**
     * Check the formatting on the Object returned by a call to VariantContext::getAttribute() and parse appropriately
     * @param integerListContainingVCField - Object returned by a call to VariantContext::getAttribute()
     * @return - array of ints
     *
     * //Note: if we're working with a integerListContainingVCField that's read directly out of the file it will be a String but
     * if it gets pulled from a VariantContext object built elsewhere in the code it will be an Integer or a List,
     */
    @SuppressWarnings("unchecked")
    public static int[] extractInts(final Object integerListContainingVCField) {
        List<Integer> mleList = null;
        if ( integerListContainingVCField instanceof List) {
            if ( ((List) integerListContainingVCField).get(0) instanceof String) {
                mleList = new ArrayList<>(((List) integerListContainingVCField).size());
                for ( final Object s : ((List)integerListContainingVCField)) {
                    mleList.add(Integer.parseInt((String) s));
                }
            } else {
                mleList = (List<Integer>) integerListContainingVCField;
            }
        } else if ( integerListContainingVCField instanceof Integer) {
            mleList = Arrays.asList((Integer) integerListContainingVCField);
        } else if ( integerListContainingVCField instanceof String) {
            mleList = Arrays.asList(Integer.parseInt((String)integerListContainingVCField));
        }
        Utils.nonNull( mleList, () -> String.format("VCF does not have properly formatted %s or %s.",
                    GATKVCFConstants.MLE_ALLELE_COUNT_KEY, VCFConstants.ALLELE_COUNT_KEY));

        final Integer firstElement = mleList.get(0);
        Utils.validate( firstElement instanceof Integer,
                () -> "BUG: The AC values should be an Integer, but was " + firstElement.getClass().getCanonicalName());

        return Ints.toArray(mleList);
    }
}
