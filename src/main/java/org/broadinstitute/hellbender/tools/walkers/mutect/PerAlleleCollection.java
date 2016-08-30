package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;

import java.util.*;

/**
 * A container for allele to value mapping.
 *
 * Each PerAlleleCollection may hold a value for each ALT allele and, optionally, a value for the REF allele.
 * For example,
 *
 *   PerAlleleCollection<Double> alleleFractions = PerAlleleCollection.createPerAltAlleleCollection()
 *
 * may be a container for allele fractions for ALT alleles in a variant context. While
 *
 *   PerAlleleCollection<Double> alleleCount = PerAlleleCollection.createPerRefAndAltAlleleCollection()
 *
 * may hold the allele counts for the REF allele and all ALT alleles in a variant context.
 *
 *
 **/
public class PerAlleleCollection<X> {
    // TODO: consider using Optional for ref allele
    private Optional<Allele> refAllele;
    private Optional<X> refValue;
    private Map<Allele, X> altAlleleValueMap;
    private Type type;

    public enum Type {ALT_ONLY, REF_AND_ALT}

    public PerAlleleCollection(final Type type){
        this.type = type;
        this.altAlleleValueMap = new HashMap<>();
        this.refAllele = Optional.empty();

    }

    /**
     * Take an allele, REF or ALT, and update its value appropriately
     *
     * @param allele : REF or ALT allele
     * @param newValue :
     */
    public void set(Allele allele, X newValue){
        if (allele == null || newValue == null){
            throw new IllegalArgumentException("allele or newValue is null");
        }
        if (allele.isReference() && type == Type.ALT_ONLY){
            throw new IllegalArgumentException("Collection stores values for alternate alleles only");
        }
        if (allele.isReference()){
            this.setRef(allele, newValue);
        } else {
            this.setAlt(allele, newValue);
        }
    }

    public void setRef(Allele refAllele, X newValue){
        if (refAllele == null || newValue == null){
            throw new IllegalArgumentException("refAllele or newValue is null");
        }
        if (refAllele.isNonReference()){
            throw new IllegalArgumentException("Setting Non-reference allele as reference");
        }

        if (this.refAllele.isPresent()){
            throw new IllegalArgumentException("Resetting the reference allele not permitted");
        }

        this.refAllele = Optional.of(refAllele);
        this.refValue = Optional.of(newValue);
    }

    public void setAlt(Allele altAllele, X newValue){
        if (altAllele == null || newValue == null){
            throw new IllegalArgumentException("altAllele or newValue is null");
        }
        if (altAllele.isReference()){
            throw new IllegalArgumentException("Setting reference allele as alt");
        }

        altAlleleValueMap.put(altAllele, newValue);
    }

    /**
     * Get the value for an allele, REF or ALT
     * @param allele
     */
    public X get(Allele allele){
        if (allele == null){
            throw new IllegalArgumentException("allele is null");
        }

        if (allele.isReference()){
            if (allele.equals(this.refAllele.get())){
                return(getRef());
            } else {
                throw new IllegalArgumentException("Requested ref allele does not match the stored ref allele");
            }
        } else {
            return(getAlt(allele));
        }
    }

    public X getRef(){
        if (type == Type.ALT_ONLY) {
            throw new IllegalStateException("Collection does not hold the REF allele");
        }

        if (this.refAllele.isPresent()){
            return(refValue.get());
        } else {
            throw new IllegalStateException("Collection's ref allele has not been set yet");
        }
    }

    public X getAlt(Allele allele){
        if (allele == null){
            throw new IllegalArgumentException("allele is null");
        }
        if (allele.isReference()){
            throw new IllegalArgumentException("allele is not an alt allele");
        }

        if (altAlleleValueMap.containsKey(allele)) {
            return(altAlleleValueMap.get(allele));
        } else {
            throw new IllegalArgumentException("Requested alt allele is not in the collection");
        }

    }


    public Set<Allele> getAltAlleles(){
        return(altAlleleValueMap.keySet());
    }
}