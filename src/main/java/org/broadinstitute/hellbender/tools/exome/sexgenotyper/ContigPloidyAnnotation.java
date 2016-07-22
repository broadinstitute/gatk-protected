package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import org.broadinstitute.hellbender.exceptions.UserException;

import javax.annotation.Nonnull;
import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * A class for storing contig ploidy annotations for sex genotyping
 * (see {@link TargetCoverageSexGenotypeCalculator)}).
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class ContigPloidyAnnotation {

    private final String contigName;

    private final ContigClass contigClass;

    private final Map<String, Integer> ploidyMap;

    private final Set<String> genotypesSet;

    /**
     * Public constructor
     *
     * @param contigName contig name
     * @param contigClass contig class (an instance of {@link ContigClass})
     * @param ploidyMap a map from genotype class identifier to the contig ploidy within that genotype class
     */
    public ContigPloidyAnnotation(@Nonnull final String contigName,
                                  @Nonnull final ContigClass contigClass,
                                  @Nonnull final Map<String, Integer> ploidyMap) {
        this.contigName = contigName;
        this.contigClass = contigClass;
        if (ploidyMap.isEmpty()) {
            throw new UserException.BadInput("The ploidy map must have at least one element");
        }
        /* if autosomal, check if all ploidy classes have the same value */
        if (contigClass == ContigClass.AUTOSOMAL) {
            if (ploidyMap.values().stream().collect(Collectors.toSet()).size() > 1) {
                throw new UserException.BadInput("Autosomal contigs must have the same ploidy for all sexes. Erroneously " +
                        "annotated contig: " + contigName);
            }
            if (ploidyMap.values().stream().filter(p -> p <= 0).count() > 0) {
                throw new UserException.BadInput("Autosomal contigs must have positive ploidies");
            }
        }

        /* if allosomal, check if all ploidy classes have >=0 values */
        if (contigClass == ContigClass.ALLOSOMAL) {
            if (ploidyMap.values().stream().filter(p -> p < 0).count() > 0) {
                throw new UserException.BadInput("Autosomal contigs must have non-negative ploidies");
            }
        }

        this.ploidyMap = Collections.unmodifiableMap(ploidyMap);
        this.genotypesSet = ploidyMap.keySet();
    }

    public String getContigName() {
        return contigName;
    }

    public ContigClass getContigClass() {
        return contigClass;
    }

    /**
     * Returns the ploidy of the contig for a genotype class
     * @param sexGenotype the string identifier of a genotype class
     * @return integer contig ploidy
     */
    public int getPloidy(@Nonnull final String sexGenotype) {
        if (!genotypesSet.contains(sexGenotype)) {
            throw new IllegalArgumentException("The ploidy class \"" + sexGenotype + "\" can not be found");
        } else {
            return ploidyMap.get(sexGenotype);
        }
    }

    /**
     * Returns the set of string identifiers for all of the ploidy-annotated genotypes
     * @return a set of strings
     */
    public Set<String> getGenotypesSet() {
        return Collections.unmodifiableSet(genotypesSet);
    }
}
