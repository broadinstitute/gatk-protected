package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

/**
 * Created by tsato on 6/21/16.
 */
public class PerAlleleCollectionUnitTest {

    @Test
    public void testSet() throws Exception {
        PerAlleleCollection<Integer> alleleCounts = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        Allele refA = Allele.create("A", true);
        Allele altT = Allele.create("T", false);
        alleleCounts.set(refA, 40);
        alleleCounts.set(altT, 10);
        assertEquals((int)alleleCounts.getRef(), 40);
        assertEquals((int)alleleCounts.getAlt(altT), 10);
    }

    @Test
    public void testGet() throws Exception {
        PerAlleleCollection<Integer> alleleCounts = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        Allele refA = Allele.create("A", true);
        Allele altT = Allele.create("T", false);
        alleleCounts.set(refA, 40);
        alleleCounts.set(altT, 10);
        assertEquals((int)alleleCounts.get(refA), 40);
        assertEquals((int)alleleCounts.get(altT), 10);
    }


    @Test
    public void testGetAltAlleles() throws Exception {
        PerAlleleCollection<Integer> alleleCounts = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        Allele altA = Allele.create("A", false);
        Allele altC = Allele.create("C", false);
        Allele altG = Allele.create("G", false);
        Allele altT = Allele.create("T", false);
        Allele[] altAlleles = {altA, altC, altG, altT};
        for (Allele altAllele : altAlleles ) {
            alleleCounts.set(altAllele, 3);
        }

        for (Allele altAllele : altAlleles ) {
            assertTrue(alleleCounts.getAltAlleles().contains(altAllele));
        }

        assertFalse(alleleCounts.getAltAlleles().contains(Allele.create("A", true)));
    }
}