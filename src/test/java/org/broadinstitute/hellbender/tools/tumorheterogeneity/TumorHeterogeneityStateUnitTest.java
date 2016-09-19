package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;

import static org.testng.Assert.*;

/**
 * Unit tests for {@link TumorHeterogeneityState}.  Checks that tests of state validity are correctly performed.
 */
public class TumorHeterogeneityStateUnitTest {
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSinglePopulation() {
        //fail if only one population (must have at least one variant and one normal)
        final double concentration = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Collections.singletonList(1.));
        final TumorHeterogeneityState.PopulationIndicators populationIndicators = new TumorHeterogeneityState.PopulationIndicators(Collections.singletonList(0));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.emptyList());
        new TumorHeterogeneityState(concentration, populationFractions, populationIndicators, variantProfiles);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testUnnormalizedPopulationFractions() {
        //fail if population fractions not normalized to unity
        final double concentration = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.1));
        final TumorHeterogeneityState.PopulationIndicators populationIndicators = new TumorHeterogeneityState.PopulationIndicators(Arrays.asList(0, 1));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.emptyList());
        new TumorHeterogeneityState(concentration, populationFractions, populationIndicators, variantProfiles);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNegativePopulationIndicators() {
        //fail if population indicators are negative
        final double concentration = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.9));
        final TumorHeterogeneityState.PopulationIndicators populationIndicators = new TumorHeterogeneityState.PopulationIndicators(Arrays.asList(0, -1));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.emptyList());
        new TumorHeterogeneityState(concentration, populationFractions, populationIndicators, variantProfiles);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInconsistentPopulationIndicators() {
        //fail if population indicators are inconsistent with total number of populations
        final double concentration = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.9));
        final TumorHeterogeneityState.PopulationIndicators populationIndicators = new TumorHeterogeneityState.PopulationIndicators(Collections.singletonList(2));
        final TumorHeterogeneityState.VariantProfile variantProfile = new TumorHeterogeneityState.VariantProfile(
                0.1,
                new TumorHeterogeneityState.VariantProfile.VariantIndicators(Collections.singletonList(true)),
                new TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators(Collections.singletonList(0)));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.singletonList(variantProfile));
        new TumorHeterogeneityState(concentration, populationFractions, populationIndicators, variantProfiles);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNoVariants() {
        //fail if number of variants is not positive
        final double concentration = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.9));
        final TumorHeterogeneityState.PopulationIndicators populationIndicators = new TumorHeterogeneityState.PopulationIndicators(Arrays.asList(0, 1));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.emptyList());
        new TumorHeterogeneityState(concentration, populationFractions, populationIndicators, variantProfiles);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadVariantSegmentFraction() {
        //fail if variant-segment fraction is not in [0, 1]
        final double concentration = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.9));
        final TumorHeterogeneityState.PopulationIndicators populationIndicators = new TumorHeterogeneityState.PopulationIndicators(Collections.singletonList(1));
        final TumorHeterogeneityState.VariantProfile variantProfile = new TumorHeterogeneityState.VariantProfile(
                -0.1,
                new TumorHeterogeneityState.VariantProfile.VariantIndicators(Collections.singletonList(true)),
                new TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators(Collections.singletonList(0)));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.singletonList(variantProfile));
        new TumorHeterogeneityState(concentration, populationFractions, populationIndicators, variantProfiles);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInconsistentNumbersOfPopulationsAndVariants() {
        //fail if number of variants + 1 is not equal to number of populations
        final double concentration = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.2, 0.7));
        final TumorHeterogeneityState.PopulationIndicators populationIndicators = new TumorHeterogeneityState.PopulationIndicators(Collections.singletonList(0));
        final TumorHeterogeneityState.VariantProfile variantProfile = new TumorHeterogeneityState.VariantProfile(
                0.1,
                new TumorHeterogeneityState.VariantProfile.VariantIndicators(Collections.singletonList(true)),
                new TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators(Collections.singletonList(0)));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.singletonList(variantProfile));
        new TumorHeterogeneityState(concentration, populationFractions, populationIndicators, variantProfiles);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDifferentNumberOfSegmentsAcrossVariants() {
        //fail if number of segments is not the same for all variants
        final double concentration = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.9));
        final TumorHeterogeneityState.PopulationIndicators populationIndicators = new TumorHeterogeneityState.PopulationIndicators(Collections.singletonList(0));
        final TumorHeterogeneityState.VariantProfile variantProfile1 = new TumorHeterogeneityState.VariantProfile(
                0.1,
                new TumorHeterogeneityState.VariantProfile.VariantIndicators(Collections.singletonList(true)),
                new TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators(Collections.singletonList(0)));
        final TumorHeterogeneityState.VariantProfile variantProfile2 = new TumorHeterogeneityState.VariantProfile(
                0.3,
                new TumorHeterogeneityState.VariantProfile.VariantIndicators(Arrays.asList(true, false)),
                new TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators(Arrays.asList(0, 1)));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Arrays.asList(variantProfile1, variantProfile2));
        new TumorHeterogeneityState(concentration, populationFractions, populationIndicators, variantProfiles);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDifferentNumberOfSegmentsWithinVariant() {
        //fail if number of segments is not the same for variant and ploidy-state indicators
        final double concentration = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.9));
        final TumorHeterogeneityState.PopulationIndicators populationIndicators = new TumorHeterogeneityState.PopulationIndicators(Collections.singletonList(0));
        final TumorHeterogeneityState.VariantProfile variantProfile = new TumorHeterogeneityState.VariantProfile(
                0.1,
                new TumorHeterogeneityState.VariantProfile.VariantIndicators(Collections.singletonList(true)),
                new TumorHeterogeneityState.VariantProfile.VariantPloidyStateIndicators(Arrays.asList(0, 1)));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.singletonList(variantProfile));
        new TumorHeterogeneityState(concentration, populationFractions, populationIndicators, variantProfiles);
    }
}