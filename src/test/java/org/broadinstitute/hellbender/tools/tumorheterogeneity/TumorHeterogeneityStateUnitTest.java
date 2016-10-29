package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Unit tests for {@link TumorHeterogeneityState}.  Checks that tests of state validity are correctly performed.
 */
public class TumorHeterogeneityStateUnitTest {
    private static final double EPSILON = 1E-10;
    private static final PosteriorSummary DUMMY_POSTERIOR_SUMMARY = new PosteriorSummary(Double.NaN, Double.NaN, Double.NaN);
    private static final DecileCollection DUMMY_DECILE_COLLECTION =
            new DecileCollection(Collections.singletonList(Double.NaN), DecileCollection.ConstructionMode.SAMPLES);
    private static final PloidyState NORMAL_PLOIDY_STATE = new PloidyState(1, 1);
    private static final PloidyStatePrior VARIANT_PLOIDY_STATE_PRIOR;
    private static final double DUMMY_HYPERPARAMETER = 1.;
    private static final TumorHeterogeneityPriorCollection PRIORS;
    private static final TumorHeterogeneityState STATE;
    private static final TumorHeterogeneityData DATA;

    static {
        DUMMY_POSTERIOR_SUMMARY.setDeciles(DUMMY_DECILE_COLLECTION);
        final Map<PloidyState, Double> unnormalizedLogProbabilityMassFunctionMap = new LinkedHashMap<>();
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 0), 0.);
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(0, 1), 0.);
        unnormalizedLogProbabilityMassFunctionMap.put(NORMAL_PLOIDY_STATE, 0.);
        unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(1, 2), 0.);
        VARIANT_PLOIDY_STATE_PRIOR = new PloidyStatePrior(unnormalizedLogProbabilityMassFunctionMap);
        PRIORS = new TumorHeterogeneityPriorCollection(
                NORMAL_PLOIDY_STATE, VARIANT_PLOIDY_STATE_PRIOR, DUMMY_HYPERPARAMETER, DUMMY_HYPERPARAMETER, DUMMY_HYPERPARAMETER, DUMMY_HYPERPARAMETER, DUMMY_HYPERPARAMETER, DUMMY_HYPERPARAMETER);

        final double concentration = 1.;
        final double copyRatioNoiseFactor = 1.;
        final double minorAlleleFractionNoiseFactor = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.2, 0.7));
        final TumorHeterogeneityState.VariantProfile variantProfile1 = new TumorHeterogeneityState.VariantProfile(
                new TumorHeterogeneityState.VariantProfile.PloidyStateIndicators(Arrays.asList(0, 1)));
        final TumorHeterogeneityState.VariantProfile variantProfile2 = new TumorHeterogeneityState.VariantProfile(
                new TumorHeterogeneityState.VariantProfile.PloidyStateIndicators(Arrays.asList(2, 0)));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Arrays.asList(variantProfile1, variantProfile2));
        STATE = new TumorHeterogeneityState(concentration, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationFractions, variantProfiles, PRIORS);

        //need valid segment-mean posterior summary to construct TumorHeterogeneityData, but it is not used in tests
        final PosteriorSummary segmentMeanPosteriorSummary = new PosteriorSummary(0., -0.1, 0.1);
        segmentMeanPosteriorSummary.setDeciles(new DecileCollection(Arrays.asList(0., -0.1, 0.1), DecileCollection.ConstructionMode.SAMPLES));
        final ACNVModeledSegment segment1 = new ACNVModeledSegment(new SimpleInterval("1", 1, 25), segmentMeanPosteriorSummary, DUMMY_POSTERIOR_SUMMARY);
        final ACNVModeledSegment segment2 = new ACNVModeledSegment(new SimpleInterval("1", 26, 100), segmentMeanPosteriorSummary, DUMMY_POSTERIOR_SUMMARY);
        DATA = new TumorHeterogeneityData(Arrays.asList(segment1, segment2));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSinglePopulation() {
        //fail if only one population (must have at least one variant and one normal)
        final double concentration = 1.;
        final double copyRatioNoiseFactor = 1.;
        final double minorAlleleFractionNoiseFactor = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Collections.singletonList(1.));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.emptyList());
        new TumorHeterogeneityState(concentration, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationFractions, variantProfiles, PRIORS);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testUnnormalizedPopulationFractions() {
        //fail if population fractions not normalized to unity
        final double concentration = 1.;
        final double copyRatioNoiseFactor = 1.;
        final double minorAlleleFractionNoiseFactor = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.1));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.emptyList());
        new TumorHeterogeneityState(concentration, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationFractions, variantProfiles, PRIORS);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNegativePopulationIndicators() {
        //fail if population indicators are negative
        final double concentration = 1.;
        final double copyRatioNoiseFactor = 1.;
        final double minorAlleleFractionNoiseFactor = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.9));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.emptyList());
        new TumorHeterogeneityState(concentration, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationFractions, variantProfiles, PRIORS);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInconsistentPopulationIndicators() {
        //fail if population indicators are inconsistent with total number of populations
        final double concentration = 1.;
        final double copyRatioNoiseFactor = 1.;
        final double minorAlleleFractionNoiseFactor = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.9));
        final TumorHeterogeneityState.VariantProfile variantProfile = new TumorHeterogeneityState.VariantProfile(
                new TumorHeterogeneityState.VariantProfile.PloidyStateIndicators(Collections.singletonList(0)));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.singletonList(variantProfile));
        new TumorHeterogeneityState(concentration, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationFractions, variantProfiles, PRIORS);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNoVariants() {
        //fail if number of variants is not positive
        final double concentration = 1.;
        final double copyRatioNoiseFactor = 1.;
        final double minorAlleleFractionNoiseFactor = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.9));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.emptyList());
        new TumorHeterogeneityState(concentration, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationFractions, variantProfiles, PRIORS);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInconsistentNumbersOfPopulationsAndVariants() {
        //fail if number of variants + 1 is not equal to number of populations
        final double concentration = 1.;
        final double copyRatioNoiseFactor = 1.;
        final double minorAlleleFractionNoiseFactor = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.2, 0.7));
        final TumorHeterogeneityState.VariantProfile variantProfile = new TumorHeterogeneityState.VariantProfile(
                new TumorHeterogeneityState.VariantProfile.PloidyStateIndicators(Collections.singletonList(0)));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.singletonList(variantProfile));
        new TumorHeterogeneityState(concentration, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationFractions, variantProfiles, PRIORS);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDifferentNumberOfSegmentsAcrossVariants() {
        //fail if number of segments is not the same for all variants
        final double concentration = 1.;
        final double copyRatioNoiseFactor = 1.;
        final double minorAlleleFractionNoiseFactor = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.2, 0.7));
        final TumorHeterogeneityState.VariantProfile variantProfile1 = new TumorHeterogeneityState.VariantProfile(
                new TumorHeterogeneityState.VariantProfile.PloidyStateIndicators(Collections.singletonList(0)));
        final TumorHeterogeneityState.VariantProfile variantProfile2 = new TumorHeterogeneityState.VariantProfile(
                new TumorHeterogeneityState.VariantProfile.PloidyStateIndicators(Arrays.asList(0, 1)));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Arrays.asList(variantProfile1, variantProfile2));
        new TumorHeterogeneityState(concentration, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationFractions, variantProfiles, PRIORS);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDifferentNumberOfSegmentsWithinVariant() {
        //fail if number of segments is not the same for variant and ploidy-state indicators
        final double concentration = 1.;
        final double copyRatioNoiseFactor = 1.;
        final double minorAlleleFractionNoiseFactor = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.9));
        final TumorHeterogeneityState.VariantProfile variantProfile = new TumorHeterogeneityState.VariantProfile(
                new TumorHeterogeneityState.VariantProfile.PloidyStateIndicators(Arrays.asList(0, 1)));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.singletonList(variantProfile));
        new TumorHeterogeneityState(concentration, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationFractions, variantProfiles, PRIORS);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInconsistentPloidyStateIndicators() {
        //fail if ploidy-state indicators are inconsistent with total number of states in prior
        final double concentration = 1.;
        final double copyRatioNoiseFactor = 1.;
        final double minorAlleleFractionNoiseFactor = 1.;
        final TumorHeterogeneityState.PopulationFractions populationFractions = new TumorHeterogeneityState.PopulationFractions(Arrays.asList(0.1, 0.9));
        final TumorHeterogeneityState.VariantProfile variantProfile = new TumorHeterogeneityState.VariantProfile(
                new TumorHeterogeneityState.VariantProfile.PloidyStateIndicators(Collections.singletonList(3)));
        final TumorHeterogeneityState.VariantProfileCollection variantProfiles = new TumorHeterogeneityState.VariantProfileCollection(Collections.singletonList(variantProfile));
        new TumorHeterogeneityState(concentration, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, populationFractions, variantProfiles, PRIORS);
    }

    @Test
    public void testCalculateFractionalSegmentLength() {
        Assert.assertEquals(DATA.fractionalLength(0), 0.25);
        Assert.assertEquals(DATA.fractionalLength(1), 0.75);
    }

    @Test
    public void testCalculatePopulationAndGenomicAveragedPloidy() {
        Assert.assertEquals(STATE.ploidy(DATA), 1.925, EPSILON);
    }
}