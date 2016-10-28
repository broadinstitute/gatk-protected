package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityState extends ParameterizedState<TumorHeterogeneityParameter> {
    private static final double POPULATION_FRACTION_NORMALIZATION_EPSILON = 1E-4;

    private final int numPopulations;   //variant populations + normal population
    private final int numSegments;
    private final int numPloidyStates;

    private final TumorHeterogeneityPriorCollection priors;

    public TumorHeterogeneityState(final double concentration,
                                   final PopulationFractions populationFractions,
                                   final VariantProfileCollection variantProfileCollection,
                                   final TumorHeterogeneityPriorCollection priors) {
        super(Arrays.asList(
                new Parameter<>(TumorHeterogeneityParameter.CONCENTRATION, concentration),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_FRACTIONS, populationFractions),
                new Parameter<>(TumorHeterogeneityParameter.VARIANT_PROFILES, variantProfileCollection)));
        Utils.validateArg(concentration > 0, "Concentration must be positive.");
        Utils.nonNull(populationFractions);
        Utils.nonNull(variantProfileCollection);
        Utils.validateArg(populationFractions.numPopulations == variantProfileCollection.numVariantPopulations + 1,
                "Number of populations must be equal to number of variant populations + 1.");
        Utils.validateArg(variantProfileCollection.stream().map(s -> Collections.max(s.ploidyStateIndicators)).allMatch(i -> i < priors.ploidyStatePrior().numPloidyStates()),
                "Ploidy-state indicators must be consistent with number of ploidy states.");
        Utils.nonNull(priors);
        numPopulations = populationFractions.numPopulations;
        numSegments = variantProfileCollection.numSegments;
        numPloidyStates = priors.ploidyStatePrior().numPloidyStates();
        this.priors = priors;
    }

    /*===============================================================================================================*
     * GETTERS                                                                                                       *
     *===============================================================================================================*/

    public int numPopulations() {
        return numPopulations;
    }

    public int numSegments() {
        return numSegments;
    }

    public double concentration() {
        return get(TumorHeterogeneityParameter.CONCENTRATION, Double.class);
    }

    public double populationFraction(final int populationIndex) {
        validatePopulationIndex(populationIndex, numPopulations);
        return get(TumorHeterogeneityParameter.POPULATION_FRACTIONS, PopulationFractions.class).get(populationIndex);
    }

    TumorHeterogeneityState.PopulationFractions populationFractions() {
        return get(TumorHeterogeneityParameter.POPULATION_FRACTIONS, PopulationFractions.class);
    }

    public boolean isNormalPopulation(final int populationIndex) {
        validatePopulationIndex(populationIndex, numPopulations);
        return populationIndex == numPopulations - 1;
    }

    public int ploidyStateIndex(final int populationIndex, final int segmentIndex) {
        validatePopulationIndex(populationIndex, numPopulations);
        validateSegmentIndex(segmentIndex, numSegments);
        if (populationIndex == numPopulations - 1) {
            throw new IllegalStateException("Attempted to get ploidy-state index for normal population.");
        }
        return get(TumorHeterogeneityParameter.VARIANT_PROFILES, VariantProfileCollection.class).get(populationIndex).ploidyStateIndex(segmentIndex);
    }

    public PloidyState ploidyState(final int populationIndex, final int segmentIndex) {
        validatePopulationIndex(populationIndex, numPopulations);
        validateSegmentIndex(segmentIndex, numSegments);
        return isNormalPopulation(populationIndex) ?
                priors.normalPloidyState() :
                priors.ploidyStatePrior().ploidyStates().get(ploidyStateIndex(populationIndex, segmentIndex));
    }

    TumorHeterogeneityState.VariantProfileCollection variantProfiles() {
        return get(TumorHeterogeneityParameter.VARIANT_PROFILES, VariantProfileCollection.class);
    }

    public TumorHeterogeneityPriorCollection priors() {
        return priors;
    }

    /*===============================================================================================================*
     * METHODS                                                                                                       *
     *===============================================================================================================*/

    public int calculateCopyNumberFunction(final int segmentIndex,
                                           final Integer populationIndex,
                                           final Function<PloidyState, Integer> copyNumberFunction) {
        validateSegmentIndex(segmentIndex, numSegments);
        validatePopulationIndex(populationIndex, numPopulations);
        return isNormalPopulation(populationIndex) ?
                copyNumberFunction.apply(priors.normalPloidyState()) :
                copyNumberFunction.apply(ploidyState(populationIndex, segmentIndex));
    }

    public double calculatePopulationAveragedCopyNumberFunction(final int segmentIndex,
                                                                final Function<PloidyState, Integer> copyNumberFunction) {
        return calculatePopulationAveragedCopyNumberFunctionExcludingPopulation(null, segmentIndex, copyNumberFunction);
    }

    public double ploidy(final TumorHeterogeneityData data) {
        validateData(data, numSegments);
        return IntStream.range(0, numSegments).mapToDouble(i -> calculatePopulationAndGenomicAveragedCopyNumberFunction(null, i, PloidyState::total, data)).sum();
    }

    public double populationPloidy(final int populationIndex, final TumorHeterogeneityData data) {
        validatePopulationIndex(populationIndex, numPopulations);
        validateData(data, numSegments);
        return isNormalPopulation(populationIndex) ?
                priors.normalPloidyState().total() :
                variantProfiles().get(populationIndex).ploidy(priors.ploidyStatePrior().ploidyStates(), data);
    }

    double calculatePopulationAveragedCopyNumberFunctionExcludingPopulation(final Integer populationIndexToExclude,
                                                                            final int segmentIndex,
                                                                            final Function<PloidyState, Integer> copyNumberFunction) {
        validateSegmentIndex(segmentIndex, numSegments);
        final int populationIndexToExcludeValue;
        if (populationIndexToExclude == null) {
            populationIndexToExcludeValue = -1; //this will not exclude any populations
        } else {
            validatePopulationIndex(populationIndexToExclude, numPopulations);
            populationIndexToExcludeValue = populationIndexToExclude;
        }

        return IntStream.range(0, numPopulations)
                .filter(populationIndex -> populationIndex != populationIndexToExcludeValue)
                .mapToDouble(populationIndex -> populationFraction(populationIndex) * calculateCopyNumberFunction(segmentIndex, populationIndex, copyNumberFunction))
                .sum();
    }

    private double calculatePopulationAndGenomicAveragedCopyNumberFunction(final Integer populationIndexToExclude,
                                                                           final int segmentIndex,
                                                                           final Function<PloidyState, Integer> copyNumberFunction,
                                                                           final TumorHeterogeneityData data) {
        validateData(data, numSegments);
        validateSegmentIndex(segmentIndex, numSegments);

        return data.fractionalLength(segmentIndex) * calculatePopulationAveragedCopyNumberFunctionExcludingPopulation(populationIndexToExclude, segmentIndex, copyNumberFunction);
    }

    /*===============================================================================================================*
     * SETTERS (SHOULD ONLY BE USED BY SAMPLERS TO MODIFY STATE FOR SAMPLING OF INDICATOR VARIABLES)                 *
     *===============================================================================================================*/

    void setPloidyStateIndex(final int populationIndex, final int segmentIndex, final int ploidyStateIndex) {
        validatePopulationIndex(populationIndex, numPopulations);
        validateSegmentIndex(segmentIndex, numSegments);
        validatePloidyStateIndex(ploidyStateIndex, numPloidyStates);
        if (populationIndex == numPopulations - 1) {
            throw new IllegalStateException("Attempted to set ploidy-state index for normal population.");
        }
        get(TumorHeterogeneityParameter.VARIANT_PROFILES, VariantProfileCollection.class).get(populationIndex).ploidyStateIndicators.set(segmentIndex, ploidyStateIndex);
    }

    void set(final TumorHeterogeneityState state) {
        update(TumorHeterogeneityParameter.CONCENTRATION, state.concentration());
        update(TumorHeterogeneityParameter.POPULATION_FRACTIONS, state.populationFractions());
        update(TumorHeterogeneityParameter.VARIANT_PROFILES, state.variantProfiles());
    }

    /*===============================================================================================================*
     * INNER CLASSES                                                                                                 *
     *===============================================================================================================*/

    public static final class PopulationFractions extends ArrayList<Double> {
        //list of doubles, size = number of populations (including normal), i-th element = fraction of population i by cell number
        //normal population is last element
        private static final long serialVersionUID = 79454L;
        private final int numPopulations;
        public PopulationFractions(final List<Double> populationFractions) {
            super(Utils.nonNull(populationFractions));
            Utils.validateArg(populationFractions.size() > 1, "Number of populations must be strictly greater than 1.");
            final double populationFractionNormalization = populationFractions.stream().mapToDouble(Double::doubleValue).sum();
            Utils.validateArg(Math.abs(1. - populationFractionNormalization) <= POPULATION_FRACTION_NORMALIZATION_EPSILON,
                    "Population fractions must sum to unity.");
            numPopulations = populationFractions.size();
        }
    }

    public static final class VariantProfileCollection extends ArrayList<VariantProfile> {
        //list of VariantProfiles, size = number of populations (excluding normal), i-th element = VariantProfile for population i
        private static final long serialVersionUID = 76498L;
        private final int numSegments;
        private final int numVariantPopulations;
        public VariantProfileCollection(final List<VariantProfile> variantProfiles) {
            super(Utils.nonNull(variantProfiles));
            Utils.validateArg(variantProfiles.size() > 0, "Number of variants must be positive.");
            final int numSegmentsForFirstVariant = variantProfiles.get(0).numSegments;
            Utils.validateArg(variantProfiles.stream().map(s -> s.numSegments).allMatch(n -> n == numSegmentsForFirstVariant),
                    "Number of segments must be same for all variants.");
            Utils.validateArg(numSegmentsForFirstVariant > 0, "Number of segments must be positive.");
            numSegments = numSegmentsForFirstVariant;
            numVariantPopulations = variantProfiles.size();
        }
    }

    /**
     * For each variant population, represents ploidy-state indicators.
     */
    static class VariantProfile {
        public static final class PloidyStateIndicators extends ArrayList<Integer> {
            //list of integers, size = number of segments, i-th element = ploidy-state index of segment i
            private static final long serialVersionUID = 78476L;
            private final int numSegments;
            public PloidyStateIndicators(final List<Integer> ploidyStateIndicators) {
                super(Utils.nonNull(ploidyStateIndicators));
                Utils.validateArg(ploidyStateIndicators.size() > 0, "Number of segments must be positive.");
                Utils.validateArg(ploidyStateIndicators.stream().allMatch(i -> i >= 0), "Population indicators must be non-negative.");
                numSegments = ploidyStateIndicators.size();
            }
        }

        private final int numSegments;
        private final PloidyStateIndicators ploidyStateIndicators;

        VariantProfile(final PloidyStateIndicators ploidyStateIndicators) {
            Utils.nonNull(ploidyStateIndicators);
            numSegments = ploidyStateIndicators.numSegments;
            this.ploidyStateIndicators = ploidyStateIndicators;
        }

        public int ploidyStateIndex(final int segmentIndex) {
            validateSegmentIndex(segmentIndex, numSegments);
            return ploidyStateIndicators.get(segmentIndex);
        }

        public int numSegments() {
            return numSegments;
        }

        public double ploidy(final List<PloidyState> ploidyStates,
                             final TumorHeterogeneityData data) {
            validateData(data, numSegments);
            return IntStream.range(0, numSegments).mapToDouble(i -> data.fractionalLength(i) * ploidyStates.get(ploidyStateIndex(i)).total()).sum();
        }
    }

    /*===============================================================================================================*
     * OTHER PRIVATE METHODS                                                                                         *
     *===============================================================================================================*/

    private static void validatePopulationIndex(final int populationIndex, final int numPopulations) {
        Utils.validateArg(0 <= populationIndex && populationIndex < numPopulations, "Population index out of range.");
    }

    private static void validateSegmentIndex(int segmentIndex, final int numSegments) {
        Utils.validateArg(0 <= segmentIndex && segmentIndex < numSegments, "Segment index out of range.");
    }

    private static void validatePloidyStateIndex(final int ploidyStateIndex, final int numPloidyStates) {
        Utils.validateArg(0 <= ploidyStateIndex && ploidyStateIndex < numPloidyStates, "Ploidy-state index out of range.");
    }

    private static void validateData(final TumorHeterogeneityData data, final int numSegments) {
        Utils.nonNull(data);
        Utils.validateArg(data.numSegments() == numSegments,
                "Tumor-heterogeneity state and data collection must have same number of segments.");
    }
}
