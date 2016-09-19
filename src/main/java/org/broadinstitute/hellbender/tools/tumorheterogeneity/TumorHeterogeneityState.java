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
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityState extends ParameterizedState<TumorHeterogeneityParameter> {
    public static final class PopulationFractions extends ArrayList<Double> {
        //list of doubles, size = number of populations (including normal), i-th element = fraction of population i by cell number
        //normal population is last element
        private static final long serialVersionUID = 79454L;
        private final int numPopulations;
        public PopulationFractions(final List<Double> populationFractions) {
            super(populationFractions);
            Utils.validateArg(populationFractions.size() > 1, "Number of populations must be strictly greater than 1.");
            final double populationFractionNormalization = populationFractions.stream().mapToDouble(Double::doubleValue).sum();
            Utils.validateArg(Math.abs(1. - populationFractionNormalization) <= POPULATION_FRACTION_NORMALIZATION_EPSILON,
                    "Population fractions must sum to unity.");
            numPopulations = populationFractions.size();
        }
    }

    public static final class PopulationIndicators extends ArrayList<Integer> {
        //list of integers, size = number of cells, i-th element = population index for cell i
        private static final long serialVersionUID = 81915L;
        private final int numCells;
        public PopulationIndicators(final List<Integer> populationIndicators) {
            super(populationIndicators);
            Utils.validateArg(populationIndicators.size() > 0, "Number of cells must be positive.");
            Utils.validateArg(populationIndicators.stream().allMatch(i -> i >= 0), "Population indicators must be non-negative.");
            numCells = populationIndicators.size();
        }
    }

    public static final class VariantProfileCollection extends ArrayList<VariantProfile> {
        //list of VariantProfiles, size = number of populations (excluding normal), i-th element = VariantProfile for population i
        private static final long serialVersionUID = 76498L;
        private final int numSegments;
        private final int numVariantPopulations;
        public VariantProfileCollection(final List<VariantProfile> variantProfiles) {
            super(variantProfiles);
            Utils.validateArg(variantProfiles.size() > 0, "Number of variants must be positive.");
            final int numSegmentsForFirstVariant = variantProfiles.get(0).numSegments;
            Utils.validateArg(variantProfiles.stream().map(s -> s.numSegments).allMatch(n -> n == numSegmentsForFirstVariant),
                    "Number of segments must be same for all variants.");
            Utils.validateArg(numSegmentsForFirstVariant > 0, "Number of segments must be positive.");
            numSegments = numSegmentsForFirstVariant;
            numVariantPopulations = variantProfiles.size();
        }
    }

    private static final double POPULATION_FRACTION_NORMALIZATION_EPSILON = 1E-3;

    private final int numPopulations;   //variant populations + normal population
    private final int numCells;
    private final int numSegments;
    private final List<Integer> populationCounts;

    private final TumorHeterogeneityPriorCollection priors;

    public TumorHeterogeneityState(final double concentration,
                                   final PopulationFractions populationFractions,
                                   final PopulationIndicators populationIndicators,
                                   final VariantProfileCollection variantProfileCollection,
                                   final TumorHeterogeneityPriorCollection priors) {
        super(Arrays.asList(
                new Parameter<>(TumorHeterogeneityParameter.CONCENTRATION, concentration),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_FRACTIONS, populationFractions),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_INDICATORS, populationIndicators),
                new Parameter<>(TumorHeterogeneityParameter.VARIANT_PROFILES, variantProfileCollection)));
        Utils.validateArg(populationFractions.numPopulations == variantProfileCollection.numVariantPopulations + 1,
                "Number of populations must be equal to number of variant populations + 1.");
        Utils.validateArg(Collections.max(populationIndicators) < populationFractions.numPopulations,
                "Population indicators must be consistent with number of populations.");
        Utils.nonNull(priors);
        numPopulations = populationFractions.numPopulations;
        numCells = populationIndicators.numCells;
        numSegments = variantProfileCollection.numSegments;
        populationCounts = Collections.unmodifiableList(sumPopulationCounts());
        this.priors = priors;
    }

    public int numPopulations() {
        return numPopulations;
    }

    public int numCells() {
        return numCells;
    }

    public int numSegments() {
        return numSegments;
    }

    public double concentration() {
        return get(TumorHeterogeneityParameter.CONCENTRATION, Double.class);
    }

    public double populationFraction(final int populationIndex) {
        return get(TumorHeterogeneityParameter.POPULATION_FRACTIONS, PopulationFractions.class).get(populationIndex);
    }

    public int populationCount(final int populationIndex) {
        return populationCounts.get(populationIndex);
    }

    public int populationIndex(final int cellIndex) {
        return get(TumorHeterogeneityParameter.POPULATION_INDICATORS, PopulationIndicators.class).get(cellIndex);
    }

    public double variantSegmentFraction(final int populationIndex) {
        if (populationIndex == numPopulations - 1) {
            throw new IllegalStateException("Attempted to get variant-segment fraction for normal population.");
        }
        return get(TumorHeterogeneityParameter.VARIANT_PROFILES, VariantProfileCollection.class).get(populationIndex).variantSegmentFraction;
    }

    public boolean isVariant(final int populationIndex, final int segmentIndex) {
        if (populationIndex == numPopulations - 1) {
            //last population is normal
            return false;
        }
        return get(TumorHeterogeneityParameter.VARIANT_PROFILES, VariantProfileCollection.class).get(populationIndex).variantIndicators.get(segmentIndex);
    }

    public int variantPloidyStateIndex(final int populationIndex, final int segmentIndex) {
        if (populationIndex == numPopulations - 1) {
            throw new IllegalStateException("Attempted to get variant-ploidy-state index for normal population.");
        }
        return get(TumorHeterogeneityParameter.VARIANT_PROFILES, VariantProfileCollection.class).get(populationIndex).variantPloidyStateIndicators.get(segmentIndex);
    }

    public TumorHeterogeneityPriorCollection priors() {
        return priors;
    }

    public double calculatePopulationWeightedCopyNumber(final TumorHeterogeneityData data, final int segmentIndex) {
        Utils.nonNull(data);
        Utils.validateArg(data.numSegments() == numSegments,
                "Tumor-heterogeneity state and data collection must have same number of segments.");

        double numCopiesInSegment = 0.;
        final List<PloidyState> variantPloidyStates = priors.variantPloidyStatePrior().ploidyStates();
        final int normalPloidy = priors.normalPloidyState().total();
        for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
            numCopiesInSegment += isVariant(populationIndex, segmentIndex) ?
                    populationCount(populationIndex) * variantPloidyStates.get(variantPloidyStateIndex(populationIndex, segmentIndex)).total() :
                    populationCount(populationIndex) * normalPloidy;
        }
        return numCopiesInSegment / numCells;
    }

    public double calculateMinorAlleleFraction(final TumorHeterogeneityData data, final int segmentIndex) {
        Utils.nonNull(data);
        Utils.validateArg(data.numSegments() == numSegments,
                "Tumor-heterogeneity state and data collection must have same number of segments.");

        double numCopiesOfFirstAlleleInSegment = 0.;
        double numCopiesOfSecondAlleleInSegment = 0.;
        final List<PloidyState> variantPloidyStates = priors.variantPloidyStatePrior().ploidyStates();
        for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
            numCopiesOfFirstAlleleInSegment += isVariant(populationIndex, segmentIndex) ?
                    populationCount(populationIndex) * variantPloidyStates.get(variantPloidyStateIndex(populationIndex, segmentIndex)).m() :
                    populationCount(populationIndex) * priors.normalPloidyState().m();
            numCopiesOfSecondAlleleInSegment += isVariant(populationIndex, segmentIndex) ?
                    populationCount(populationIndex) * variantPloidyStates.get(variantPloidyStateIndex(populationIndex, segmentIndex)).n() :
                    populationCount(populationIndex) * priors.normalPloidyState().n();
        }
        return Math.min(numCopiesOfFirstAlleleInSegment, numCopiesOfSecondAlleleInSegment) / (numCopiesOfFirstAlleleInSegment + numCopiesOfSecondAlleleInSegment);
    }

    public double calculatePopulationWeightedGenomicAveragedPloidy(final TumorHeterogeneityData data) {
        Utils.nonNull(data);
        Utils.validateArg(data.numSegments() == numSegments,
                "Tumor-heterogeneity state and data collection must have same number of segments.");

        double numCopiesWeightedByGenomicLength = 0.;
        final List<PloidyState> variantPloidyStates = priors.variantPloidyStatePrior().ploidyStates();
        final int normalPloidy = priors.normalPloidyState().total();
        for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
            double numCopiesInSegment = 0.;
            for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
                numCopiesInSegment += isVariant(populationIndex, segmentIndex) ?
                        populationCount(populationIndex) * variantPloidyStates.get(variantPloidyStateIndex(populationIndex, segmentIndex)).total() :
                        populationCount(populationIndex) * normalPloidy;
            }
            numCopiesWeightedByGenomicLength += numCopiesInSegment * data.segmentLength(segmentIndex);
        }
        return numCopiesWeightedByGenomicLength / (numCells * data.totalLength());
    }

    /**
     * For each variant population, represents variant-segment fraction and per-segment variant and ploidy indicators.
     */
    static class VariantProfile {
        public static final class VariantIndicators extends ArrayList<Boolean> {
            //list of booleans, size = number of segments, i-th element = true if segment i is variant
            private static final long serialVersionUID = 35746L;
            private final int numSegments;
            public VariantIndicators(final List<Boolean> variantIndicators) {
                super(variantIndicators);
                Utils.validateArg(variantIndicators.size() > 0, "Number of segments must be positive.");
                numSegments = variantIndicators.size();
            }
        }

        public static final class VariantPloidyStateIndicators extends ArrayList<Integer> {
            //list of integers, size = number of segments, i-th element = variant-ploidy-state index of segment i
            private static final long serialVersionUID = 78476L;
            private final int numSegments;
            public VariantPloidyStateIndicators(final List<Integer> variantPloidyStateIndicators) {
                super(variantPloidyStateIndicators);
                Utils.validateArg(variantPloidyStateIndicators.size() > 0, "Number of segments must be positive.");
                numSegments = variantPloidyStateIndicators.size();
            }
        }

        private final int numSegments;
        private final double variantSegmentFraction;
        private final VariantIndicators variantIndicators;
        private final VariantPloidyStateIndicators variantPloidyStateIndicators;

        VariantProfile(final double variantSegmentFraction,
                       final VariantIndicators variantIndicators,
                       final VariantPloidyStateIndicators variantPloidyStateIndicators) {
            Utils.validateArg(0. <= variantSegmentFraction && variantSegmentFraction <= 1.,
                    "Variant-segment fraction must be in [0, 1].");
            Utils.validateArg(variantIndicators.numSegments == variantPloidyStateIndicators.numSegments,
                    "Number of segments must be same for variant and ploidy-state indicators.");
            numSegments = variantIndicators.numSegments;
            this.variantSegmentFraction = variantSegmentFraction;
            this.variantIndicators = variantIndicators;
            this.variantPloidyStateIndicators = variantPloidyStateIndicators;
        }
    }

    private List<Integer> sumPopulationCounts() {
        final List<MutableInt> populationCounts = IntStream.range(0, numPopulations).boxed()
                .map(j -> new MutableInt(0)).collect(Collectors.toList());
        for (int cellIndex = 0; cellIndex < numCells; cellIndex++) {
            final int populationIndex = populationIndex(cellIndex);
            populationCounts.get(populationIndex).increment();
        }
        return populationCounts.stream().map(MutableInt::intValue).collect(Collectors.toList());
    }
}
