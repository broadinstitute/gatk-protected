package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.Arrays;
import java.util.Collections;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityState extends ParameterizedState<TumorHeterogeneityParameter> {
    public TumorHeterogeneityState(final double concentration,
                                   final double copyRatioNoiseConstant,
                                   final double copyRatioNoiseFactor,
                                   final double minorAlleleFractionNoiseFactor,
                                   final double initialPloidy,
                                   final double ploidy,
                                   final PopulationMixture populationMixture) {
        super(Arrays.asList(
                new Parameter<>(TumorHeterogeneityParameter.CONCENTRATION, concentration),
                new Parameter<>(TumorHeterogeneityParameter.COPY_RATIO_NOISE_CONSTANT, copyRatioNoiseConstant),
                new Parameter<>(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FACTOR, copyRatioNoiseFactor),
                new Parameter<>(TumorHeterogeneityParameter.MINOR_ALLELE_FRACTION_NOISE_FACTOR, minorAlleleFractionNoiseFactor),
                new Parameter<>(TumorHeterogeneityParameter.INITIAL_PLOIDY, initialPloidy),
                new Parameter<>(TumorHeterogeneityParameter.PLOIDY, ploidy),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_MIXTURE,
                        new PopulationMixture(populationMixture.populationFractions(), populationMixture.variantProfileCollection(), populationMixture.normalPloidyState()))));
        Utils.validateArg(concentration > 0., "Concentration must be positive.");
        Utils.validateArg(copyRatioNoiseConstant >= 0., "Copy-ratio noise constant must be non-negative.");
        Utils.validateArg(copyRatioNoiseFactor > 0., "Copy-ratio noise factor must be positive.");
        Utils.validateArg(minorAlleleFractionNoiseFactor > 0., "Minor-allele-fraction noise factor must be positive.");
        Utils.validateArg(initialPloidy >= 0, "Initial ploidy must be non-negative.");
        Utils.validateArg(ploidy >= 0, "Ploidy must be non-negative.");
        Utils.nonNull(populationMixture);
    }

    public double concentration() {
        return get(TumorHeterogeneityParameter.CONCENTRATION, Double.class);
    }

    public double copyRatioNoiseConstant() {
        return get(TumorHeterogeneityParameter.COPY_RATIO_NOISE_CONSTANT, Double.class);
    }

    public double copyRatioNoiseFactor() {
        return get(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FACTOR, Double.class);
    }

    public double minorAlleleFractionNoiseFactor() {
        return get(TumorHeterogeneityParameter.MINOR_ALLELE_FRACTION_NOISE_FACTOR, Double.class);
    }

    public double initialPloidy() {
        return get(TumorHeterogeneityParameter.INITIAL_PLOIDY, Double.class);
    }

    public double ploidy() {
        return get(TumorHeterogeneityParameter.PLOIDY, Double.class);
    }

    public PopulationMixture populationMixture() {
        return get(TumorHeterogeneityParameter.POPULATION_MIXTURE, PopulationMixture.class);
    }

    /**
     * Initialize state with evenly distributed population fractions and normal variant profiles.
     */
    static TumorHeterogeneityState initializeState(final TumorHeterogeneityPriorCollection priors,
                                                   final int numSegments,
                                                   final int numPopulations) {
        final double concentration = priors.concentrationPriorHyperparameterValues().getAlpha() / priors.concentrationPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseConstant = priors.copyRatioNoiseConstantPriorHyperparameterValues().getAlpha() / priors.copyRatioNoiseConstantPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseFactor = priors.copyRatioNoiseFactorPriorHyperparameterValues().getAlpha() /
                (priors.copyRatioNoiseFactorPriorHyperparameterValues().getAlpha() + priors.copyRatioNoiseFactorPriorHyperparameterValues().getBeta());
        final double minorAlleleFractionNoiseFactor = priors.minorAlleleFractionNoiseFactorPriorHyperparameterValues().getAlpha() /
                (priors.minorAlleleFractionNoiseFactorPriorHyperparameterValues().getAlpha() + priors.minorAlleleFractionNoiseFactorPriorHyperparameterValues().getBeta());
        //initialize population fractions to be evenly distributed
        final PopulationMixture.PopulationFractions populationFractions =
                new PopulationMixture.PopulationFractions(Collections.nCopies(numPopulations, 1. / numPopulations));
        //initialize variant profiles to normal
        final PloidyState normalPloidyState = priors.normalPloidyState();
        final double normalPloidy = normalPloidyState.total();
        final int numVariantPopulations = numPopulations - 1;
        final PopulationMixture.VariantProfileCollection variantProfileCollection =
                initializeNormalProfiles(numVariantPopulations, numSegments, normalPloidyState);
        final PopulationMixture populationMixture = new PopulationMixture(populationFractions, variantProfileCollection, normalPloidyState);
        return new TumorHeterogeneityState(
                concentration, copyRatioNoiseConstant, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, normalPloidy, normalPloidy, populationMixture);
    }

    /**
     * Initialize variant profiles to normal.
     */
    private static PopulationMixture.VariantProfileCollection initializeNormalProfiles(final int numVariantPopulations,
                                                                                       final int numSegments,
                                                                                       final PloidyState normalPloidyState) {
        return new PopulationMixture.VariantProfileCollection(
                Collections.nCopies(numVariantPopulations, initializeNormalProfile(numSegments, normalPloidyState)));
    }

    /**
     * Initialize a variant profile to normal.
     */
    private static PopulationMixture.VariantProfile initializeNormalProfile(final int numSegments,
                                                                            final PloidyState normalPloidyState) {
        final PopulationMixture.VariantProfile ploidyStateIndicators =
                new PopulationMixture.VariantProfile(Collections.nCopies(numSegments, normalPloidyState));
        return new PopulationMixture.VariantProfile(ploidyStateIndicators);
    }
}
