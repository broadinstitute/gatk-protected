package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.Arrays;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityState extends ParameterizedState<TumorHeterogeneityParameter> {
    private final TumorHeterogeneityPriorCollection priors;

    public TumorHeterogeneityState(final double concentration,
                                   final double copyRatioNoiseFloor,
                                   final double copyRatioNoiseFactor,
                                   final double minorAlleleFractionNoiseFactor,
                                   final PopulationMixture populationMixture,
                                   final TumorHeterogeneityPriorCollection priors) {
        super(Arrays.asList(
                new Parameter<>(TumorHeterogeneityParameter.CONCENTRATION, concentration),
                new Parameter<>(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FLOOR, copyRatioNoiseFloor),
                new Parameter<>(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FACTOR, copyRatioNoiseFactor),
                new Parameter<>(TumorHeterogeneityParameter.MINOR_ALLELE_FRACTION_NOISE_FACTOR, minorAlleleFractionNoiseFactor),
                new Parameter<>(TumorHeterogeneityParameter.POPULATION_MIXTURE,
                        new PopulationMixture(populationMixture.populationFractions(), populationMixture.variantProfileCollection(), priors.normalPloidyState()))));
        Utils.validateArg(concentration > 0, "Concentration must be positive.");
        Utils.validateArg(copyRatioNoiseFactor >= 0, "Copy-ratio noise floor must be non-negative.");
        Utils.validateArg(copyRatioNoiseFactor >= 0, "Copy-ratio noise factor must be non-negative.");
        Utils.validateArg(minorAlleleFractionNoiseFactor >= 1, "Minor-allele-fraction noise factor must be >= 1.");
        Utils.nonNull(populationMixture);
        Utils.nonNull(priors);
        this.priors = priors;
    }

    public double concentration() {
        return get(TumorHeterogeneityParameter.CONCENTRATION, Double.class);
    }

    public double copyRatioNoiseFloor() {
        return get(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FLOOR, Double.class);
    }

    public double copyRatioNoiseFactor() {
        return get(TumorHeterogeneityParameter.COPY_RATIO_NOISE_FACTOR, Double.class);
    }

    public double minorAlleleFractionNoiseFactor() {
        return get(TumorHeterogeneityParameter.MINOR_ALLELE_FRACTION_NOISE_FACTOR, Double.class);
    }

    public PopulationMixture populationMixture() {
        return get(TumorHeterogeneityParameter.POPULATION_MIXTURE, PopulationMixture.class);
    }

    public TumorHeterogeneityPriorCollection priors() {
        return priors;
    }
}
