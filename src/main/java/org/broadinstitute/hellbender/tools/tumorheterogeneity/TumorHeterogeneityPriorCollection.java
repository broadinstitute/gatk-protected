package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityPriorCollection {
    private final PloidyState normalPloidyState;
    private final PloidyStatePrior variantPloidyStatePrior;
    private final double concentrationPriorAlpha;
    private final double concentrationPriorBeta;
    private final double variantSegmentFractionPriorAlpha;
    private final double variantSegmentFractionPriorBeta;

    public TumorHeterogeneityPriorCollection(final PloidyState normalPloidyState,
                                             final PloidyStatePrior variantPloidyStatePrior,
                                             final double concentrationPriorAlpha,
                                             final double concentrationPriorBeta,
                                             final double variantSegmentFractionPriorAlpha,
                                             final double variantSegmentFractionPriorBeta) {
        Utils.nonNull(normalPloidyState);
        Utils.nonNull(variantPloidyStatePrior);
        Utils.validateArg(Double.isNaN(variantPloidyStatePrior.logProbability(normalPloidyState)),
                "Variant-ploidy state prior should not be specified for normal ploidy state.");
        Utils.validateArg(concentrationPriorAlpha > 0, "Hyperparameter for concentration prior must be positive.");
        Utils.validateArg(concentrationPriorBeta > 0, "Hyperparameter for concentration prior must be positive.");
        Utils.validateArg(variantSegmentFractionPriorAlpha > 0, "Hyperparameter for variant-segment fraction must be positive.");
        Utils.validateArg(variantSegmentFractionPriorBeta > 0, "Hyperparameter for variant-segment fraction must be positive.");
        this.normalPloidyState = normalPloidyState;
        this.variantPloidyStatePrior = variantPloidyStatePrior;
        this.concentrationPriorAlpha = concentrationPriorAlpha;
        this.concentrationPriorBeta = concentrationPriorBeta;
        this.variantSegmentFractionPriorAlpha = variantSegmentFractionPriorAlpha;
        this.variantSegmentFractionPriorBeta = variantSegmentFractionPriorBeta;
    }

    public PloidyState normalPloidyState() {
        return normalPloidyState;
    }

    public PloidyStatePrior variantPloidyStatePrior() {
        return variantPloidyStatePrior;
    }

    public double concentrationPriorAlpha() {
        return concentrationPriorAlpha;
    }

    public double concentrationPriorBeta() {
        return concentrationPriorBeta;
    }

    public double variantSegmentFractionPriorAlpha() {
        return variantSegmentFractionPriorAlpha;
    }

    public double variantSegmentFractionPriorBeta() {
        return variantSegmentFractionPriorBeta;
    }
}
