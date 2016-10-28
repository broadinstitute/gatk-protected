package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityPriorCollection {
    private final PloidyState normalPloidyState;
    private final int normalPloidyStateIndex;
    private final PloidyStatePrior ploidyStatePrior;
    private final double concentrationPriorAlpha;
    private final double concentrationPriorBeta;

    public TumorHeterogeneityPriorCollection(final PloidyState normalPloidyState,
                                             final PloidyStatePrior ploidyStatePrior,
                                             final double concentrationPriorAlpha,
                                             final double concentrationPriorBeta) {
        Utils.nonNull(normalPloidyState);
        Utils.nonNull(ploidyStatePrior);
        Utils.validateArg(ploidyStatePrior.ploidyStates().contains(normalPloidyState),
                "Ploidy-state prior must contain normal ploidy state.");
        Utils.validateArg(concentrationPriorAlpha > 0, "Hyperparameter for concentration prior must be positive.");
        Utils.validateArg(concentrationPriorBeta > 0, "Hyperparameter for concentration prior must be positive.");
        this.normalPloidyState = normalPloidyState;
        this.ploidyStatePrior = ploidyStatePrior;
        normalPloidyStateIndex = ploidyStatePrior().ploidyStates().indexOf(normalPloidyState);
        this.concentrationPriorAlpha = concentrationPriorAlpha;
        this.concentrationPriorBeta = concentrationPriorBeta;
    }

    public PloidyState normalPloidyState() {
        return normalPloidyState;
    }

    public int normalPloidyStateIndex() {
        return normalPloidyStateIndex;
    }

    public PloidyStatePrior ploidyStatePrior() {
        return ploidyStatePrior;
    }

    public double concentrationPriorAlpha() {
        return concentrationPriorAlpha;
    }

    public double concentrationPriorBeta() {
        return concentrationPriorBeta;
    }
}
