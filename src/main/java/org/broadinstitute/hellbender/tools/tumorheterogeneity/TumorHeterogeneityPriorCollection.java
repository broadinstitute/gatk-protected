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
    private final double proposalWidthFactor;
    private final double concentrationPriorAlpha;
    private final double concentrationPriorBeta;
    private final double copyRatioNoiseFloorPriorAlpha;
    private final double copyRatioNoiseFloorPriorBeta;
    private final double copyRatioNoiseFactorPriorAlpha;
    private final double copyRatioNoiseFactorPriorBeta;
    private final double minorAlleleFractionNoiseFactorPriorAlpha;
    private final double minorAlleleFractionNoiseFactorPriorBeta;

    public TumorHeterogeneityPriorCollection(final PloidyState normalPloidyState,
                                             final PloidyStatePrior ploidyStatePrior,
                                             final double proposalWidthFactor,
                                             final double concentrationPriorAlpha,
                                             final double concentrationPriorBeta,
                                             final double copyRatioNoiseFloorPriorAlpha,
                                             final double copyRatioNoiseFloorPriorBeta,
                                             final double copyRatioNoiseFactorPriorAlpha,
                                             final double copyRatioNoiseFactorPriorBeta,
                                             final double minorAlleleFractionNoiseFactorPriorAlpha,
                                             final double minorAlleleFractionNoiseFactorPriorBeta) {
        Utils.nonNull(normalPloidyState);
        Utils.nonNull(ploidyStatePrior);
        Utils.validateArg(ploidyStatePrior.ploidyStates().contains(normalPloidyState),
                "Ploidy-state prior must contain normal ploidy state.");
        Utils.validateArg(proposalWidthFactor > 0, "Proposal-width factor must be positive.");
        Utils.validateArg(concentrationPriorAlpha > 0, "Hyperparameter for concentration prior must be positive.");
        Utils.validateArg(concentrationPriorBeta > 0, "Hyperparameter for concentration prior must be positive.");
        Utils.validateArg(copyRatioNoiseFloorPriorAlpha > 0, "Hyperparameter for copy-ratio noise-floor prior must be positive.");
        Utils.validateArg(copyRatioNoiseFloorPriorBeta > 0, "Hyperparameter for copy-ratio noise-floor prior must be positive.");
        Utils.validateArg(copyRatioNoiseFactorPriorAlpha > 0, "Hyperparameter for copy-ratio noise-factor prior must be positive.");
        Utils.validateArg(copyRatioNoiseFactorPriorBeta > 0, "Hyperparameter for copy-ratio noise-factor prior must be positive.");
        Utils.validateArg(minorAlleleFractionNoiseFactorPriorAlpha > 0, "Hyperparameter for minor-allele-fraction noise-factor prior must be positive.");
        Utils.validateArg(minorAlleleFractionNoiseFactorPriorBeta > 0, "Hyperparameter for minor-allele-fraction noise-factor prior must be positive.");
        this.normalPloidyState = normalPloidyState;
        this.ploidyStatePrior = ploidyStatePrior;
        normalPloidyStateIndex = ploidyStatePrior().ploidyStates().indexOf(normalPloidyState);
        this.proposalWidthFactor = proposalWidthFactor;
        this.concentrationPriorAlpha = concentrationPriorAlpha;
        this.concentrationPriorBeta = concentrationPriorBeta;
        this.copyRatioNoiseFloorPriorAlpha = copyRatioNoiseFloorPriorAlpha;
        this.copyRatioNoiseFloorPriorBeta = copyRatioNoiseFloorPriorBeta;
        this.copyRatioNoiseFactorPriorAlpha = copyRatioNoiseFactorPriorAlpha;
        this.copyRatioNoiseFactorPriorBeta = copyRatioNoiseFactorPriorBeta;
        this.minorAlleleFractionNoiseFactorPriorAlpha = minorAlleleFractionNoiseFactorPriorAlpha;
        this.minorAlleleFractionNoiseFactorPriorBeta = minorAlleleFractionNoiseFactorPriorBeta;
    }

    public PloidyState normalPloidyState() {
        return normalPloidyState;
    }

    public int normalPloidyStateIndex() {
        return normalPloidyStateIndex;
    }

    public double proposalWidthFactor() {
        return proposalWidthFactor;
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

    public double copyRatioNoiseFloorPriorAlpha() {
        return copyRatioNoiseFloorPriorAlpha;
    }

    public double copyRatioNoiseFloorPriorBeta() {
        return copyRatioNoiseFloorPriorBeta;
    }

    public double copyRatioNoiseFactorPriorAlpha() {
        return copyRatioNoiseFactorPriorAlpha;
    }

    public double copyRatioNoiseFactorPriorBeta() {
        return copyRatioNoiseFactorPriorBeta;
    }

    public double minorAlleleFractionNoiseFactorPriorAlpha() {
        return minorAlleleFractionNoiseFactorPriorAlpha;
    }

    public double minorAlleleFractionNoiseFactorPriorBeta() {
        return minorAlleleFractionNoiseFactorPriorBeta;
    }
}
