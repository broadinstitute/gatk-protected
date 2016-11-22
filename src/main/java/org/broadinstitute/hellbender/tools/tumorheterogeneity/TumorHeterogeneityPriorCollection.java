package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityPriorCollection {
    private final PloidyState normalPloidyState;
    private final PloidyStatePrior ploidyStatePrior;
    private final int maxTotalCopyNumber;
    private final double proposalWidthFactor;

    private final HyperparameterValues concentrationPriorHyperparameterValues;
    private final HyperparameterValues copyRatioNoiseFloorPriorHyperparameterValues;
    private final HyperparameterValues copyRatioNoiseFactorPriorHyperparameterValues;
    private final HyperparameterValues minorAlleleFractionNoiseFactorPriorHyperparameterValues;

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
        this.normalPloidyState = normalPloidyState;
        this.ploidyStatePrior = ploidyStatePrior;
        maxTotalCopyNumber = Collections.max(ploidyStatePrior().ploidyStates().stream().map(PloidyState::n).collect(Collectors.toList()));
        this.proposalWidthFactor = proposalWidthFactor;
        concentrationPriorHyperparameterValues = new HyperparameterValues(concentrationPriorAlpha, concentrationPriorBeta);
        copyRatioNoiseFloorPriorHyperparameterValues = new HyperparameterValues(copyRatioNoiseFloorPriorAlpha, copyRatioNoiseFloorPriorBeta);
        copyRatioNoiseFactorPriorHyperparameterValues = new HyperparameterValues(copyRatioNoiseFactorPriorAlpha, copyRatioNoiseFactorPriorBeta);
        minorAlleleFractionNoiseFactorPriorHyperparameterValues = new HyperparameterValues(minorAlleleFractionNoiseFactorPriorAlpha, minorAlleleFractionNoiseFactorPriorBeta);
    }

    public PloidyState normalPloidyState() {
        return normalPloidyState;
    }

    public PloidyStatePrior ploidyStatePrior() {
        return ploidyStatePrior;
    }

    public int maxTotalCopyNumber() {
        return maxTotalCopyNumber;
    }

    public double proposalWidthFactor() {
        return proposalWidthFactor;
    }

    public HyperparameterValues concentrationPriorHyperparameterValues() {
        return concentrationPriorHyperparameterValues;
    }

    public HyperparameterValues copyRatioNoiseFloorPriorHyperparameterValues() {
        return copyRatioNoiseFloorPriorHyperparameterValues;
    }

    public HyperparameterValues copyRatioNoiseFactorPriorHyperparameterValues() {
        return copyRatioNoiseFactorPriorHyperparameterValues;
    }

    public HyperparameterValues minorAlleleFractionNoiseFactorPriorHyperparameterValues() {
        return minorAlleleFractionNoiseFactorPriorHyperparameterValues;
    }

    public static class HyperparameterValues {
        private final double alpha;
        private final double beta;

        public HyperparameterValues(final double alpha, final double beta) {
            Utils.validateArg(alpha > 0, "Hyperparameter must be positive.");
            Utils.validateArg(beta > 0, "Hyperparameter must be positive.");
            this.alpha = alpha;
            this.beta = beta;
        }

        public double getAlpha() {
            return alpha;
        }

        public double getBeta() {
            return beta;
        }
    }
}
