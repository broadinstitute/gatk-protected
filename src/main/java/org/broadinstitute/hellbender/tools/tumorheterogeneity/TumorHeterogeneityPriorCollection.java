package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityPriorCollection {
    private final PloidyState normalPloidyState;
    private final PloidyStatePrior ploidyStatePrior;

    private final HyperparameterValues concentrationPriorHyperparameterValues;
    private final HyperparameterValues copyRatioNoiseConstantPriorHyperparameterValues;
    private final HyperparameterValues copyRatioNoiseFactorPriorHyperparameterValues;
    private final HyperparameterValues minorAlleleFractionNoiseFactorPriorHyperparameterValues;

    private final double ploidyMismatchPenalty;

    public TumorHeterogeneityPriorCollection(final PloidyState normalPloidyState,
                                             final PloidyStatePrior ploidyStatePrior,
                                             final double concentrationPriorAlpha,
                                             final double concentrationPriorBeta,
                                             final double copyRatioNoiseConstantPriorAlpha,
                                             final double copyRatioNoiseConstantPriorBeta,
                                             final double copyRatioNoiseFactorPriorAlpha,
                                             final double copyRatioNoiseFactorPriorBeta,
                                             final double minorAlleleFractionNoiseFactorPriorAlpha,
                                             final double minorAlleleFractionNoiseFactorPriorBeta,
                                             final double ploidyMismatchPenalty) {
        Utils.nonNull(normalPloidyState);
        Utils.nonNull(ploidyStatePrior);
        Utils.validateArg(ploidyStatePrior.ploidyStates().contains(normalPloidyState),
                "Ploidy-state prior must contain normal ploidy state.");
        Utils.validateArg(ploidyMismatchPenalty >= 0, "Ploidy-mismatch penalty must be non-negative.");
        this.normalPloidyState = normalPloidyState;
        this.ploidyStatePrior = ploidyStatePrior;
        concentrationPriorHyperparameterValues = new HyperparameterValues(concentrationPriorAlpha, concentrationPriorBeta);
        copyRatioNoiseConstantPriorHyperparameterValues = new HyperparameterValues(copyRatioNoiseConstantPriorAlpha, copyRatioNoiseConstantPriorBeta);
        copyRatioNoiseFactorPriorHyperparameterValues = new HyperparameterValues(copyRatioNoiseFactorPriorAlpha, copyRatioNoiseFactorPriorBeta);
        minorAlleleFractionNoiseFactorPriorHyperparameterValues = new HyperparameterValues(minorAlleleFractionNoiseFactorPriorAlpha, minorAlleleFractionNoiseFactorPriorBeta);
        this.ploidyMismatchPenalty = ploidyMismatchPenalty;
    }

    public PloidyState normalPloidyState() {
        return normalPloidyState;
    }

    public PloidyStatePrior ploidyStatePrior() {
        return ploidyStatePrior;
    }

    public HyperparameterValues concentrationPriorHyperparameterValues() {
        return concentrationPriorHyperparameterValues;
    }

    public HyperparameterValues copyRatioNoiseConstantPriorHyperparameterValues() {
        return copyRatioNoiseConstantPriorHyperparameterValues;
    }

    public HyperparameterValues copyRatioNoiseFactorPriorHyperparameterValues() {
        return copyRatioNoiseFactorPriorHyperparameterValues;
    }

    public HyperparameterValues minorAlleleFractionNoiseFactorPriorHyperparameterValues() {
        return minorAlleleFractionNoiseFactorPriorHyperparameterValues;
    }

    public double ploidyMismatchPenalty() {
        return ploidyMismatchPenalty;
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
