package org.broadinstitute.hellbender.tools.exome.allelefractionrevision;

import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public final class AlleleFractionState extends ParameterizedState<AlleleFractionParameter> {
    public static final double MIN_MINOR_FRACTION = 0.0;   //by definition!
    public static final double MAX_MINOR_FRACTION = 0.5;   //by definition!

    public static final class MinorFractions extends ArrayList<Double> {
        private static final long serialVersionUID = 1029384756L;
        public MinorFractions(final int numSegments) { super(numSegments); }
        public MinorFractions(final List<Double> other) {
            super(new ArrayList<>(other));
        }
    }

    public AlleleFractionState(final double sampleDepth,
                               final double sampleBias,
                               final double outlierProbability,
                               final MinorFractions minorFractions) {
        super(Arrays.asList(
                new Parameter<>(AlleleFractionParameter.SAMPLE_DEPTH, sampleDepth),
                new Parameter<>(AlleleFractionParameter.SAMPLE_BIAS, sampleBias),
                new Parameter<>(AlleleFractionParameter.OUTLIER_PROBABILITY, outlierProbability),
                new Parameter<>(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, minorFractions)));
    }


    public double sampleDepth() {
        return get(AlleleFractionParameter.SAMPLE_DEPTH, Double.class);
    }

    public double sampleBias() {
        return get(AlleleFractionParameter.SAMPLE_BIAS, Double.class);
    }

    public double outlierProbability() {
        return get(AlleleFractionParameter.OUTLIER_PROBABILITY, Double.class);
    }

    public double segmentMinorFraction(final int segment) {
        return get(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, MinorFractions.class).get(segment);
    }

    public AlleleFractionGlobalParameters globalParameters() {
        return new AlleleFractionGlobalParameters(sampleDepth(), sampleBias(), outlierProbability());
    }

    public MinorFractions minorFractions() {
        return get(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, MinorFractions.class);
    }
}
