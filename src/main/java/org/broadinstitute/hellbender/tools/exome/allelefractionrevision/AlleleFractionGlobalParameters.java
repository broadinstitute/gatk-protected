package org.broadinstitute.hellbender.tools.exome.allelefractionrevision;

/**
 * Encapsulates the global parameters of the allele fraction model: the mean and variance of the common prior on
 * allelic biases and the outlier probability.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionGlobalParameters {
    private final double sampleDepth;
    private final double sampleBias;
    private final double outlierProbability;

    public AlleleFractionGlobalParameters(final double sampleDepth, final double sampleBias, final double outlierProbability) {
        this.sampleDepth = sampleDepth;
        this.sampleBias = sampleBias;
        this.outlierProbability = outlierProbability;
    }

    public double getSampleDepth() {
        return sampleDepth;
    }

    public double getSampleBias() {
        return sampleBias;
    }

    public double getOutlierProbability() {
        return outlierProbability;
    }

    public AlleleFractionGlobalParameters copyWithNewSampleDepth(final double newSampleDepth) {
        return new AlleleFractionGlobalParameters(newSampleDepth, sampleBias, outlierProbability);
    }

    public AlleleFractionGlobalParameters copyWithNewSampleBias(final double newSampleBias) {
        return new AlleleFractionGlobalParameters(sampleDepth, newSampleBias, outlierProbability);
    }

    public AlleleFractionGlobalParameters copyWithNewOutlierProbability(final double newOutlierProbability) {
        return new AlleleFractionGlobalParameters(sampleDepth, sampleBias, newOutlierProbability);
    }
}
