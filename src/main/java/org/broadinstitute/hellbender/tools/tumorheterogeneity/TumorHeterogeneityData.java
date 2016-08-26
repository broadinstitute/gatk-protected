package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;

import java.util.ArrayList;
import java.util.List;

/**
 * {@link DataCollection} for the tumor-heterogeneity model.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityData implements DataCollection {
    private final List<ACNVModeledSegment> segments;
    private final int numBinsCopyRatio;
    private final int numBinsAlleleFraction;

    public TumorHeterogeneityData(final List<ACNVModeledSegment> segments,
                                  final int numBinsCopyRatio,
                                  final int numBinsAlleleFraction) {
        Utils.nonNull(segments);
        Utils.validateArg(numBinsCopyRatio > 0, "Number of copy-ratio bins must be positive.");
        Utils.validateArg(numBinsAlleleFraction > 0, "Number of allele-fraction bins must be positive.");
        this.segments = new ArrayList<>(segments);
        this.numBinsCopyRatio = numBinsCopyRatio;
        this.numBinsAlleleFraction = numBinsAlleleFraction;

    }
}