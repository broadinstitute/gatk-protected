package org.broadinstitute.hellbender.tools.tumorheterogeneity;

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
    private final List<Double> points;
    private final int numPoints;

    public TumorHeterogeneityData(final List<Double> points) {
        Utils.nonNull(points);
        this.points = new ArrayList<>(points);
        this.numPoints = points.size();
    }

    public int numPoints() {
        return numPoints;
    }

    public double getPoint(final int index) {
        return points.get(index);
    }
}