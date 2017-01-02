package org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Represents an allelic copy number state (m, n).
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class PloidyState {
    public static final PloidyState OUTLIER_PLOIDY_STATE = new PloidyState(Integer.MAX_VALUE, Integer.MAX_VALUE);

    private final int m;
    private final int n;
    private final int total;

    public PloidyState(final int m, final int n) {
        Utils.validateArg(m >= 0, "Number of allele copies must be non-negative.");
        Utils.validateArg(n >= 0, "Number of allele copies must be non-negative.");
        this.m = m;
        this.n = n;
        total = m + n;
    }

    public int m() {
        if (isOutlier()) {
            throw new IllegalStateException("Allele count for outlier state is undefined.");
        }
        return m;
    }

    public int n() {
        if (isOutlier()) {
            throw new IllegalStateException("Allele count for outlier state is undefined.");
        }
        return n;
    }

    public int total() {
        if (isOutlier()) {
            throw new IllegalStateException("Total count for outlier state is undefined.");
        }
        return total;
    }

    public boolean isOutlier() {
        return this == OUTLIER_PLOIDY_STATE;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final PloidyState that = (PloidyState) o;

        return m == that.m && n == that.n;

    }

    @Override
    public int hashCode() {
        return 31 * m + n;
    }

    @Override
    public String toString() {
        return String.format("%d_%d", m, n);
    }
}
