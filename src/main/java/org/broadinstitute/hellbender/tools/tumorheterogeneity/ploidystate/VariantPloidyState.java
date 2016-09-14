package org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class VariantPloidyState {
    private final int m;
    private final int n;
    private final int total;

    public VariantPloidyState(final int m, final int n) {
        Utils.validateArg(m >= 0, "Number of major-allele copies must be non-negative.");
        Utils.validateArg(n >= 0, "Number of minor-allele copies must be non-negative.");
        Utils.validateArg(m >= n, "Number of major-allele copies must be greater than or equal to number of minor-allele copies.");
        this.m = m;
        this.n = n;
        total = m + n;
    }

    public int m() {
        return m;
    }

    public int n() {
        return n;
    }

    public int total() {
        return total;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final VariantPloidyState that = (VariantPloidyState) o;

        return m == that.m && n == that.n;

    }

    @Override
    public int hashCode() {
        return 31 * m + n;
    }

    @Override
    public String toString() {
        return "(m = " + m +", n = " + n + ')';
    }
}
