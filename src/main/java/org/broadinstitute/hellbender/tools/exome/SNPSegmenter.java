package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.segmenter.RCBSSegmenter;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Finds segments based on allelic counts at SNP sites.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SNPSegmenter {
    /**
     * Writes (as a side effect) and returns segment file based on allelic counts at SNP sites.  Converts allelic counts
     * to target coverages, which are written to a file and then passed to {@link RCBSSegmenter}.
     * @param snpCounts             list of allelic counts at SNP sites
     * @param allelicFractionSkew   allelic fraction skew (1.0 in the ideal case)
     * @param sampleName            sample name
     * @param outputFile            segment file to write to and return
     * @param minLogValue           values of log2 copy ratio below this minimum value are set to it
     * @return                      segment file
     * @throws IOException
     */
    public static File findSegments(final List<AllelicCount> snpCounts, final double allelicFractionSkew,
                                    final String sampleName, final File outputFile, final float minLogValue)
            throws IOException {
        final File targetsFromSNPCountsFile = File.createTempFile("targets-from-snps", ".tsv");

        List<TargetCoverage> targetsFromSNPCounts = snpCounts.stream()
                .map(count -> count.asTargetCoverage("snp-target", allelicFractionSkew)).collect(Collectors.toList());

        TargetCoverage.writeTargetsWithCoverage(targetsFromSNPCountsFile, sampleName, targetsFromSNPCounts);

        return RCBSSegmenter.segment(sampleName, targetsFromSNPCountsFile.getAbsolutePath(),
                outputFile.getAbsolutePath(), minLogValue);
    }
}
