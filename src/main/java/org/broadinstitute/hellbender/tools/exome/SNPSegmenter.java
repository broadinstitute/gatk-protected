package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.segmenter.RCBSSegmenter;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;


public final class SNPSegmenter {
    //@staticmethod
//def create_snp_based_segments(dataframe, allelic_fraction_skew, sample_name):
//        """
//        This function creates segments from SNPs. Log copy ratios are simulated using allelic fractions.
//
//        :param dataframe: data frame with SNPs
//        :param allelic_fraction_skew: initial allelic fraction skew
//        :param sample_name: sample name
//        :return: SNPs based segments
//        """
//        dataframe["stop"] = dataframe["Start_position"]
//        dataframe = dataframe.rename(columns={"Chromosome": "contig", "Start_position": "start"})  # rename
//        dataframe[sample_name] = SNPSegmenter._transform_allelic_fractions(reference_counts=dataframe["t_ref_count"],
//        alternate_counts=dataframe["t_alt_count"],
//        allelic_fraction_skew=allelic_fraction_skew)
//        column_names = ["contig", "start", "stop", sample_name]
//        dataframe = dataframe[column_names]
//        segmenter = Segmenter()
//
//        segments = segmenter.input_and_segment_data(df=dataframe, sample_name=sample_name)
//        return segments
    public static File findSegments(final List<AllelicCount> snpCounts, final double allelicFractionSkew,
                                    final String sampleName, final File outputFile, final float minLogValue)
            throws IOException {
        final File targetsFromSNPAllelicCountsFile = File.createTempFile("targets-from-snps", ".tsv");

        return RCBSSegmenter.segment(sampleName, targetsFromSNPAllelicCountsFile.getAbsolutePath(),
                outputFile.getAbsolutePath(), minLogValue);
    }


>>>>>>> b05c904... Refactored Pulldown as a list of AllelicCounts. Made minor fixes to tests for GetHetCoverage and HetPulldownCalculator.
}
