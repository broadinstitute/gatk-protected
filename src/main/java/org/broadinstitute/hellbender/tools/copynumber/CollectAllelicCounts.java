package org.broadinstitute.hellbender.tools.copynumber;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollector;
import org.broadinstitute.hellbender.utils.Nucleotide;

import java.io.File;
import java.util.List;

/**
 * Collects reference/alternate allele counts at sites.  The alt count is defined as the total count minus the ref count,
 * and the alt nucleotide is defined as the non-ref base with the highest count, with ties broken by the order of the
 * bases in {@link AllelicCountCollector#BASES}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Collects ref/alt counts at sites.",
        oneLineSummary = "Collects ref/alt counts at sites",
        programGroup = CopyNumberProgramGroup.class
)
public final class CollectAllelicCounts extends LocusWalker {

    private static final Logger logger = LogManager.getLogger(CollectAllelicCounts.class);

    @Argument(
            doc = "Output allelic-counts file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false
    )
    protected File outputAllelicCountsFile;

    private static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 30;

    @Argument(
            doc = "Minimum base quality; base calls with lower quality will be filtered out of pileup.",
            fullName = "minimumBaseQuality",
            shortName = "minBQ",
            optional = true,
            minValue = 0
    )
    protected int minimumBaseQuality = 20;

    private AllelicCountCollector allelicCountCollector = new AllelicCountCollector();

    @Override
    public void onTraversalStart() {
        logger.info("Collecting allelic counts...");
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> initialReadFilters = super.getDefaultReadFilters();
        initialReadFilters.add(new MappingQualityReadFilter(DEFAULT_MINIMUM_MAPPING_QUALITY));
        return initialReadFilters;
    }

    @Override
    public Object onTraversalSuccess() {

        allelicCountCollector.getAllelicCounts().write(outputAllelicCountsFile);
        logger.info("Allelic counts written to " + outputAllelicCountsFile.toString());
        return("SUCCESS");
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {

        final byte refAsByte = referenceContext.getBase();
        allelicCountCollector.collectAtLocus(Nucleotide.valueOf(refAsByte), alignmentContext.getBasePileup(), alignmentContext.getLocation(), minimumBaseQuality);
    }
}