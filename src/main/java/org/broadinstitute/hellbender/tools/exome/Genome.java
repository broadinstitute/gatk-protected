package org.broadinstitute.hellbender.tools.exome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Represents copy-number data from a single individual for exome analysis.  Contains coverage data from targets
 * and ref/alt allele counts at normal germline het SNP sites.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class Genome {

    private final Logger logger = LogManager.getLogger(Genome.class);

    private final TargetCollection<ReadCountRecord.SingleSampleRecord> targets;
    private final TargetCollection<AllelicCount> snps;
    private final String sampleName;

    /**
     * Constructs a genome from lists containing linear target-coverage and SNP-allele-count data.
     * @param targets       list of linear target coverages, cannot be {@code null}
     * @param snps          list of SNP allele counts, cannot be {@code null}
     */
    public Genome(final ReadCountCollection targets, final List<AllelicCount> snps) {
        Utils.nonNull(targets, "The list of targets cannot be null.");
        Utils.nonNull(snps, "The list of SNPs cannot be null.");
        this.targets = new HashedListTargetCollection<>(targets.records().stream().map(ReadCountRecord::asSingleSampleRecord).collect(Collectors.toList()));
        this.snps = new HashedListTargetCollection<>(snps);
        sampleName = ReadCountCollectionUtils.getSampleNameFromReadCounts(targets);
    }

    /**
     * Constructs a genome from files containing linear target-coverage and SNP-allele-count data.
     * @param tangentNormalizedCoverageFile    linear target-coverage file
     * @param snpFile       SNP-allele-count file
     */
    public Genome(final File tangentNormalizedCoverageFile, final File snpFile) {
        logger.info("Constructing genome...");
        sampleName = ReadCountCollectionUtils.getSampleNameForCLIsFromReadCountsFile(tangentNormalizedCoverageFile);
        try {
            final ReadCountCollection rcc = ReadCountCollectionUtils.parse(tangentNormalizedCoverageFile);
            logger.info("Target counts from read count collection: " + rcc.targets().size());
            targets = new HashedListTargetCollection<>(rcc.records()
                    .stream().map(ReadCountRecord::asSingleSampleRecord).collect(Collectors.toList()));
            logger.info("Target counts from hashed list: " + targets.targets().size());
            snps = new HashedListTargetCollection<>(new AllelicCountCollection(snpFile).getCounts());
        } catch (final IOException e) {
            throw new UserException.BadInput("Could not read normalized coverage file");
        }
    }

    public final TargetCollection<ReadCountRecord.SingleSampleRecord> getTargets() {  return targets; }

    public final TargetCollection<AllelicCount> getSNPs() {  return snps; }

    public final String getSampleName() {   return sampleName;  }
}
