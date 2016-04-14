package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.acnvconversion.ACNVModeledSegmentConversionUtils;
import org.broadinstitute.hellbender.tools.exome.acsconversion.ACSModeledSegmentUtils;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionModeller;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents an ACNV segmented model for copy ratio and allele fraction.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ACNVModeller {
    //use 95% HPD interval to construct {@link PosteriorSummary} for segment means and minor allele fractions
    private static final double CREDIBLE_INTERVAL_ALPHA = 0.05;

    public static final Logger logger = LogManager.getLogger(ACNVModeller.class);

    private final SegmentedGenome targetSegmentedGenome;
    private final SegmentedGenome snpSegmentedGenome;
    private final Genome genome;
    private final AllelicPanelOfNormals allelicPON;
    private final List<ACNVModeledSegment> segments = new ArrayList<>();

    private final int numSamplesCopyRatio;
    private final int numBurnInCopyRatio;
    private final int numSamplesAlleleFraction;
    private final int numBurnInAlleleFraction;
    private final JavaSparkContext ctx;

    public List<ACNVModeledSegment> getACNVModeledSegments() {
        return Collections.unmodifiableList(segments);
    }

    /**
     * Constructs a copy-ratio and allele-fraction modeller for a {@link SegmentedGenome},
     * specifying number of total samples and number of burn-in samples for Markov-Chain Monte Carlo model fitting.
     * An initial model fit is performed.
     *
     * @param targetSegmentedGenome     contains target segments, target coverages, and SNP counts for modelling copy ratio
     * @param snpSegmentedGenome        contains SNP segments, target coverages, and SNP counts for modelling allele fraction
     * @param numSamplesCopyRatio       number of total samples for copy-ratio model MCMC
     * @param numBurnInCopyRatio        number of burn-in samples to discard for copy-ratio model MCMC
     * @param numSamplesAlleleFraction  number of total samples for allele-fraction model MCMC
     * @param numBurnInAlleleFraction   number of burn-in samples to discard for allele-fraction model MCMC
     * @param ctx                       JavaSparkContext, used for kernel density estimation in {@link PosteriorSummary}
     */
    public ACNVModeller(final SegmentedGenome targetSegmentedGenome, final SegmentedGenome snpSegmentedGenome,
                        final int numSamplesCopyRatio, final int numBurnInCopyRatio,
                        final int numSamplesAlleleFraction, final int numBurnInAlleleFraction,
                        final JavaSparkContext ctx) {
        this(targetSegmentedGenome, snpSegmentedGenome, AllelicPanelOfNormals.EMPTY_PON, numSamplesCopyRatio, numBurnInCopyRatio, numSamplesAlleleFraction, numBurnInAlleleFraction, ctx);
    }

    /**
     * Constructs a copy-ratio and allele-fraction modeller for a {@link SegmentedGenome},
     * specifying number of total samples and number of burn-in samples for Markov-Chain Monte Carlo model fitting.
     * An initial model fit is performed.
     *
     * @param targetSegmentedGenome     contains target segments, target coverages, and SNP counts for modelling copy ratio
     * @param snpSegmentedGenome        contains SNP segments, target coverages, and SNP counts for modelling allele fraction
     * @param allelicPON                allelic-bias panel of normals
     * @param numSamplesCopyRatio       number of total samples for copy-ratio model MCMC
     * @param numBurnInCopyRatio        number of burn-in samples to discard for copy-ratio model MCMC
     * @param numSamplesAlleleFraction  number of total samples for allele-fraction model MCMC
     * @param numBurnInAlleleFraction   number of burn-in samples to discard for allele-fraction model MCMC
     * @param ctx                       JavaSparkContext, used for kernel density estimation in {@link PosteriorSummary}
     */
    public ACNVModeller(final SegmentedGenome targetSegmentedGenome, final SegmentedGenome snpSegmentedGenome, final AllelicPanelOfNormals allelicPON,
                        final int numSamplesCopyRatio, final int numBurnInCopyRatio,
                        final int numSamplesAlleleFraction, final int numBurnInAlleleFraction,
                        final JavaSparkContext ctx) {
        this.targetSegmentedGenome = targetSegmentedGenome;
        this.snpSegmentedGenome = snpSegmentedGenome;
        if (!targetSegmentedGenome.getGenome().equals(snpSegmentedGenome.getGenome())) {
            throw new IllegalArgumentException("Target-segmented genome and SNP-segmented genome must contain identical genome data.");
        }
        this.genome = targetSegmentedGenome.getGenome();
        this.allelicPON = allelicPON;
        this.numSamplesCopyRatio = numSamplesCopyRatio;
        this.numBurnInCopyRatio = numBurnInCopyRatio;
        this.numSamplesAlleleFraction = numSamplesAlleleFraction;
        this.numBurnInAlleleFraction = numBurnInAlleleFraction;
        this.ctx = ctx;
        logger.info("Fitting initial model...");
        fitModel();
    }

    /**
     * Performs Markov-Chain Monte Carlo model fitting using the
     * number of total samples and number of burn-in samples specified at construction.
     */
    public void fitModel() {
        //perform MCMC to generate posterior samples
        logger.info("Fitting copy-ratio model...");
        final ACNVCopyRatioModeller copyRatioModeller = new ACNVCopyRatioModeller(targetSegmentedGenome);
        copyRatioModeller.fitMCMC(numSamplesCopyRatio, numBurnInCopyRatio);
        logger.info("Fitting allele-fraction model...");
        final AlleleFractionModeller alleleFractionModeller = new AlleleFractionModeller(snpSegmentedGenome, allelicPON);
        alleleFractionModeller.fitMCMC(numSamplesAlleleFraction, numBurnInAlleleFraction);

        //update list of ACNVModeledSegment with new PosteriorSummaries
        final List<PosteriorSummary> segmentMeanPosteriorSummaries =
                copyRatioModeller.getSegmentMeanPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);
        final List<PosteriorSummary> minorAlleleFractionPosteriorSummaries =
                alleleFractionModeller.getMinorAlleleFractionPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);
        
        final TargetCollection<LocatablePosteriorSummary> segmentMeanPosteriorSummariesMap =
                makePosteriorSummariesMap(targetSegmentedGenome.getSegments(), segmentMeanPosteriorSummaries);
        final TargetCollection<LocatablePosteriorSummary> minorAlleleFractionPosteriorSummariesMap =
                makePosteriorSummariesMap(snpSegmentedGenome.getSegments(), minorAlleleFractionPosteriorSummaries);
        
        updateACNVModeledSegments(targetSegmentedGenome.getSegments(), snpSegmentedGenome.getSegments(), 
                segmentMeanPosteriorSummariesMap, minorAlleleFractionPosteriorSummariesMap, genome);
    }
    /**
     * Writes the list of {@link ACNVModeledSegment} held internally to file.
     * See {@link SegmentUtils#writeACNVModeledSegmentFile}.
     * @param outFile   output file
     */
    public void writeACNVModeledSegmentFile(final File outFile) {
        SegmentUtils.writeACNVModeledSegmentFile(outFile, segments, genome);
    }

    /**
     * Writes the list of {@link ACNVModeledSegment} held internally to file as a list of {@link ModeledSegment}.
     * See {@link SegmentUtils#writeModeledSegmentFile}.
     * @param outFile   output file
     */
    public void writeModeledSegmentFile(final File outFile) {
        SegmentUtils.writeModeledSegmentFile(outFile,
                ACNVModeledSegmentConversionUtils.convertACNVModeledSegmentsToModeledSegments(segments, genome),
                genome.getSampleName());
    }

    public void writeAllelicCapSegFile(final File outFile) {
        ACSModeledSegmentUtils.writeACNVModeledSegmentFileAsAllelicCapSegFile(outFile, segments, genome);
    }


    private static TargetCollection<LocatablePosteriorSummary> makePosteriorSummariesMap(final List<SimpleInterval> segments,
                                                                                         final List<PosteriorSummary> posteriorSummaries) {
        final List<LocatablePosteriorSummary> locatablePosteriorSummaries =
                IntStream.range(0, segments.size()).boxed()
                        .map(i -> new LocatablePosteriorSummary(segments.get(i), posteriorSummaries.get(i)))
                        .collect(Collectors.toList());
        return new HashedListTargetCollection<>(locatablePosteriorSummaries);
    }

    private void updateACNVModeledSegments(final List<SimpleInterval> targetSegments,
                                           final List<SimpleInterval> snpSegments,
                                           final TargetCollection<LocatablePosteriorSummary> segmentMeanPosteriorSummariesMap,
                                           final TargetCollection<LocatablePosteriorSummary> minorAlleleFractionPosteriorSummariesMap,
                                           final Genome genome) {
        //create unioned segments from naive combination of target and SNP breakpoints
        final List<SimpleInterval> unionedUnmodeledSegments = SegmentUtils.unionSegmentsNaively(targetSegments, snpSegments, genome);
        final List<ACNVModeledSegment> updatedSegments =
                unionedUnmodeledSegments.stream()
                        .map(segment -> {
                            //for each naive unioned segment, take posterior summaries of segment mean and minor-allele fraction when available,
                            //otherwise take posterior summary with all NaN entries
                            final LocatablePosteriorSummary segmentMeanLocatablePosteriorSummary = segmentMeanPosteriorSummariesMap.target(segment);
                            final LocatablePosteriorSummary minorAlleleFractionLocatablePosteriorSummary = minorAlleleFractionPosteriorSummariesMap.target(segment);
                            return new ACNVModeledSegment(
                                    segment,
                                    segmentMeanLocatablePosteriorSummary == null ? PosteriorSummary.NAN_POSTERIOR_SUMMARY : segmentMeanLocatablePosteriorSummary.posteriorSummary,
                                    minorAlleleFractionLocatablePosteriorSummary == null ? PosteriorSummary.NAN_POSTERIOR_SUMMARY : minorAlleleFractionLocatablePosteriorSummary.posteriorSummary);
                        })
                        .collect(Collectors.toList());
        segments.clear();
        segments.addAll(updatedSegments);
    }

    private static class LocatablePosteriorSummary implements Locatable {
        final SimpleInterval segment;
        final PosteriorSummary posteriorSummary;
        
        public LocatablePosteriorSummary(final SimpleInterval segment, final PosteriorSummary posteriorSummary) {
            this.segment = segment;
            this.posteriorSummary = posteriorSummary;
        }

        @Override
        public String getContig() {
            return segment.getContig();
        }

        @Override
        public int getStart() {
            return segment.getStart();
        }

        @Override
        public int getEnd() {
            return segment.getEnd();
        }
    }
}
