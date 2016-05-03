package org.broadinstitute.hellbender.tools.exome;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.NotPrimaryAlignmentFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegrator;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegratorFactory;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A Bayesian heterozygous SNP pulldown calculator. Base qualities are taken into account
 * to increase precision (see CNV-methods.pdf for details).
 *
 * TODO
 *
 * - The quadrature order can be adaptively chosen based on the pileup size and the accuracy
 *   required for likelihood estimation. In theory, a Gaussian quadrature of order N yields
 *   the exact result for a pileup of size N.
 *
 * - Include the possibility to correct for reference bias for the HetPriorType.BALANCED prior
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */

public final class BayesianHetPulldownCalculator {

    private final Logger logger = LogManager.getLogger(BayesianHetPulldownCalculator.class);

    /**
     * a simple class to handle base read and mapping error probabilitites
     */
    @VisibleForTesting
    public static final class BaseQuality {

        final double readErrorProbability, mappingErrorProbability;

        public BaseQuality(final double readErrorProbability, final double mappingErrorProbability) {
            this.readErrorProbability = readErrorProbability;
            this.mappingErrorProbability = mappingErrorProbability;
        }

        public double getReadErrorProbability() { return readErrorProbability; }
        public double getMappingErrorProbability() { return mappingErrorProbability; }

    }

    /**
     * prior model for Het site calling
     *
     *  - BALANCED: minor allele fraction is assumed to be 1/2. this should be the default choice for paired
     *    normal-tumor samples
     *
     *  - HETEROGENEOUS: minor allele fraction can deviate from 1/2 due to CNV events, ploidy, subclonality,
     *    contamination, etc. this should be the default choice for tumor-only samples
     *
     */
    private enum HetPriorType {
        BALANCED, HETEROGENEOUS
    }

    private HetPriorType hetPriorType;

    private static final Nucleotide[] BASES = {Nucleotide.A, Nucleotide.C, Nucleotide.T, Nucleotide.G};

    private final File refFile;
    private final IntervalList snpIntervals;

    private final int readDepthThreshold;
    private final int minMappingQuality;
    private final int minBaseQuality;
    private final ValidationStringency validationStringency;

    /* experimental */
    private final double errorProbabilityAdjustmentFactor;


    /* these are for building a prior for allele fraction */
    private double minAbnormalFraction;
    private double maxAbnormalFraction;
    private double maxCopyNumber;
    private double minHetAlleleFraction;
    private double breakpointHetAlleleFraction;

    /* integration quadrature */
    @VisibleForTesting
    final ArrayList<Double> gaussIntegrationWeights = new ArrayList<>();
    @VisibleForTesting
    final ArrayList<Double> gaussIntegrationLogWeights = new ArrayList<>();
    @VisibleForTesting
    final ArrayList<Double> gaussIntegrationAbscissas = new ArrayList<>();

    /* allele fraction prior for Het sites */
    @VisibleForTesting
    final ArrayList<Double> alleleFractionPriors = new ArrayList<>();
    @VisibleForTesting
    final ArrayList<Double> alleleFractionLogPriors = new ArrayList<>();

    /* minimum order of the integration quadrature */
    private static final int MIN_QUADRATURE_ORDER = 50;

    /* interval threshold for indexing for SamLocusIterator */
    private static final int MAX_INTERVALS_FOR_INDEX = 25000;

    /* default priors */
    private static final double DEFAULT_PRIOR_REF_HOM = 0.5; /* a homozygous site being the ref allele */
    private static final double DEFAULT_PRIOR_HET = 0.5; /* a site being heterozygous */

    /* a third of the minimum sequencing error (for safeguarding log likelihood calculations) */
    private static final double DEFAULT_MIN_BASE_ERROR_THIRD = 1e-6;

    /* approximate number of status updates printed to log */
    private static final int NUMBER_OF_SITES_PER_LOGGED_STATUS_UPDATE = 10000;

    /**
     * Constructor of {@link BayesianHetPulldownCalculator} object
     *
     * NOTE: The default {@link HetPriorType} is BALANCED. Make a call to
     *       {@link BayesianHetPulldownCalculator#useHeterogeneousHetPrior} to switch to the HETEROGENEOUS prior.
     *
     * @param refFile the reference genome file
     * @param snpIntervals {@link IntervalList} of common SNPs
     * @param minMappingQuality minimum phred mapping quality
     * @param minBaseQuality minimum phred base quality
     * @param readDepthThreshold minimum read depth
     * @param validationStringency validation stringency
     * @param errorProbabilityAdjustmentFactor (experimental) multiplicative factor for read and mapping error
     *                                         probabilities
     */
    public BayesianHetPulldownCalculator(final File refFile, final IntervalList snpIntervals,
                                         final int minMappingQuality, final int minBaseQuality,
                                         final int readDepthThreshold, final ValidationStringency validationStringency,
                                         final double errorProbabilityAdjustmentFactor) {

        ParamUtils.isPositiveOrZero(minMappingQuality, "Minimum mapping quality must be nonnegative.");
        ParamUtils.isPositiveOrZero(minBaseQuality, "Minimum base quality must be nonnegative.");

        this.refFile = refFile;
        this.snpIntervals = snpIntervals;
        this.minMappingQuality = ParamUtils.isPositive(minMappingQuality, "Minimum mapping quality must be a positive integer");
        this.minBaseQuality = ParamUtils.isPositive(minBaseQuality, "Minimum base quality must be a positive integer");
        this.readDepthThreshold = ParamUtils.isPositive(readDepthThreshold, "Read depth threshold must be a positive integer");
        this.validationStringency = validationStringency;
        this.errorProbabilityAdjustmentFactor = ParamUtils.isPositive(errorProbabilityAdjustmentFactor,
                "Error adjustment factor must be positive.");

        /* the default Het prior */
        useBalancedHetPrior();
    }

    /**
     * use the BALANCED prior for Het sites.
     */
    public void useBalancedHetPrior() {
        hetPriorType = HetPriorType.BALANCED;
        logger.info("Using the BALANCED Het prior.");
    }

    /**
     * use the HETEROGENEOUS prior for Het sites.
     * @param minAbnormalFraction estimated minimum fraction of non-germline cells in the sample
     * @param maxAbnormalFraction estimated maximum fraction of non-germline cells in the sample
     * @param maxCopyNumber estimated maximum copy number in non-germline events (note: we use a flat probability
     *                      function for the copy number. using a large value of maxCopyNumber will result in a
     *                      significant spread of the minor allele fraction prior around 1/2. it is recommended not
     *                      to use values > 4).
     * @param quadratureOrder the order of quadrature used in numerical integrations
     */
    public void useHeterogeneousHetPrior(final double minAbnormalFraction, final double maxAbnormalFraction,
                                         final int maxCopyNumber, final int quadratureOrder) {

        hetPriorType = HetPriorType.HETEROGENEOUS;

        /* parameters for building the allele ratio prior at Het sites */
        this.minAbnormalFraction = ParamUtils.inRange(minAbnormalFraction, 0.0, 1.0, "Minimum fraction of abnormal" +
                " cells must be between 0 and 1.");
        this.maxAbnormalFraction = ParamUtils.inRange(maxAbnormalFraction, this.minAbnormalFraction, 1.0, "Maximum fraction of abnormal" +
                " cells must be greater than the provided minimum and less than 1.");
        this.maxCopyNumber = ParamUtils.isPositive(maxCopyNumber, "Maximum copy number must be positive");

        /* auxiliary member functions */
        this.minHetAlleleFraction = (1 - this.maxAbnormalFraction) / (this.maxCopyNumber * this.maxAbnormalFraction +
                2 * (1 - this.maxAbnormalFraction));
        this.breakpointHetAlleleFraction = (1 - this.minAbnormalFraction) / (this.maxCopyNumber * this.minAbnormalFraction +
                2 * (1 - this.minAbnormalFraction));

        /* initialize the integration quadrature and calculate the allele fraction prior on the abscissas */
        initializeIntegrationQuadrature(ParamUtils.isPositive(quadratureOrder - MIN_QUADRATURE_ORDER,
                "Quadrature order must be greater than " + MIN_QUADRATURE_ORDER) + MIN_QUADRATURE_ORDER);
        initializeHetAlleleFractionPrior();
        logger.info("Using the HETEROGENEOUS Het prior.");
    }

    /******************************
     * core computational methods *
     ******************************/

    /**
     * Initilizes the quadrature for calculating allele ratio integrals in getHetLogLikelihood
     * @param numIntegPoints  number of points in the quadrature
     */
    private void initializeIntegrationQuadrature(final int numIntegPoints) {

        /* get Gauss-Legendre quadrature factory of order @numIntegPoints */
        GaussIntegratorFactory integratorFactory = new GaussIntegratorFactory();
        GaussIntegrator gaussIntegrator = integratorFactory.legendre(numIntegPoints,
                minHetAlleleFraction, 1.0 - minHetAlleleFraction);

        /* abscissas */
        gaussIntegrationAbscissas.clear();
        gaussIntegrationAbscissas.addAll(IntStream.range(0, numIntegPoints).
                mapToDouble(gaussIntegrator::getPoint).boxed().collect(Collectors.toList()));

        /* weights */
        gaussIntegrationWeights.clear();
        gaussIntegrationWeights.addAll(IntStream.range(0, numIntegPoints).
                mapToDouble(gaussIntegrator::getWeight).boxed().collect(Collectors.toList()));

        /* log of weights */
        gaussIntegrationLogWeights.clear();
        gaussIntegrationLogWeights.addAll(gaussIntegrationWeights.stream().
                mapToDouble(FastMath::log).boxed().collect(Collectors.toList()));
    }

    /**
     * (advanced) calculate a simple prior (and log prior) for the allele fraction based on
     * (1) minimum purity of the sample, and (2) maximum copy number. See CNV-methods.pdf for details.
     */
    private void initializeHetAlleleFractionPrior() {

        alleleFractionPriors.clear();
        alleleFractionPriors.addAll(gaussIntegrationAbscissas.stream()
                .mapToDouble(alleleFraction -> {
                    double minorAlleleFraction = (alleleFraction < 0.5) ? alleleFraction : 1 - alleleFraction;
                    if (minorAlleleFraction < minHetAlleleFraction) {
                        return 0;
                    } else if (minorAlleleFraction < breakpointHetAlleleFraction) {
                        double denom = 2 * FastMath.pow((1 - minorAlleleFraction) * minorAlleleFraction * maxCopyNumber, 2) *
                                maxAbnormalFraction * (maxAbnormalFraction - minAbnormalFraction);
                        double num = (-1 + (-1 + minorAlleleFraction * maxCopyNumber) * maxAbnormalFraction) *
                                (-1 + maxAbnormalFraction + minorAlleleFraction * (2 + (-2 + maxCopyNumber) * maxAbnormalFraction)) +
                                2 * (1 + minorAlleleFraction * (-2 + minorAlleleFraction * maxCopyNumber)) * maxAbnormalFraction *
                                        FastMath.log(FastMath.abs(((1 + minorAlleleFraction * (-2 + maxCopyNumber)) * maxAbnormalFraction)) /
                                                (1 - 2 * minorAlleleFraction));
                        return num / denom;
                    } else { /* breakpointHetAlleleFraction < minorAlleleFraction < 1/2 */
                        double denom = 2 * FastMath.pow((1 - minorAlleleFraction) * minorAlleleFraction * maxCopyNumber, 2) *
                                maxAbnormalFraction * minAbnormalFraction * (maxAbnormalFraction - minAbnormalFraction);
                        double num = (maxAbnormalFraction - minAbnormalFraction) * (-1 + 2 * minorAlleleFraction +
                                (1 + minorAlleleFraction * (-2 + maxCopyNumber)) * (-1 + minorAlleleFraction * maxCopyNumber) *
                                maxAbnormalFraction * minAbnormalFraction) + 2 * (1 + minorAlleleFraction * (-2 + minorAlleleFraction *
                                maxCopyNumber)) * maxAbnormalFraction * minAbnormalFraction * FastMath.log(maxAbnormalFraction / minAbnormalFraction);
                        return num / denom;
                    }
                })
                .boxed().collect(Collectors.toList()));

        /* calculate the log prior */
        alleleFractionLogPriors.clear();
        alleleFractionLogPriors.addAll(alleleFractionPriors.stream()
                .mapToDouble(FastMath::log).boxed().collect(Collectors.toList()));
    }

    /**
     * calculate the log probability, safeguarded with a minimum probability DEFAULT_MIN_BASE_ERROR_THIRD
     * @param prob probability
     * @return safeguarded-log of probability
     */
    private double getSafeguardedLogProbability(final double prob) {
        return FastMath.log(FastMath.max(DEFAULT_MIN_BASE_ERROR_THIRD, prob));
    }

    /**
     * calculate the log likelihood of a SNP site being homozygous for a given read pileup
     * (see CNV-method.pdf for details)
     * @param baseQualities map of bases to list of their calling error probabilities at the SNP site
     * @param alleleRef the ref allele base
     * @param alleleAlt the alt allele base
     * @param homRefPrior the prior probability of the ref allele given that the site is homozygous
     * @return the log likelihood
     */
    @VisibleForTesting
    public double getHomLogLikelihood(final Map<Nucleotide, ArrayList<BaseQuality>> baseQualities,
                                      final Nucleotide alleleRef, final Nucleotide alleleAlt,
                                      final double homRefPrior) {

        /* initilize the log likehoods of ref and alt with the priors ... */
        double homRefLogLikelihood = FastMath.log(Utils.nonNull(homRefPrior));
        double homAltLogLikelihood = FastMath.log(1.0 - Utils.nonNull(homRefPrior));

        /* ... and add the log likelihood of the reads */
        for (Nucleotide base : baseQualities.keySet()) {
            for (BaseQuality currentBaseQuality : baseQualities.get(base)) {
                homRefLogLikelihood += getSafeguardedLogProbability(currentBaseQuality.getReadErrorProbability() / 3 +
                        currentBaseQuality.getMappingErrorProbability() / 4 +
                        ((base == alleleRef) ? (1 - 4 * currentBaseQuality.getReadErrorProbability() / 3
                                - currentBaseQuality.getMappingErrorProbability()) : 0));
                homAltLogLikelihood += getSafeguardedLogProbability(currentBaseQuality.getReadErrorProbability() / 3 +
                                currentBaseQuality.getMappingErrorProbability() / 4 +
                                ((base == alleleAlt) ? (1 - 4 * currentBaseQuality.getReadErrorProbability() / 3
                                        - currentBaseQuality.getMappingErrorProbability()) : 0));
            }
        }

        /* return the sum of |hom,ref) and |hom,alt) likelihoods */
        return GATKProtectedMathUtils.naturalLogSumExp(homRefLogLikelihood, homAltLogLikelihood);
    }

    /**
     * [internal helper function]
     * calculate the log likelihood of hetrozygosity from just alt and ref pileup for a given @alleleFraction
     * (see CNV-method.pdf for details)
     * @param alleleFraction ref-to-alt allele fraction
     * @param alphaList [internal to getHetLogLikelihood]
     * @param betaList [internal to getHetLogLikelihood]
     * @return log likelihood
     */
    private double getRefAltHetLogLikelihoodFixedAlleleFraction(final double alleleFraction,
                                                                final ArrayList<Double> alphaList,
                                                                final ArrayList<Double> betaList) {
        double logLikelihood = 0;
        for (int i=0; i<alphaList.size(); i++) {
            logLikelihood += getSafeguardedLogProbability(alphaList.get(i) + alleleFraction * betaList.get(i));
        }

        return logLikelihood;
    }

    /**
     * calculate the log likelihood of a SNP site being heterozygous for a given read pileup
     * (see CNV-method.pdf for details)
     * @param baseQualities map of bases to list of their calling error probabilities at the SNP site
     * @param alleleRef the ref allele base
     * @param alleleAlt the alt allele base
     * @return the log likelihood
     */
    @VisibleForTesting
    public double getHetLogLikelihood(final Map<Nucleotide, ArrayList<BaseQuality>> baseQualities,
                                       final Nucleotide alleleRef, final Nucleotide alleleAlt) {

        double errorLogLikelihood = 0.0;
        double refAltLogLikelihood;
        ArrayList<Double> alphaList = new ArrayList<>();
        ArrayList<Double> betaList = new ArrayList<>();

        for (Nucleotide base : baseQualities.keySet()) {
            if (base == alleleRef) {
                for (BaseQuality currentBaseQuality : baseQualities.get(base)) {
                    alphaList.add(currentBaseQuality.getReadErrorProbability() / 3 +
                        currentBaseQuality.getMappingErrorProbability() / 4);
                    betaList.add(1 - 4 * currentBaseQuality.getReadErrorProbability() / 3 -
                        currentBaseQuality.getMappingErrorProbability() / 4);
                }
            }
            else if (base == alleleAlt) {
                for (BaseQuality currentBaseQuality : baseQualities.get(base)) {
                    alphaList.add(1 - currentBaseQuality.getReadErrorProbability() -
                        3 * currentBaseQuality.getMappingErrorProbability() / 4);
                    betaList.add(-1 + 4 * currentBaseQuality.getReadErrorProbability() / 3 +
                        currentBaseQuality.getMappingErrorProbability());
                }
            }
            else {
                for (BaseQuality currentBaseQuality : baseQualities.get(base)) {
                    errorLogLikelihood += getSafeguardedLogProbability(currentBaseQuality.getReadErrorProbability() / 3 +
                        currentBaseQuality.getMappingErrorProbability() / 4);
                }
            }
        }

        switch (hetPriorType) {

            case BALANCED:
                /* set the allele fraction to 1/2 */
                refAltLogLikelihood = getRefAltHetLogLikelihoodFixedAlleleFraction(0.5, alphaList, betaList);
                break;

            case HETEROGENEOUS:
                /* marginalize allele fraction by integrating {@link getRefAltHetLogLikelihoodFixedAlleleFraction} over
                   the precomputed allele fraction prior */
                ArrayList<Double> refAltLogLikelihoodList = new ArrayList<>(gaussIntegrationAbscissas.size());
                refAltLogLikelihoodList.addAll(
                        gaussIntegrationAbscissas.stream()
                                .mapToDouble(Double::doubleValue)
                                .map(f -> getRefAltHetLogLikelihoodFixedAlleleFraction(f, alphaList, betaList))
                                .boxed().collect(Collectors.toList()));
                ArrayList<Double> logLikelihoodIntegrandWithPriorAndWeights = new ArrayList<>(gaussIntegrationAbscissas.size());
                logLikelihoodIntegrandWithPriorAndWeights.addAll(
                        IntStream.range(0, gaussIntegrationAbscissas.size())
                                .mapToDouble(i -> refAltLogLikelihoodList.get(i) + gaussIntegrationLogWeights.get(i) +
                                        alleleFractionLogPriors.get(i))
                                .boxed().collect(Collectors.toList()));
                refAltLogLikelihood = GATKProtectedMathUtils.naturalLogSumExp(logLikelihoodIntegrandWithPriorAndWeights);
                break;

            default: /* we shouldn't be here anyway */
                refAltLogLikelihood = 0.0;

        }

        return errorLogLikelihood + refAltLogLikelihood;
    }

    /**
     * wrapper methods
     */

    /**
     * Returns map of base-pair to error probabilities at a given locus. All reads are considered (not just ACTG)
     * @param locus locus
     * @return map of base-pair to error probabilities
     */
    private Map<Nucleotide, ArrayList<BaseQuality>> getPileupBaseQualities(final SamLocusIterator.LocusInfo locus) {

        Map<Nucleotide, ArrayList<BaseQuality>> baseQualities = new HashMap<>();

        /* make sure that we have keys at least for A, C, T, and G */
        for (Nucleotide base : BASES) {
            baseQualities.put(base, new ArrayList<>());
        }

        for (final SamLocusIterator.RecordAndOffset rec : locus.getRecordAndPositions()) {
            Nucleotide currentBase = Nucleotide.valueOf(rec.getReadBase());
            ArrayList<BaseQuality> currentBaseQualityList = baseQualities.get(currentBase);
            if (currentBaseQualityList == null) { /* this can occur if a non-ACTG read is found */
                currentBaseQualityList = new ArrayList<>();
                baseQualities.put(currentBase, currentBaseQualityList);
            }
            currentBaseQualityList.add(new BaseQuality(
                    errorProbabilityAdjustmentFactor * QualityUtils.qualToErrorProb(rec.getBaseQuality()),
                    errorProbabilityAdjustmentFactor * QualityUtils.qualToErrorProb(rec.getRecord().getMappingQuality()))
            );
        }

        return baseQualities;
    }

    /**
     * Returns base-pair counts at a given locus.
     * @param locus locus
     * @return      base-pair counts
     */
    static Nucleotide.Counter getPileupBaseCounts(final SamLocusIterator.LocusInfo locus) {
        final Nucleotide.Counter result = new Nucleotide.Counter();
        for (final SamLocusIterator.RecordAndOffset rec : locus.getRecordAndPositions()) {
            result.add(rec.getReadBase());
        }
        return result;
    }

    /**
     * get base counts map from base quality map
     * @param baseQualities map from bases to list of {@link BaseQuality}
     * @return base counts map
     */
    private Map<Nucleotide, Integer> getBaseCountsFromBaseQualities(final Map<Nucleotide, ArrayList<BaseQuality>> baseQualities) {
        Map<Nucleotide, Integer> baseCounts = new HashMap<>();
        for (Nucleotide base : BASES) {
            if (baseQualities.containsKey(base)) {
                baseCounts.put(base, baseQualities.get(base).size());
            } else {
                baseCounts.put(base, 0);
            }
        }
        return baseCounts;
    }

    /**
     * Infer the alt allele base from the pileup. We choose the base with the highest frequency
     * that is not Ref (this coincides with the maximum likelihood estimation).
     *
     * * Remark: while this is definitely not the perfect way to do it, it has no serious pitfalls:
     *
     * - if the input is balanced between betweet ref and alt, it correctly chooses alt.
     * - if the input is all ref, then alt is chosen randomly, but it's OK since there are no alt counts
     *   in the pileup anyway.
     * - if the input is all alt, we are good.
     * - if the input in all ref + a few errors, then alt will be choosen as the most frequent erroneous
     *   base; again, this is OK because Hom likelihood is far greater than Het likelihood in this case
     * - if the input is all alt + a few errors, we are good.

     * @param baseQualities map from bases to error probabilities
     * @param refBase the ref allele base
     * @return the likely alt allele
     */
    @VisibleForTesting
    public static Nucleotide inferAltFromPileup(final Map<Nucleotide, ArrayList<BaseQuality>> baseQualities,
                                                final Nucleotide refBase) {

        /* sort the bases in the descending order by their frequency */
        Nucleotide[] bases = {Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.T};
        Arrays.sort(bases, (L, R) -> Integer.compare(baseQualities.get(R).size(), baseQualities.get(L).size()));
        /* pick the base with highest frequency, skip over ref */
        for (Nucleotide base : bases) {
            if (base != refBase) {
                return base;
            }
        }

        /* we shouldn't be here unless the baseErrorProbabilities is malformed */
        return Nucleotide.N;
    }

    /**
     * For a normal or tumor sample, returns a data structure giving (intervals, reference counts, alternate counts),
     * where intervals give positions of likely heterozygous SNP sites.
     *
     * The analysis for normal and tumor is identical. The IntervalList gives common SNP sites in 1-based
     * format.
     *
     * The @hetCallingStrigency parameters sets the threshold Het posterior for calling:
     *
     *      hetPosteriorThreshold = 1 - 10^{-hetCallingStringency}
     *      hetThresholdLogOdds = log(hetPosteriorThreshold/(1-hetPosteriorThreshold))
     *                          = log(10^{hetCallingStringency} - 1)
     *
     * @param bamFile sorted BAM file for sample
     * @param hetCallingStringency strigency for calling a Het site
     * @return Pulldown of heterozygous SNP sites in 1-based format
     */
    @VisibleForTesting
    public Pulldown getHetPulldown(final File bamFile, final double hetCallingStringency) {

        /* log odds from stringency */
        final double hetThresholdLogOdds = FastMath.log(FastMath.pow(10, hetCallingStringency) - 1);

        try (final SamReader bamReader = SamReaderFactory.makeDefault().validationStringency(validationStringency)
                .referenceSequence(refFile).open(bamFile);
             final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(this.refFile)) {
            if (bamReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new UserException.BadInput("BAM file " + bamFile.toString() + " must be coordinate sorted.");
            }

            final Pulldown hetPulldown = new Pulldown(bamReader.getFileHeader());

            final int totalNumberOfSNPs = snpIntervals.size();
            final SamLocusIterator locusIterator = new SamLocusIterator(bamReader, snpIntervals,
                    totalNumberOfSNPs < MAX_INTERVALS_FOR_INDEX);

            /* set read and locus filters */
            final List<SamRecordFilter> samFilters = Arrays.asList(new NotPrimaryAlignmentFilter(),
                    new DuplicateReadFilter());
            locusIterator.setSamFilters(samFilters);
            locusIterator.setEmitUncoveredLoci(false);
            locusIterator.setIncludeNonPfReads(false);
            locusIterator.setMappingQualityScoreCutoff(minMappingQuality);
            locusIterator.setQualityScoreCutoff(minBaseQuality);

            logger.info("Examining " + totalNumberOfSNPs + " sites in total...");
            int locusCount = 0;
            for (final SamLocusIterator.LocusInfo locus : locusIterator) {
                if (locusCount % NUMBER_OF_SITES_PER_LOGGED_STATUS_UPDATE == 0) {
                    logger.info("Examined " + locusCount + " covered sites.");
                }
                locusCount++;

                final int totalReadCount = locus.getRecordAndPositions().size();
                if (totalReadCount <= readDepthThreshold) {
                    continue;
                }

                final Map<Nucleotide, ArrayList<BaseQuality>> baseQualities = getPileupBaseQualities(locus);
                final Nucleotide refBase = Nucleotide.valueOf(refWalker.get(locus.getSequenceIndex()).getBases()[locus.getPosition() - 1]);
                final Nucleotide altBase = inferAltFromPileup(baseQualities, refBase);

                /* calculate Het log odds */
                final double hetLogLikelihood = getHetLogLikelihood(baseQualities, refBase, altBase);
                final double homLogLikelihood = getHomLogLikelihood(baseQualities, refBase, altBase,
                        DEFAULT_PRIOR_REF_HOM);
                final double hetLogOdds = (hetLogLikelihood + FastMath.log(DEFAULT_PRIOR_HET)) -
                        (homLogLikelihood + FastMath.log(1 - DEFAULT_PRIOR_HET));

                if (hetLogOdds > hetThresholdLogOdds) {
                    hetPulldown.add(new SimpleInterval(locus.getSequenceName(), locus.getPosition(), locus.getPosition()),
                            baseQualities.get(refBase).size(), baseQualities.get(altBase).size(),
                            refBase, altBase);
                }
            }

            logger.info(locusCount + " covered sites out of " + totalNumberOfSNPs + " total sites were examined.");

            return hetPulldown;

        } catch (final IOException | SAMFormatException e) {
            throw new UserException(e.getMessage());
        }
    }

    public Pulldown getTumorHetPulldownFromNormalPulldown(final File tumorBamFile, final Pulldown normalHetPulldown) {

        try (final SamReader bamReader = SamReaderFactory.makeDefault().validationStringency(validationStringency)
                .referenceSequence(refFile).open(tumorBamFile)) {
            if (bamReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new UserException.BadInput("BAM file " + tumorBamFile.toString() + " must be coordinate sorted.");
            }

            final Pulldown tumorHetPulldown = new Pulldown(bamReader.getFileHeader());

            final int totalNumberOfSNPs = snpIntervals.size();
            final SamLocusIterator locusIterator = new SamLocusIterator(bamReader, snpIntervals,
                    totalNumberOfSNPs < MAX_INTERVALS_FOR_INDEX);

            /* set read and locus filters */
            final List<SamRecordFilter> samFilters = Arrays.asList(new NotPrimaryAlignmentFilter(),
                    new DuplicateReadFilter());
            locusIterator.setSamFilters(samFilters);
            locusIterator.setEmitUncoveredLoci(false);
            locusIterator.setIncludeNonPfReads(false);
            locusIterator.setMappingQualityScoreCutoff(minMappingQuality);
            locusIterator.setQualityScoreCutoff(minBaseQuality);

            /* get a map of SimpleIntervals in the pulldown to their index */
            Map<SimpleInterval, Integer> normalPulldownIndexMap = normalHetPulldown.getSimpleIntervalToIndexMap();

            logger.info("Examining " + totalNumberOfSNPs + " sites in total...");
            int locusCount = 0;
            for (final SamLocusIterator.LocusInfo locus : locusIterator) {
                if (locusCount % NUMBER_OF_SITES_PER_LOGGED_STATUS_UPDATE == 0) {
                    logger.info("Examined " + locusCount + " covered sites.");
                }
                locusCount++;

                final int totalReadCount = locus.getRecordAndPositions().size();
                if (totalReadCount <= readDepthThreshold) {
                    continue;
                }

                /* find the AllelicCount from the normal pulldown */
                int indexInNormalPulldown;
                try {
                    indexInNormalPulldown = normalPulldownIndexMap.get(
                            new SimpleInterval(locus.getSequenceName(), locus.getPosition(), locus.getPosition()));
                } catch (NullPointerException e) {
                    throw new GATKException.ShouldNeverReachHereException("Can not find the required AllelicCount " +
                            "object in the normal pulldown. Stopping.");
                }

                /* just count the alt and ref nucleotide and add to the tumor pulldown */
                final Nucleotide.Counter baseCounts = getPileupBaseCounts(locus);
                tumorHetPulldown.add(new SimpleInterval(locus.getSequenceName(), locus.getPosition(), locus.getPosition()),
                        (int) baseCounts.get(normalHetPulldown.getCounts().get(indexInNormalPulldown).getRefNucleotide()),
                        (int) baseCounts.get(normalHetPulldown.getCounts().get(indexInNormalPulldown).getAltNucleotide()));
            }

            logger.info(locusCount + " covered sites out of " + totalNumberOfSNPs + " total sites were examined.");

            return tumorHetPulldown;

        } catch (final IOException | SAMFormatException e) {
            throw new UserException(e.getMessage());
        }
    }

}