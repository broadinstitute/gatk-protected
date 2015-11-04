package org.broadinstitute.hellbender.tools.exome;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.BroadcastJoinReadsWithRefBases;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Gets alt/ref counts at SNP sites from the input BAM file",
        oneLineSummary = "Gets alt/ref counts at SNP sites from a BAM file",
        programGroup = SparkProgramGroup.class)
public final class GetHetPulldownSpark extends GATKSparkTool {
    protected static final String MODE_FULL_NAME = "mode";
    protected static final String MODE_SHORT_NAME = "m";

    protected static final String NORMAL_MODE_ARGUMENT = "normal";
    protected static final String TUMOR_MODE_ARGUMENT = "tumor";

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public boolean requiresReference() { return true; }

    @Argument(
            doc = "Mode for filtering SNP sites (normal mode performs binomial test to check for heterozygosity).",
            fullName = MODE_FULL_NAME,
            shortName = MODE_SHORT_NAME,
            optional = false
    )
    protected String mode;

    @Argument(
            doc = "Interval-list file of common SNPs (normal-mode only).",
            fullName = GetHetCoverage.SNP_FILE_FULL_NAME,
            shortName = GetHetCoverage.SNP_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpFile;

    @Argument(
        doc = "Normal-mode heterozygous allele fraction.",
        fullName = GetHetCoverage.HET_ALLELE_FRACTION_FULL_NAME,
        shortName = GetHetCoverage.HET_ALLELE_FRACTION_SHORT_NAME,
        optional = false
    )
    protected static double hetAlleleFraction = 0.5;

    @Argument(
            doc = "Normal-mode p-value threshold for binomial test for heterozygous SNPs.",
            fullName = GetHetCoverage.PVALUE_THRESHOLD_FULL_NAME,
            shortName = GetHetCoverage.PVALUE_THRESHOLD_SHORT_NAME,
            optional = false
    )
    protected static double pvalThreshold = 0.05;

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    protected String out;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        if (!(mode.equals(NORMAL_MODE_ARGUMENT) || mode.equals(TUMOR_MODE_ARGUMENT))) {
            throw new UserException.BadArgumentValue(MODE_FULL_NAME, mode, "Mode must be normal or tumor.");
        }
        final TargetCollection<Locatable> snps = getSNPs();

        final IntervalsSkipList<Locatable> isl = new IntervalsSkipList<>(snps.targets());
        final Broadcast<IntervalsSkipList<Locatable>> islBroad = ctx.broadcast(isl);

        final JavaRDD<GATKRead> rawReads = getReads();
        final JavaRDD<GATKRead> reads = rawReads.filter(read -> !read.isUnmapped() && read.getStart() <= read.getEnd() && !read.isDuplicate());
        final JavaPairRDD<GATKRead, ReferenceBases> readsWithRefBases = BroadcastJoinReadsWithRefBases.addBases(getReference(), reads);

        Function<Tuple2<Locatable, Tuple2<Integer, Integer>>, Boolean> filter;
        if (mode.equals(NORMAL_MODE_ARGUMENT)) {
            filter = keyValue -> isRefAltCountHetCompatible(keyValue._2());
        } else {
            filter = keyValue -> isOverCountThresholds(keyValue._2());
        }

        final Map<Locatable, Tuple2<Integer, Integer>> byKey = readsWithRefBases.flatMapToPair(readWithRefBases ->
                islBroad.getValue().getOverlapping(new SimpleInterval(readWithRefBases._1())).stream()
                        .map(over -> new Tuple2<>(over, getRefAltTuple2(readWithRefBases, over)))
                        .collect(Collectors.toList()))
                .reduceByKey((t1, t2) -> new Tuple2<>(t1._1() + t2._1(), t1._2() + t2._2()))
                .filter(filter)
                .collectAsMap();

        final SortedMap<Locatable, Tuple2<Integer, Integer>> byKeySorted = new TreeMap<>(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
        byKeySorted.putAll(byKey);

        print(byKeySorted, System.out);

        if (out != null){
            final File file = new File(out);
            try(final OutputStream outputStream = BucketUtils.createFile(file.getPath(), (PipelineOptions) null);
                final PrintStream ps = new PrintStream(outputStream)) {
                print(byKeySorted, ps);
            } catch(final IOException e){
                throw new UserException.CouldNotCreateOutputFile(file, e);
            }
        }
    }

    private Tuple2<Integer, Integer> getRefAltTuple2(final Tuple2<GATKRead, ReferenceBases> readWithRefBases, final Locatable snpSite) {
        final GATKRead read = readWithRefBases._1();
        final ReferenceBases referenceBases = readWithRefBases._2();
        final int offset = snpSite.getStart() - read.getStart();
        final byte readBase = read.getBases()[offset];
        final byte refBase = referenceBases.getBases()[offset];
        if (readBase == refBase) {
            return new Tuple2<>(1, 0);
        }
        return new Tuple2<>(0, 1);
    }

    private void print(final SortedMap<Locatable, Tuple2<Integer, Integer>> byKeySorted, final PrintStream ps) {
        ps.println(AllelicCountCollection.CONTIG_COLUMN_NAME + "\t" + AllelicCountCollection.POSITION_COLUMN_NAME
                + "\t" + AllelicCountCollection.REF_COUNT_COLUMN_NAME + "\t" + AllelicCountCollection.ALT_COUNT_COLUMN_NAME);
        for (final Locatable loc : byKeySorted.keySet()){
            ps.println(loc.getContig() + "\t" + loc.getStart() + "\t" + byKeySorted.get(loc)._1() + "\t" + byKeySorted.get(loc)._2());
        }
    }

    private TargetCollection<Locatable> getSNPs() {
        final IntervalList snps = IntervalList.fromFile(snpFile);
        //NOTE: java does not allow conversion of List<Interval> to List<Locatable>.
        //So we do this trick
        return new HashedListTargetCollection<>(snps.getIntervals().stream().collect(Collectors.toList()));
    }

    private static boolean isRefAltCountHetCompatible(final Tuple2<Integer, Integer> refAltCounts) {
        final int refCount = refAltCounts._1();
        final int altCount = refAltCounts._2();
        final int totalCount = refCount + altCount;
        final int majorReadCount = Math.max(refCount, altCount);
        final int minorReadCount = Math.min(refCount, altCount);

        if (majorReadCount == 0 || minorReadCount == 0 || totalCount <= HetPulldownCalculator.READ_DEPTH_THRESHOLD) {
            return false;
        }

        final double pval = new BinomialTest().binomialTest(totalCount, majorReadCount, hetAlleleFraction,
                AlternativeHypothesis.TWO_SIDED);

        return pval >= pvalThreshold;
    }

    private static boolean isOverCountThresholds(final Tuple2<Integer, Integer> refAltCounts) {
        final int refCount = refAltCounts._1();
        final int altCount = refAltCounts._2();
        final int totalCount = refCount + altCount;

        return refCount > 0 && altCount > 0 && totalCount > HetPulldownCalculator.READ_DEPTH_THRESHOLD;
    }
}

