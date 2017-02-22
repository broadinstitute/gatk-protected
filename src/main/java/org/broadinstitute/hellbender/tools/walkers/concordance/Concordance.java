package org.broadinstitute.hellbender.tools.walkers.concordance;


import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

@CommandLineProgramProperties(
        summary = "Count ",
        oneLineSummary = "Count PASS variants",
        programGroup = VariantProgramGroup.class
)

/**
 * Created by tsato on 1/30/17.
 */
// TODO: maybe this tool should be a variant walker. Or a multiple variant walker.
public class Concordance extends VariantWalker {
    @Argument(doc = "truth vcf (tool assumes all sites in truth are PASS)", fullName= "truth", shortName = "T", optional = false)
    protected File truth;

    // TODO: TO BE IMPLEMENTED
    @Argument(doc = "???", fullName= "confidence", shortName = "C", optional = true)
    protected File confidence_region;

    @Argument(doc = "A table of variants. Prints out for each variant (row) its basic annotations, tumor alt allele fraction, and truth status", fullName="table", shortName = "O", optional = false)
    protected File table;

    @Argument(doc = "A table of summary statistics (true positives, sensitivity, etc.)", fullName="summary", shortName = "S", optional = false)
    protected File summary;

    @Argument(doc = "sample name of the tumor", fullName ="tumorSampleName", shortName = "tumor", optional = false)
    protected String tumorSampleName;

    // TODO: output a vcf of false positives (and negative) sites?
    // TODO: take a low confidence region where a false positive there is not counted (masked in dream?)
    // TODO: ideas. for stratifying, output a table of annotaitons (AF, DP, SNP/INDEL, and truth status)
    // @assumes that all variants in the truth vcf are REAL
    // TODO: make allelesMatch function extendable, such that a user inherits the class to write his/her own evaluator (but get the iteration working first)

    // use MutableLong?
    private long truePositives = 0;
    private long falsePositives = 0;
    private long falseNegatives = 0;

    private VariantStatusRecord.Writer variantStatusWriter;
    private ConcordanceSummaryRecord.Writer concordanceSummaryWriter;

    // TODO: want to make this final but can't
    private Iterator<VariantContext> truthIterator;
    private VariantContext currentTruthVariant;

    private boolean exhaustedTruthVariants = false;

    @Override
    public void onTraversalStart() {
        // TODO: do we still need tumor sample name?
        // TODO: once contamination is rebased do this instead:
        // final String tumorSample = getHeaderForVariants().getMetaDataLine(Mutect2Engine.TUMOR_SAMPLE_KEY_IN_VCF_HEADER).getValue();
        // and get rid of tumor sample name argument
        Utils.regularReadableUserFile(truth);

        variantStatusWriter = VariantStatusRecord.getWriter(table);
        concordanceSummaryWriter = ConcordanceSummaryRecord.getWriter(summary);

        if (! getHeaderForVariants().getSampleNamesInOrder().contains(tumorSampleName)){
            // TODO: do I need to close files?
            throw new IllegalArgumentException(String.format("the tumor sample %s does not exit in vcf", tumorSampleName));
        }

        final FeatureDataSource<VariantContext> truthSource = new FeatureDataSource<>(truth);
        truthIterator = truthSource.iterator();
        currentTruthVariant = truthIterator.next();
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ) {
        // no more truth variants; eval variant is a FP unless filtered
        if (exhaustedTruthVariants && variant.isNotFiltered()){
            falsePositives++;
            return;
        }

        // TODO: use a comparator?
        // TODO: need more asserts to ensure that we aren't breaking invariants
        // case 1: the position of the two variants match
        if (comparePositions(variant, currentTruthVariant) == 0) {
            // do the genotypes match too?
            if (variant.isFiltered()) {
                falseNegatives++;
                writeThisRecord(new VariantStatusRecord(variant, FALSE_NEGATIVE, tumorSampleName), variant);
            } else if (allelesMatch(variant, currentTruthVariant)) {
                truePositives++;
                writeThisRecord(new VariantStatusRecord(variant, TRUE_POSITIVE, tumorSampleName), variant);
            } else {
                // the eval variant matches truth's position but their alleles don't match
                // e.g. genotype or alleles don't match
                falsePositives++;
                falseNegatives++;
                writeThisRecord(new VariantStatusRecord(variant, FALSE_POSITIVE_AND_FALSE_NEGATIVE, tumorSampleName), variant);
            }

            if (truthIterator.hasNext()) {
                currentTruthVariant = truthIterator.next();
            } else {
                exhaustedTruthVariants = true;
            }

            return;
        }

        // case 2: truth got ahead of eval; the eval variant must be a false positive
        if (comparePositions(variant, currentTruthVariant) < 0) {
            if (variant.isNotFiltered()) {
                falsePositives++;
                writeThisRecord(new VariantStatusRecord(variant, FALSE_POSITIVE, tumorSampleName), variant);
            }
            return;
        }

        // case 3: eval leapfrogged truth; we missed at least one variant
        // move the truth cursor forward until it catches up to eval
        while (comparePositions(variant, currentTruthVariant) > 0) {
            falseNegatives++;
            if (truthIterator.hasNext()) {
                currentTruthVariant = truthIterator.next();
            } else {
                exhaustedTruthVariants = true;
                return;
            }
        }
    }

    private void writeThisRecord(final VariantStatusRecord record, final VariantContext variant){
        try {
            variantStatusWriter.writeRecord(record);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception writing a record at chrom %d, pos ",
                    variant.getContig(), variant.getStart()), e);
        }
    }

    @Override
    public Object onTraversalSuccess(){
        if (! exhaustedTruthVariants){
            // the remaining variants in truth are false negatives
            falseNegatives++;
            while (truthIterator.hasNext()){
                falseNegatives++;
            }
        }

        try {
            concordanceSummaryWriter.writeRecord(new ConcordanceSummaryRecord(truePositives, falsePositives, falseNegatives));
        } catch (IOException e){
            throw new UserException("Encountered an IO exception writing the concordance summary table", e);
        }

        return "SUCCESS";
    }

    // Possible statuses
    private static String FALSE_NEGATIVE = "FALSE_NEGATIVE";
    private static String TRUE_POSITIVE = "TRUE_POSITIVE";
    private static String FALSE_POSITIVE_AND_FALSE_NEGATIVE = "FALSE_POSITIVE_AND_FALSE_NEGATIVE";
    private static String FALSE_POSITIVE = "FALSE_POSITIVE";

    // TODO: make variantcontext comparable?
    protected static int comparePositions(final VariantContext truth, final VariantContext eval){
        if (! truth.getContig().equals(eval.getContig()) ) {
            // TODO: make sure x, y, and mt work as expected
            return truth.getContig().compareTo(eval.getContig());
        } else {
            // same chromosome, compare the start positions
            if (truth.getStart() > eval.getStart()) {
                return 1;
            } else if (truth.getStart() < eval.getStart()){
                return -1;
            } else {
                return 0;
            }
        }
    }

    // TODO: eventually this should be an abstract method to be overridden by a subclass
    // TODO: make protected to support subclassing
    protected static boolean allelesMatch(final VariantContext truth, final VariantContext eval){
        final boolean sameContig = truth.getContig().equals(eval.getContig());
        final boolean sameStartPosition = truth.getStart() == eval.getStart();
        final boolean sameRefAllele = truth.getReference().equals(eval.getReference());
        // TODO: assume single alt allele in truth?
        final boolean containsAltAllele = eval.getAlternateAlleles().contains(truth.getAlternateAllele(0));
        return sameContig && sameStartPosition && sameRefAllele && containsAltAllele;
    }


}
