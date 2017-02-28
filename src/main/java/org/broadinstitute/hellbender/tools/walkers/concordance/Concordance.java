package org.broadinstitute.hellbender.tools.walkers.concordance;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Optional;

@CommandLineProgramProperties(
        summary = "Evaluate a vcf against a vcf of validated (true) variants",
        oneLineSummary = "Evaluate a vcf against a vcf of validated (true) variants",
        programGroup = VariantProgramGroup.class
)

/**
 *
 * This tool evaluates a vcf against a validated truth vcf. We assume that the truth vcf only contains PASS variants.
 * The summary statistics (# true positives, # false positives, # false negatives, sensitivity, precision)
 * are reported as a summary tsv (--summary). Also reported as a tsv (--table) are evaluated variants of interest
 * (i.e. true positives, false positives, false negatives) and their annotations (ref allele, alt allele, minor allele fraction, etc.)
 * for further statistical analysis e.g. with an R script
 *
 * java -jar gatk.jar Concordance -V na12878-eval.vcf --truth na12878-truth.vcf --table table.tsv --summary summary.tsv --evalSampleName NA12878
 *
 * Created by tsato on 1/30/17.
 */
public class Concordance extends VariantWalker {
    @Argument(doc = "truth vcf (tool assumes all sites in truth are PASS)", fullName= "truth", shortName = "T")
    protected File truth;

    @Argument(doc = "TO BE IMPLEMENTED", fullName= "confidence", shortName = "C", optional = true)
    protected File highConfidenceRegion;

    @Argument(doc = "A table of evaluation variants of interest and their basic annotations",
            fullName="table", shortName = "O")
    protected File table;

    @Argument(doc = "A table of summary statistics (true positives, sensitivity, etc.)",
            fullName="summary", shortName = "S")
    protected File summary;

    @Argument(doc = "sample name of the eval variants", fullName ="evalSampleName", shortName = "sample")
    protected String evalSampleName;

    // TODO: output a vcf of false positives (and negative) sites?
    // TODO: take a low confidence region where a false positive there is not counted (masked in dream?)
    // TODO: what about other variant types besides SNP and INDEL?

    // we count evaluation status in these variables
    private long snpTruePositives = 0;
    private long snpFalsePositives = 0;
    private long snpFalseNegatives = 0;
    private long indelTruePositives = 0;
    private long indelFalsePositives = 0;
    private long indelFalseNegatives = 0;

    private VariantStatusRecord.Writer variantStatusWriter;

    private Iterator<VariantContext> truthIterator;
    private Optional<VariantContext> currentTruthVariant;

    private boolean exhaustedTruthVariants = false;

    private VariantContextComparator variantContextComparator;

    // Possible statuses
    private static String FALSE_NEGATIVE = "FALSE_NEGATIVE";
    private static String TRUE_POSITIVE = "TRUE_POSITIVE";
    private static String FALSE_POSITIVE = "FALSE_POSITIVE";

    @Override
    public void onTraversalStart() {
        Utils.regularReadableUserFile(truth);

        variantStatusWriter = VariantStatusRecord.getWriter(table);

        SAMSequenceDictionary dictionary = getSequenceDictionaryForDrivingVariants();
        variantContextComparator = new VariantContextComparator(dictionary);

        if (! getHeaderForVariants().getSampleNamesInOrder().contains(evalSampleName)){
            throw new IllegalArgumentException(String.format("the eval sample %s does not exist in vcf", evalSampleName));
        }

        final FeatureDataSource<VariantContext> truthSource = new FeatureDataSource<>(truth);
        truthIterator = truthSource.iterator();
        currentTruthVariant = getNextNonSVTruthVariant();
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ) {
        // TODO: When we give the tool an interval List (-L), variant walker subsets the range of driving variants to these intervals
        // for you behind the scene. But it doesn't do the same for the truth vcf i.e. truth vcf starts starts at chr 1 when -L 21.

        // case 0: We've reached the end of truth variants; eval variant is a FP unless filtered
        // note exhaustedTruthVariants is not always equal to !currentTruthVariant.isPresent();
        // in other words, iterator.hasNext() can give you false
        // when currentTruthVariant is non-empty (the last truth variant is either SNP or INDEL)
        if (exhaustedTruthVariants){
            if (variant.isNotFiltered()){
                if (variant.isSNP()) { snpFalsePositives++; } else { indelFalsePositives++; }
                writeThisRecord(new VariantStatusRecord(variant, FALSE_POSITIVE, evalSampleName), variant);
            }

            return;
        }

        // case 1: Eval leapfrogged truth; we missed at least one variant.
        // Move the truth cursor forward until it catches up to eval.
        // Note we must check for this first, *then* check for the other two conditions;
        // otherwise we move the eval forward without giving it a diagnosis
        while (variantContextComparator.compare(variant, currentTruthVariant.get()) > 0) {
            if (currentTruthVariant.get().isSNP()) { snpFalseNegatives++; } else { indelFalseNegatives++; }
            if (truthIterator.hasNext()) {
                currentTruthVariant = getNextNonSVTruthVariant();
                if (! currentTruthVariant.isPresent()) {
                    return;
                }
            } else {
                exhaustedTruthVariants = true;
                return;
            }
        }

        // case 2: the position of the two variants match
        if (variantContextComparator.compare(variant, currentTruthVariant.get()) == 0) {
            // do the genotypes match too?
            if (variant.isFiltered()) {
                if (variant.isSNP()) { snpFalseNegatives++; } else { indelFalseNegatives++; }
                writeThisRecord(new VariantStatusRecord(variant, FALSE_NEGATIVE, evalSampleName), variant);
            } else if (allelesMatch(variant, currentTruthVariant.get())) {
                if (variant.isSNP()) { snpTruePositives++; } else { indelTruePositives++; }
                writeThisRecord(new VariantStatusRecord(variant, TRUE_POSITIVE, evalSampleName), variant);
            } else {
                // the eval variant matches truth's position but their alleles don't match
                // e.g. genotype or alleles don't match
                if (variant.isSNP()) {
                    snpFalsePositives++;
                    snpFalseNegatives++;
                } else {
                    indelFalsePositives++;
                    indelFalseNegatives++;
                }
                writeThisRecord(new VariantStatusRecord(variant, FALSE_POSITIVE, evalSampleName), variant);
            }

            if (truthIterator.hasNext()) {
                currentTruthVariant = getNextNonSVTruthVariant();
            } else {
                exhaustedTruthVariants = true;
            }

            return;
        }

        // case 3: truth got ahead of eval; the eval variant must be a false positive if not filtered
        if (variantContextComparator.compare(variant, currentTruthVariant.get()) < 0) {
            if (variant.isNotFiltered()) {
                if (variant.isSNP()) { snpFalsePositives++; } else { indelFalsePositives++; }
                writeThisRecord(new VariantStatusRecord(variant, FALSE_POSITIVE, evalSampleName), variant);
            }

            // we don't care about true negatives - don't record them in the table
            return;
        }
    }

    // ugly code alert
    // this method just moves the iterator forward but also ensures that it skips SV's
    // assumes that the caller already checked the truth iterator has next
    // TODO: extend the truth iterator such that the iterator.next() automatically skips SVs
    private Optional<VariantContext> getNextNonSVTruthVariant(){
        VariantContext nextVariant;
        do {
            // if the truth vcf ends with an SV, truthIterator.next() may return null
            // so make sure that we didn't reach the end of truth
            if (! truthIterator.hasNext()){
                exhaustedTruthVariants = true;
                return Optional.empty();
            }
            nextVariant = truthIterator.next();
        } while (nextVariant.getStructuralVariantType() != null);

        return Optional.of(nextVariant);
    }

    private void writeThisRecord(final VariantStatusRecord record, final VariantContext variant){
        try {
            variantStatusWriter.writeRecord(record);
        } catch (Exception e) {
            throw new UserException(String.format("Encountered an IO exception writing a record at chrom %d, pos ",
                    variant.getContig(), variant.getStart()), e);
        }
    }

    @Override
    public Object onTraversalSuccess(){
        if (! exhaustedTruthVariants){
            // the remaining variants in truth are false negatives
            if (currentTruthVariant.get().isSNP()) { snpFalseNegatives++; } else { indelFalseNegatives++; }
            while (truthIterator.hasNext()){
                currentTruthVariant = getNextNonSVTruthVariant();
                if (! currentTruthVariant.isPresent()) { break; }
                if (currentTruthVariant.get().isSNP()) { snpFalseNegatives++; } else { indelFalseNegatives++; }
            }
        }

        try ( ConcordanceSummaryRecord.Writer concordanceSummaryWriter = ConcordanceSummaryRecord.getWriter(summary) ){
            // TODO: write all the summary
            concordanceSummaryWriter.writeRecord(new ConcordanceSummaryRecord(VariantContext.Type.SNP, snpTruePositives, snpFalsePositives, snpFalseNegatives));
            concordanceSummaryWriter.writeRecord(new ConcordanceSummaryRecord(VariantContext.Type.INDEL, indelTruePositives, indelFalsePositives, indelFalseNegatives));

            // We must close the variant status writer to flush the buffer.
            // Summary writer will be closed implicitly
            variantStatusWriter.close();
        } catch (IOException e){
            throw new UserException("Encountered an IO exception writing the concordance summary table", e);
        }

        return "SUCCESS";
    }

    // TODO: eventually this should be an abstract method to be overridden by a subclass
    // such that the user can adjust the stringency of allele comparison
    protected static boolean allelesMatch(final VariantContext eval, final VariantContext truth){
        final boolean sameContig = truth.getContig().equals(eval.getContig());
        final boolean sameStartPosition = truth.getStart() == eval.getStart();
        final boolean sameRefAllele = truth.getReference().equals(eval.getReference());
        // we assume that the truth has a single alternate allele
        final boolean containsAltAllele = eval.getAlternateAlleles().contains(truth.getAlternateAllele(0));

        return sameContig && sameStartPosition && sameRefAllele && containsAltAllele;
    }


}
