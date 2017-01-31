package org.broadinstitute.hellbender.tools.walkers.mutect;


import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import scala.collection.immutable.Map;
import scala.collection.immutable.Stream;

import java.io.File;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

@CommandLineProgramProperties(
        summary = "Count ",
        oneLineSummary = "Count PASS variants",
        programGroup = VariantProgramGroup.class

)

/**
 * Created by tsato on 1/30/17.
 */
public class Concordance extends CommandLineProgram {
    @Argument(doc = "", fullName= "output", shortName = "V", optional = false)
    protected File eval;

    @Argument(doc = "vcf file of true variants", fullName= "", shortName = "T", optional = false)
    protected File truth;

    @Argument(doc = "", fullName= "", shortName = "", optional = true)
    protected File confidence_region;


    private final int TRUE_POSITIVE_INDEX = 0;
    private final int FALSE_POSITIVE_INDEX = 1;
    private final int FALSE_NEGATIVE_INDEX = 2;


    @Override
    public Object doWork() {
        // Need to calculate true positive rate
        // False positive rate
        // maybe take in a low confidence region where a false positive there is not counted (masked in dream?)
        // precision
        //
        // @assumes that all variants in the truth vcf are REAL
        // TODO: bams to test are here: /humgen/gsa-hpprojects/dev/tsato/wdl/cromwell-executions/M2DreamChallenge/b840bb65-0cd2-48d1-8880-259dcece8036/call-EvaluateD3/execution
        // TODO: ideas. for stratifying, output a table of annotaitons (AF, DP, SNP/INDEL, and truth status)
        // make match function extendable, such that a user inherits the class to write his/her own evaluator (but get the iteration working first)
        final int[] counts = walkAndCount(eval, truth);
        final int truePositives = counts[TRUE_POSITIVE_INDEX];
        final int falsePositives = counts[FALSE_POSITIVE_INDEX];
        final int falseNegatives = counts[FALSE_NEGATIVE_INDEX];

        // sensitivity aka recall
        final int sensitivity = truePositives / (truePositives + falseNegatives);
        final int precision = truePositives / (truePositives + falsePositives);
        System.out.println(String.format("TP: %d", truePositives));
        System.out.println(String.format("FP: %d", falsePositives));
        System.out.println(String.format("FN: %d", falseNegatives));
        System.out.println(String.format("Sensitivity: %d", sensitivity));
        System.out.println(String.format("Precision: %d", precision));

        return "SUCCESS";
    }

    /**
     *
     * Calculate rough concordance between two vcfs, comparing only the positions, alleles, and the first genotype.
     * Borrows heavily from the same function in HaplotypeCallerIntegrationTest
     *
     * @param eval
     * @param truth
     * @return
     */
    private int[] walkAndCount(final File eval, final File truth ) {
        int truePositives = 0;
        int falsePositives = 0;
        int falseNegatives = 0;

        // TODO: to reduce the memory burden it may make more sense to divide each section by chromosomes

        try (final FeatureDataSource<VariantContext> truthSource = new FeatureDataSource<>(truth);
             final FeatureDataSource<VariantContext> evalSource = new FeatureDataSource<>(eval)) {

            Iterator<VariantContext> truthIterator = truthSource.iterator();
            Iterator<VariantContext> evalIterator = evalSource.iterator();

            // TODO: is this how you start an iterator?
            VariantContext truthvc = truthIterator.next();
            VariantContext evalvc = evalIterator.next();

            // TODO: check the FILTER flag
            while (true) {
                if (comparePositions(truthvc, evalvc) == 0) {
                    if (match(truthvc, evalvc)){
                        truePositives++;
                    } else {
                        falsePositives++;
                        falseNegatives++;
                    }

                    if (truthIterator.hasNext() && evalIterator.hasNext()){
                        truthvc = truthIterator.next();
                        evalvc = evalIterator.next();
                        continue;
                    } else {
                        break;
                    }
                }

                if (comparePositions(truthvc, evalvc) > 0) {
                    // truth is ahead of eval
                    falsePositives++;
                    if (evalIterator.hasNext()){
                        evalvc = evalIterator.next();
                        continue;
                    } else {
                        break;
                    }
                }

                if (comparePositions(truthvc, evalvc) < 0) {
                    // eval leapfrogged truth
                    falseNegatives++;
                    if (truthIterator.hasNext()){
                        truthvc = truthIterator.next();
                        continue;
                    } else {
                        break;
                    }
                }
            }

            if (! truthIterator.hasNext()){
                // exhausted the truth vcf first; the rest of variants in eval are false positives
                while (evalIterator.hasNext()){
                    evalvc = evalIterator.next();
                    if (evalvc.isNotFiltered()){
                        falsePositives++;
                    }
                }
            } else {
                // exhausted the eval vcf first; the rest of
                while (truthIterator.hasNext()){
                    truthvc = truthIterator.next();
                    if (truthvc.isNotFiltered()){
                        falseNegatives++;
                    }
                }
            }
        }

        final int[] counts = new int[3];
        counts[TRUE_POSITIVE_INDEX] = truePositives;
        counts[FALSE_POSITIVE_INDEX] = falsePositives;
        counts[FALSE_NEGATIVE_INDEX] = falseNegatives;
        return counts;
    }

    // TODO: make variantcontext comparable?
    private int comparePositions(final VariantContext vc1, final VariantContext vc2){
        // TODO: make sure x, y, and mt work as expected
        if (! vc1.getContig().equals(vc2.getContig()) ) {
            return vc1.getContig().compareTo(vc2.getContig());
        } else {
            // same chromosome, compare the start positions
            if (vc1.getStart() > vc2.getStart()) {
                return 1;
            } else if (vc1.getStart() < vc2.getStart()){
                return -1;
            } else {
                return 0;
            }
        }
    }

    // TODO: eventually this should be an abstract method to be overridden by a subclass
    private boolean match(final VariantContext vc1, final VariantContext vc2){
        final boolean sameContig = vc1.getContig().equals(vc2.getContig());
        final boolean sameStartPosition = vc1.getStart() == vc2.getStart();
        final boolean sameRefAllele = vc1.getReference().equals(vc2.getReference());
        final boolean containsAltAllele = vc2.getAlternateAlleles().contains(vc1.getAlternateAllele(0));
        // TODO: assume single alt allele in truth?
        return sameContig && sameStartPosition && sameRefAllele && containsAltAllele;
    }
}
