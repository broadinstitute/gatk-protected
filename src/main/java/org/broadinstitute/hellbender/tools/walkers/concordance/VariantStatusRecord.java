package org.broadinstitute.hellbender.tools.walkers.concordance;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.bdgenomics.formats.avro.Variant;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.io.IOException;

/**
 * Created by tsato on 2/8/17.
 */
public class VariantStatusRecord {
    private static final String CHROMOSOME_COLUMN_NAME = "CHROM";
    private static final String START_POSITION_COLUMN_NAME = "START";
    private static final String END_POSITION_COLUMN_NAME = "END";
    private static final String REF_ALLELE_COLUMN_NAME = "REF";
    private static final String ALT_ALLELE_COLUMN_NAME = "ALT";
    private static final String VARIANT_TYPE_COLUMN_NAME = "TYPE";
    private static final String ALLELE_FRACTION_COLUMN_NAME = "TUMOR_ALT_AF";
    private static final String TRUTH_STATUS_COLUMN_NAME = "TRUTH_STATUS";
    private static final String[] VARIANT_TABLE_COLUMN_HEADERS = {CHROMOSOME_COLUMN_NAME, START_POSITION_COLUMN_NAME, END_POSITION_COLUMN_NAME,
            REF_ALLELE_COLUMN_NAME, ALT_ALLELE_COLUMN_NAME, VARIANT_TYPE_COLUMN_NAME, ALLELE_FRACTION_COLUMN_NAME, TRUTH_STATUS_COLUMN_NAME};

    String truthStatus;
    VariantContext variantContext;
    // TODO: should be an array of doubles
    String tumorAlleleFraction;

    public VariantStatusRecord(final VariantContext variantContext, final String truthStatus, final String tumorSampleName){
        this.truthStatus = truthStatus;
        this.variantContext = variantContext;
        this.tumorAlleleFraction = (String) variantContext.getGenotype(tumorSampleName).getAnyAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY);
    }

    public String getTruthStatus() {
        return truthStatus;
    }

    public VariantContext getVariantContext() {
        return variantContext;
    }

    public String getTumorAlleleFraction(){ return tumorAlleleFraction; }

    public static class Writer extends TableWriter<VariantStatusRecord> {
        private Writer(final File output) throws IOException {
            super(output, new TableColumnCollection(VariantStatusRecord.VARIANT_TABLE_COLUMN_HEADERS));
        }

        @Override
        protected void composeLine(final VariantStatusRecord record, final DataLine dataLine) {
            dataLine.set(VariantStatusRecord.CHROMOSOME_COLUMN_NAME, record.getVariantContext().getContig())
                    .set(VariantStatusRecord.START_POSITION_COLUMN_NAME, record.getVariantContext().getStart())
                    .set(VariantStatusRecord.END_POSITION_COLUMN_NAME, record.getVariantContext().getEnd())
                    .set(VariantStatusRecord.REF_ALLELE_COLUMN_NAME, record.getVariantContext().getReference().toString())
                    .set(VariantStatusRecord.ALT_ALLELE_COLUMN_NAME, record.getVariantContext().getAlternateAlleles().toString())
                    .set(VariantStatusRecord.VARIANT_TYPE_COLUMN_NAME, record.getVariantContext().getType().toString())
                    .set(VariantStatusRecord.ALLELE_FRACTION_COLUMN_NAME, record.getTumorAlleleFraction()) // TODO: must be able to retrieve the allele fraction from tumor
                    .set(VariantStatusRecord.TRUTH_STATUS_COLUMN_NAME, record.getTruthStatus());
        }
    }

    public static Writer getWriter(final File outputTable){
        try {
            Writer writer = new Writer(outputTable);
            return writer;
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception opening %s.", outputTable), e);
        }
    }


}