package org.broadinstitute.hellbender.tools.walkers.concordance;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

/**
 * Created by tsato on 2/8/17.
 */

public class ConcordanceSummaryRecord {
    private static final String TRUE_POSITIVE_COLUMN_NAME = "true-positive";
    private static final String FALSE_POSITIVE_COLUMN_NAME = "false-positive";
    private static final String FALSE_NEGATIVE_COLUMN_NAME = "false-negative";
    private static final String SENSITIVITY_COLUMN_NAME = "sensitivity";
    private static final String PRECISION_COLUMN_NAME = "precision";
    private static final String[] SUMMARY_TABLE_COLUMN_HEADER =
            {TRUE_POSITIVE_COLUMN_NAME, FALSE_POSITIVE_COLUMN_NAME, FALSE_NEGATIVE_COLUMN_NAME, SENSITIVITY_COLUMN_NAME, PRECISION_COLUMN_NAME};

    final long truePositives;
    final long falsePositives;
    final long falseNegatives;

    public ConcordanceSummaryRecord(final long truePositives, final long falsePositives, final long falseNegatives){
        this.truePositives = truePositives;
        this.falsePositives = falsePositives;
        this.falseNegatives = falseNegatives;
    }

    public long getTruePositives(){
        return truePositives;
    }

    public long getFalsePositives() {
        return falsePositives;
    }

    public long getFalseNegatives(){
        return falseNegatives;
    }

    public double getSensitivity(){
        return (double) truePositives / (truePositives + falseNegatives);
    }

    public double getPrecision(){
        return (double) truePositives / (truePositives + falsePositives);
    }


    private static class ConcordanceSummaryTableWriter extends TableWriter<ConcordanceSummaryRecord> {
        private ConcordanceSummaryTableWriter(final File output) throws IOException {
            super(output, new TableColumnCollection(ConcordanceSummaryRecord.SUMMARY_TABLE_COLUMN_HEADER));
        }

        @Override
        protected void composeLine(final ConcordanceSummaryRecord record, final DataLine dataLine) {
            dataLine.set(ConcordanceSummaryRecord.TRUE_POSITIVE_COLUMN_NAME, record.getTruePositives())
                    .set(ConcordanceSummaryRecord.FALSE_POSITIVE_COLUMN_NAME, record.getFalsePositives())
                    .set(ConcordanceSummaryRecord.FALSE_NEGATIVE_COLUMN_NAME, record.getFalseNegatives())
                    .set(ConcordanceSummaryRecord.SENSITIVITY_COLUMN_NAME, record.getSensitivity())
                    .set(ConcordanceSummaryRecord.PRECISION_COLUMN_NAME, record.getPrecision());
        }
    }

    public static ConcordanceSummaryTableWriter getConcordanceSummaryTableWriter(final File outputTable){
        try (ConcordanceSummaryTableWriter writer = new ConcordanceSummaryTableWriter(outputTable)) {
            return writer;
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", outputTable));
        }
    }
}

