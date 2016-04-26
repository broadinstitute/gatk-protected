package org.broadinstitute.hellbender.tools.exome;

import java.util.stream.Stream;

/**
 * Created by davidben on 11/30/15.
 */
public enum AllelicCountTableColumns {
    CONTIG("CONTIG"), POSITION("POS"), REF_COUNT("REF_COUNT"), ALT_COUNT("ALT_COUNT"), REF_NUCLEOTIDE("REF_NUCLEOTIDE"),
    ALT_NUCLEOTIDE("ALT_NUCLEOTIDE"), READ_DEPTH("READ_DEPTH"), HET_LOG_ODDS("HET_LOG_ODDS");

    private String columnName;

    AllelicCountTableColumns(String columnName) {this.columnName = columnName; }

    @Override
    public String toString() {
        return columnName;
    }

    public static final String[] BASIC_COLUMN_NAME_ARRAY = {
            CONTIG.toString(),
            POSITION.toString(),
            REF_COUNT.toString(),
            ALT_COUNT.toString()};

    public static final String[] INTERMEDIATE_COLUMN_NAME_ARRAY = {
            CONTIG.toString(),
            POSITION.toString(),
            REF_COUNT.toString(),
            ALT_COUNT.toString(),
            REF_NUCLEOTIDE.toString(),
            ALT_NUCLEOTIDE.toString(),
            READ_DEPTH.toString()
    };

    public static final String[] FULL_COLUMN_NAME_ARRAY = {
            CONTIG.toString(),
            POSITION.toString(),
            REF_COUNT.toString(),
            ALT_COUNT.toString(),
            REF_NUCLEOTIDE.toString(),
            ALT_NUCLEOTIDE.toString(),
            READ_DEPTH.toString(),
            HET_LOG_ODDS.toString()
    };

}