package org.broadinstitute.hellbender.utils.segmenter;

import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

/**
 * Calls an R script to perform segmentation
 */
public final class RCBSSegmenter {
    private static final String R_SCRIPT = "CBS.R";

    public RCBSSegmenter(String sample_name, String tn_file, String output, Float min_log_value) {
        final RScriptExecutor exectutor = new RScriptExecutor();
        exectutor.addScript(new Resource(R_SCRIPT, RCBSSegmenter.class));
        exectutor.addArgs(sample_name, tn_file, output, min_log_value);
        exectutor.exec();
    }

    public void main(String[] args, String sample_name, String tn_file, String output, Float min_log_value) {
        new RCBSSegmenter(sample_name, tn_file, output, min_log_value);
    }
}