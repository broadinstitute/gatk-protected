package org.broadinstitute.hellbender.utils.segmenter;

import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

/**
 * Calls and Rscript to perform segmentation
 */
public class RCBSSegmenter extends RScriptExecutor {
    private static final String R_SCRIPT = "CBS.R";

    public RCBSSegmenter(String sample_name, String tn_file, String output) {
        addScript(new Resource(R_SCRIPT, RCBSSegmenter.class));
        addArgs(sample_name, tn_file, output);
    }

    public void main(String[] args, String sample_name, String tn_file, String output) {
        new RCBSSegmenter(sample_name, tn_file, output);
    }
}