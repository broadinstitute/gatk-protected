package org.broadinstitute.hellbender.utils.segmenter;

import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;
import java.io.File;

/**
 * Calls an R script to perform segmentation
 */
public final class RCBSSegmenter {
    private static final String R_SCRIPT = "CBS.R";

    private RCBSSegmenter() {
    }

    public static File segment(String sample_name, String tn_file, String output_file, Float min_log_value) {
        final RScriptExecutor exectutor = new RScriptExecutor();
        exectutor.addScript(new Resource(R_SCRIPT, RCBSSegmenter.class));
        exectutor.addArgs(sample_name, tn_file, output_file, min_log_value);
        exectutor.exec();
        return new File(output_file);
    }
}