package org.broadinstitute.hellbender.utils.gc;

import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.util.List;

public class GcCorrector {

    public static ReadCountCollection correctGc(final ReadCountCollection readCountCollection, final File fasta){
        final double[] readCounts = readCountCollection.counts().getColumn(0);
        final List<Target> targetList = readCountCollection.targets();
        final double[] gc = new double[targetList.size()];
        SimpleInterval interval;
        for(int i = 0; i < targetList.size(); i++) {
            interval = targetList.get(i).getInterval();

        }

//        TODO: Calculate Curve-fitting function

//        TODO: gcNormalizedReadCounts=ReadCounts/(cff(ReadCounts.gc))
        return readCountCollection;
    }
}
