package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.barclay.argparser.Argument;

/**
 * Arguments for native PairHMM implementations
 */
public class PairHMMNativeArgumentCollection {

    @Argument(fullName = "nativePairHmmThreads", shortName = "nativePairHmmThreads", doc="Number of threads native pairHMM implementations should use", optional = true)
    private int pairHmmNativeThreads = 1;

    @Argument(fullName = "nativePairHmmDoublePrecision", shortName = "nativePairHmmDoublePrecision", doc="Use double precision in the native pairHmm. " +
            "This is slower but matches the java implementation better", optional = true)
    private boolean useDoublePrecision = false;

    public PairHMMNativeArguments getPairHMMArgs(){
        final PairHMMNativeArguments args = new PairHMMNativeArguments();
        args.maxNumberOfThreads = pairHmmNativeThreads;
        args.useDoublePrecision = useDoublePrecision;
        return args;
    }

}
