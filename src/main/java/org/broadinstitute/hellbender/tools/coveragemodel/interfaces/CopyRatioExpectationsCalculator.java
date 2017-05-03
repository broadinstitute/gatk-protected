package org.broadinstitute.hellbender.tools.coveragemodel.interfaces;

import org.broadinstitute.hellbender.tools.coveragemodel.CopyRatioCallingMetadata;
import org.broadinstitute.hellbender.tools.coveragemodel.CopyRatioExpectations;
import org.broadinstitute.hellbender.tools.coveragemodel.CopyRatioHiddenMarkovModelResults;
import org.broadinstitute.hellbender.tools.exome.Target;

import javax.annotation.Nonnull;
import java.util.List;

/**
 * An interface for calculating copy ratio prior and posterior expectations
 *
 * @param <D> emission data type
 * @param <S> hidden state type
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public interface CopyRatioExpectationsCalculator<D, S> {
    /**
     * Calculates copy ratio posterior expectations for a given list of targets and emission data
     *
     * @param copyRatioCallingMetadata sample metadata (sex genotype, depth of coverage, ...)
     * @param targetsList list of targets corresponding to the list of emission data
     * @param emissionData list of emission probability calculation data
     * @return posterior expectations
     */
    CopyRatioExpectations getCopyRatioPosteriorExpectations(@Nonnull final CopyRatioCallingMetadata copyRatioCallingMetadata,
                                                            @Nonnull final List<Target> targetsList,
                                                            @Nonnull final List<D> emissionData);

    /**
     * Calculates copy ratio prior expectations for a given list of targets and emission data
     *
     * @param copyRatioCallingMetadata sample metadata (sex genotype, depth of coverage, ...)
     * @param targetsList list of targets corresponding to the list of emission data
     * @return prior expectations
     */
    CopyRatioExpectations getCopyRatioPriorExpectations(@Nonnull final CopyRatioCallingMetadata copyRatioCallingMetadata,
                                                        @Nonnull final List<Target> targetsList);

    /**
     * Performs forward-backward and Viterbi algorithm
     *
     * @param copyRatioCallingMetadata sample metadata (sex genotype, depth of coverage, ...)
     * @param targetsList list of targets corresponding to the list of emission data
     * @param emissionData list of emission probability calculation data
     * @return an instance of {@link CopyRatioHiddenMarkovModelResults}
     */
    CopyRatioHiddenMarkovModelResults<D, S> getCopyRatioHiddenMarkovModelResults(@Nonnull final CopyRatioCallingMetadata copyRatioCallingMetadata,
                                                                                 @Nonnull final List<Target> targetsList,
                                                                                 @Nonnull final List<D> emissionData);

    void initializeCaches(@Nonnull final List<Target> allTargets);

    void clearCaches();
}