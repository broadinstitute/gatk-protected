package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link ImmutableComputableGraphUtils}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class ImmutableComputableGraphUtilsUnitTest extends BaseTest {

    @Test(expectedExceptions = ImmutableComputableGraphUtils.ImmutableComputableGraphBuilder.DuplicateNodeKeyException.class)
    public void testDuplicatePrimitiveNode() {
        ImmutableComputableGraph.builder()
                .primitiveNode("x", new String[] {}, new DuplicableNDArray())
                .primitiveNode("x", new String[] {}, new DuplicableNDArray())
                .build();
    }

    @Test(expectedExceptions = ImmutableComputableGraphUtils.ImmutableComputableGraphBuilder.DuplicateNodeKeyException.class)
    public void testDuplicateComputableNode() {
        ImmutableComputableGraph.builder()
                .computableNode("x", new String[] {}, new String[] {}, null, true)
                .computableNode("x", new String[] {}, new String[] {}, null, true)
                .build();
    }
}
