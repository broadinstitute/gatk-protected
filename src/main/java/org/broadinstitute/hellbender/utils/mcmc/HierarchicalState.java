package org.broadinstitute.hellbender.utils.mcmc;

import java.util.ArrayList;
import java.util.List;

public class HierarchicalState extends ParameterizedState {
    private final List<ParameterizedState> subModelStates;

    public HierarchicalState(final List<Parameter> parameters, final List<ParameterizedState> subModelStates) {
        super(parameters);
        this.subModelStates = new ArrayList<>(subModelStates);
    }
}
