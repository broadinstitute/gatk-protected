package org.broadinstitute.hellbender.utils.mcmc.coordinates;

import org.apache.commons.math3.analysis.function.Logit;
import org.apache.commons.math3.analysis.function.Sigmoid;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * TODO
 */
public final class CoordinateUtils {
    private static final double EPSILON = 1E-10;

    private CoordinateUtils() {}

    public static double transformWalkerCoordinateToBoundedVariable(final double walkerCoordinate,
                                                                    final double min,
                                                                    final double max) {
        Utils.validateArg(min < max, "Minimum bound must be strictly less than maximum bound.");
        return min + (max - min) * new Sigmoid().value(walkerCoordinate);
    }

    public static double transformBoundedVariableToWalkerCoordinate(final double boundedVariable,
                                                                    final double min,
                                                                    final double max) {
        Utils.validateArg(min < max, "Minimum bound must be strictly less than maximum bound.");
        Utils.validateArg(min <= boundedVariable && boundedVariable <= max, "Variable value " + boundedVariable + " is not within variable bounds [" + min + ", " + max + "].");
        return new Logit().value((boundedVariable - min) / (max - min));
    }

    public static double calculateLogJacobianFactor(final double boundedVariable,
                                                    final double min,
                                                    final double max) {
        Utils.validateArg(min < max, "Minimum bound must be strictly less than maximum bound.");
        return boundedVariable < min || boundedVariable > max ?
                Double.NEGATIVE_INFINITY :
                Math.log(boundedVariable - min) + Math.log(max - boundedVariable) - Math.log(max - min);
    }
}
