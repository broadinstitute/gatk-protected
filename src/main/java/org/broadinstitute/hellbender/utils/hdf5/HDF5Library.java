package org.broadinstitute.hellbender.utils.hdf5;

import ncsa.hdf.hdf5lib.H5;

/**
 * Utils to handle HDF5 files.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class HDF5Library {


    private

    /**
     * Check whether HDF5 is supported given the presence of the necessary native
     * libraries.
     * @return {@code true} iff HDF5 is supported in the current host.
     */
    public static boolean isHDF5Supported() {
        try {

            return H5.H5open() >= 0;
        } catch (final Exception ex) {
            return false;
        }
    }
}
