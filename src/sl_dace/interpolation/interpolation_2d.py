import numpy as np
import dace
from sl_dace.dims import I, J, K, H

# to dace sdfg
def interpolate_lin_2d(
    psi: dace.float32[I, J, K],
    lx: dace.float32[I, J, K],
    ly: dace.float32[I, J, K],
    i_dep: dace.int32[I + H, J + H, K],
    j_dep: dace.int32[I + H, J + H, K],
    psi_dep: dace.float32[I + H, J + H, K],
):
    """Perform a 2d linear interpolation on a regular horizontal plane.

    Note: H stands for the size of the halo (corresponding to max CFL)

    Args:
        psi (np.ndarray): field to interpolate
        lx (np.ndarray): fractional shift from departure point on axis x
        ly (np.ndarray): fractional shift from departure point on axis y
        i_dep (np.ndarray): index of departure point on axis x
        j_dep (np.ndarray): index of departure point on axis y
    """

    # Lookup
    for i, j, k in dace.map[0:I, 0:J, 0:K]:

        with dace.tasklet:

            id_lookup << i_dep[i, j, k]
            jd_lookup << j_dep[i, j, k]
            wx << lx[i, j, k]
            wy << ly[i, j, k]

            # 4 points to interpolate
            lower_left << psi[id_lookup, jd_lookup, k]
            lower_right << psi[id_lookup + 1, jd_lookup, k]
            upper_left << psi[id_lookup, jd_lookup + 1, k]
            upper_right << psi[id_lookup + 1, jd_lookup + 1, k]

            interp_value= (1 - wy) * (
                (1 - wx) * lower_left + wx * lower_right
        ) + wy * (
                (1 - wx) * upper_left + wx * upper_right
        )
            interp_value = psi_dep[i, j, k]

