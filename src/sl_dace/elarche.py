import dace
from typing import Tuple
import numpy as np
from sl_dace.interpolation.interpolation_2d import interpolate_lin_2d
from sl_dace.diagnostics import diagnostic_lipschitz

@dace.program
def lagrangian_search(
    dx: dace.float64,
    dy: dace.float64,
    dth: dace.float64,
    bcx_kind: dace.int32,
    bcy_kind: dace.int32,
    I: np.ndarray,
    J: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    nitmp: dace.int32 = 4,
) -> Tuple[np.ndarray]:
    """Research departure point for a given grid and velocity field.
    Terminates on nsiter iterations.

    Args:
        x (np.ndarray): grid of arrival points
        v (np.ndarray): velocity fields
        nsiter (int, optional): number of iterations. Defaults to 10.

    Returns:
        np.ndarray: departure point
    """

    #todo: move to gt4py stencil
    vx_tmp = vx.copy()
    vy_tmp = vy.copy()

    # Array declaration
    for l in range(nitmp):
        lx, i_d = dep_search_1d(I, vx_e, vx_tmp, dx, dth)
        ly, j_d = dep_search_1d(J, vy_e, vy_tmp, dy, dth)

        # todo: move to gt4py stencil
        lipschitz = diagnostic_lipschitz(
            vx_tmp, vy_tmp, dx, dy, dth
        )

        ####### Interpolation for fields ########
        vx_tmp = interpolate_lin_2d(
            vx, lx, ly, i_d, j_d, bcx_kind, bcy_kind, I, J
        )

        vy_tmp = interpolate_lin_2d(
            vy, lx, ly, i_d, j_d, bcx_kind, bcy_kind, I, J
        )

    return lx, ly, i_d, j_d