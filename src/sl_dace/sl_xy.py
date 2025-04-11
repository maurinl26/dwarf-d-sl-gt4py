from typing import Tuple
import numpy as np
import logging

from config import Config
from sl_dace.diagnostics import diagnostic_lipschitz
from sl_dace.interpolation.interpolation_2d import (
    interpolate_lin_2d,
)
import dace

logging.getLogger(__name__)

# todo: Link with dym symbols from gt4py
I = dace.symbol("I")
J = dace.symbol("J")
K = dace.symbol("K")



# ELARCHE
@dace.program
def sl_xy(
    bcx_kind: dace.int32,
    bcy_kind: dace.int32,
    dx: dace.float64,
    dy: dace.float64,
    dth: dace.float64,
    I: dace.int32[I,J,K],
    J: dace.int32[I,J,K],
    vx: dace.float32[I, J, K],
    vy: dace.float32[I, J, K],
    vx_e: dace.float32[I, J, K],
    vy_e: dace.float32[I, J, K],
    tracer: dace.float32[I, J, K],
    nitmp: dace.int32,
) -> np.ndarray:
    """Performs tracer advection with 2D semi lagrangian.
    1: search for departure point
    2: interpolate tracer field

    Args:
        config (Config): grid configuration
        vx (np.ndarray): velocity on x
        vy (np.ndarray): velocity on y
        vx_e (np.ndarray): velocity on x at t + dt
        vy_e (np.ndarray): velocity on y at t + dt
        tracer (np.ndarray): tracer field
        tracer_e (np.ndarray): tracer at t + dt (ebauche)
        nitmp (int): number of iterations for departure search

    Returns:
        np.ndarray: tracer outline (ebauche) at t + dt
    """

    # Recherche semi lag
    lx_d, ly_d, i_d, j_d = lagrangian_search(
        I=I,
        J=J,
        dx=dx,
        dy=dy,
        dth=dth,
        bcx_kind=bcx_kind,
        bcy_kind=bcy_kind,
        vx_e=vx_e,
        vy_e=vy_e,
        vx=vx,
        vy=vy,
        nitmp=nitmp,
    )

    # Interpolate
    tracer_e = interpolate_lin_2d(
        tracer,
        lx_d,
        ly_d,
        i_d,
        j_d,
        bcx_kind,
        bcy_kind,
        I,
        J,
    )

    # todo: add filter
    return tracer_e

