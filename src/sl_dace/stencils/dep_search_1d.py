from ifs_physics_common.framework.stencil import stencil_collection
from gt4py.cartesian.gtscript import (
    stencil,
    Field,
    IJK,
    floor,
    computation,
    PARALLEL,
    interval,
)

@stencil_collection("dep_search_1d")
def dep_search_1d(
    vx_e: Field[IJK, float],
    vx_tmp: Field[IJK, float],
    i_a: Field[IJK, int],
    i_d: Field[IJK, int],
    lx: Field[IJK, float],
    dx: float,
    dth: float,
):
    with computation(PARALLEL), interval(...):
        trajx = -dth * (vx_e + vx_tmp) / dx  # dth = dt / 2
        i_d = i_a + floor(trajx)
        lx = trajx - floor(trajx)
