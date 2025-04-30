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
import numpy as np

@stencil_collection("copy")
def copy(
    vx_tmp: Field[IJK, np.float32],
    vy_tmp: Field[IJK, np.float32],
    vx: Field[IJK, np.float32],
    vy: Field[IJK, np.float32]
):

    with computation(PARALLEL), interval(...):
        vx_tmp = vx
        vy_tmp = vy

@stencil_collection("backup")
def backup(
    tracer: Field[IJK, np.float32],
    tracer_e: Field[IJK, np.float32]
):
    with computation(PARALLEL), interval(...):
        tracer = tracer_e