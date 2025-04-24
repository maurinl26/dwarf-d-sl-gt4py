import numpy as np
from gt4py.cartesian.gtscript import (
    computation,
    PARALLEL,
    interval,
    Field,
    IJK
)

def settls_init(
        vx_e: Field[IJK, np.float32],
        vy_e: Field[IJK, np.float32],
        vx: Field[IJK, np.float32],
        vy: Field[IJK, np.float32],
        vx_p: Field[IJK, np.float32],
        vy_p: Field[IJK, np.float32]
):

    with computation(PARALLEL), interval(...):
        vx_e = 2 * vx - vx_p
        vy_e = 2 * vy - vy_p

def nesc_init(
        vx_e: Field[IJK, np.float32],
        vy_e: Field[IJK, np.float32],
        vx: Field[IJK, np.float32],
        vy: Field[IJK, np.float32]
):

    with computation(PARALLEL), interval(...):
        vx_e = vx
        vy_e = vy