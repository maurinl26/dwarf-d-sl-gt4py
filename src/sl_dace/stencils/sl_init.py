import dace
from sl_dace.utils.typingx import dtype_float, dtype_int
from sl_dace.utils.dims import I, J, K, KH

def settls_init(
        vx_e: dtype_float[I, J, K],
        vy_e: dtype_float[I, J, K],
        vx: dtype_float[I, J, K],
        vy: dtype_float[I, J, K],
        vx_p: dtype_float[I, J, K],
        vy_p: dtype_float[I, J, K]
):

    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        vx_e[i, j, k] = 2 * vx[i, j, k] - vx_p[i, j, k]
        vy_e[i, j, k] = 2 * vy[i, j, k] - vy_p[i, j, k]

def nesc_init(
        vx_e: dtype_float[I, J, K],
        vy_e: dtype_float[I, J, K],
        vx: dtype_float[I, J, K],
        vy: dtype_float[I, J, K]
):

    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        vx_e[i, j, k] = vx[i, j, k]
        vy_e[i, j, k] = vy[i, j, k]