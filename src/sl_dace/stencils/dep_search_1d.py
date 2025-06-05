import dace
from sl_dace.utils.dims import I, J, K
from sl_dace.utils.typingx import dtype_float, dtype_int

import numpy as np

@dace.program
def dep_search_1d(
    vx_e: dtype_float[I, J, K],
    vx_tmp: dtype_float[I, J, K],
    i_a: dtype_int[I, J, K],
    i_d: dtype_int[I, J, K],
    lx: dtype_float[I, J, K],
    dx: dtype_float,
    dth: dtype_float,
):
    # Temporary
    cx = np.ndarray([I, J, K], dtype=dtype_float)

    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        # cx is a cfl
        cx[i, j, k] = -dth * (vx_e[i, j, k] + vx_tmp[i, j, k]) / dx  # dth = dt / 2
        i_d[i, j, k]= i_a[i, j, k] + floor(cx[i, j, k])
        lx[i, j, k] = cx[i, j, k] - floor(cx[i, j, k])
