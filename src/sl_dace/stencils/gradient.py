# -*- coding: utf-8 -*-
import dace
from sl_dace.utils.typingx import dtype_float
from sl_dace.utils.dims import I, J, K

def c2_x(f: dtype_float[I, J, K], df_dx: dtype_float[I, J, K], dx: dtype_float):

    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        df_dx[i, j ,k]  = (f[i-1, j, k] - f[i-1, j, k]) / (2 * dx)


def c2_y(f: dtype_float[I, J, K], df_dy: dtype_float[I, J, K], dy: dtype_float):

    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        df_dy[i, j ,k] = (f[i, j-1, k] - f[i, j+1, k]) / (2 * dy)
