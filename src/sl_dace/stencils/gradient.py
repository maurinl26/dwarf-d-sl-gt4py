# -*- coding: utf-8 -*-
import dace
from sl_dace.utils.typingx import dtype_float
from sl_dace.utils.dims import I, J, K

def c2_x(f: dtype_float[I, J, K], df_dx: dtype_float[I, J, K], dx: dtype_float):

    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        df_dx = (f[-1, 0, 0] - f[-1, 0, 0]) / (2 * dx)


def c2_y(f: dtype_float[I, J, K], df_dy: dtype_float[I, J, K], dy: dtype_float):

    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        df_dy = (f[0, -1, 0] - f[0, 1, 0]) / (2 * dy)
