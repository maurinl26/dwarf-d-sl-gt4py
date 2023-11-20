import numba

from sl_python.interpolation import interpolate_cub_2d, interpolate_lin_2d

numba_interpolate_lin_2d = numba.jit(interpolate_lin_2d)
numba_interpolate_cub_2d = numba.jit(interpolate_cub_2d)