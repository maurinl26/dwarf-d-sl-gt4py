import numpy as np
import dace
from typing import Tuple
from config import Config

from sl_dace.utils.sdfg import build_sdfg
from sl_dace.utils.typingx import dtype_float, dtype_int
from sl_dace.stencils.gradient import c2_x, c2_y


class LipschitzDiag:

    def __init__(self,
                 grid: Tuple[int],
                 config: Config,
                 ):

        # todo: a grid class which takes
        self.config = config
        self.grid = grid

        self.c2_x = build_sdfg(c2_x)
        self.c2_y = build_sdfg(c2_y)

    def __call__(self,
                 u: np.ndarray,
                 v: np.ndarray,
                 du_dx: np.ndarray,
                 du_dy: np.ndarray,
                 dv_dx: np.ndarray,
                 dv_dy: np.ndarray,
                 dx: dtype_float,
                 dy: dtype_float,
                 dth: dtype_float
                 ):

        self.c2_x(u, du_dx, dx, origin=(1, 0, 0), domain=self.grid)
        self.c2_x(v, dv_dx, dx, origin=(1, 0, 0), domain=self.grid)

        self.c2_y(u, du_dy, dy, origin=(0, 1, 0), domain=self.grid)
        self.c2_y(v, dv_dy, dy, origin=(0, 1, 0), domain=self.grid)

        # dace.reduce for performances
        return dth * np.maximum(np.maximum(du_dx, du_dy), np.maximum(dv_dx, dv_dy))

