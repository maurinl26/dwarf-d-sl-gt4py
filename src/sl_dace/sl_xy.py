import logging
from functools import cached_property

import numpy as np
from config import Config
from sl_dace.interpolation.interpolation_2d import (
    interpolate_lin_2d,
)
from sl_dace.elarche import Elarche
from sl_dace.dims import I, J, K
import dace

from typing import Tuple

logging.getLogger(__name__)


class SmiLagXY:

    def __init__(self, grid: Tuple[int], nitmp: int = 4):
        self.elarche = Elarche(grid, nitmp)

    @cached_property
    def _temporaries(self):
        return {
            "lx": {"grid": (I, J, K), "dtype": float},
            "ly": {"grid": (I, J, K), "dtype": float},
            "i_dep": {"grid": (I, J, K), "dtype": int},
            "j_dep": {"grid": (I, J, K), "dtype": int}
        }

    def __call__(self,
                 state: dict,
                 grid_spacings: Tuple[float],
                 boundaries: Tuple[int],
                 x_axis: np.ndarray,
                 y_axis: np.ndarray,
                 dth: float
                 ):

        #todo context manager for temporaries

        self.elarche(
            idx=x_axis,
            jdx=y_axis,
            dx=grid_spacings[0],
            dy=grid_spacings[1],
            dth=dth,
            bcx_kind=boundaries[0],
            bcy_kind=boundaries[1],
            vx_e=state["vx_e"],
            vy_e=state["vy_e"],
            vx=state["vx"],
            vy=state["vy"],
            lx=tmps["lx"],
            ly=tmps["ly"],
            i_dep=tmps["i_dep"],
            j_dep=tmps["j_dep"]
        )

        # Interpolate
        self.elarche.d_interpolate_lin_2d(
            psi=state["tracer_e"],
            psi_dep=state["tracer"],
            lx=tmps["lx"],
            ly=tmps["ly"],
            i_dep=tmps["i_dep"],
            j_dep=tmps["j_dep"],
            **self.elarche.symbol_mapping
        )

