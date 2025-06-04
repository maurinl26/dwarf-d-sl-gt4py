import logging
from functools import cached_property

import numpy as np
from config import Config
from sl_dace.interpolation.interpolation_2d import (
    interpolate_lin_2d,
)
from sl_dace.elarche import Elarche
from sl_dace.utils.dims import I, J, K
import dace
from contextlib import contextmanager
from ifs_physics_common.framework.grid import ComputationalGrid, I, J, K

from typing import Tuple, Dict

logging.getLogger(__name__)

# note : similar to ifs_physics_common
@contextmanager
def temporaries(
        temporary_fields: dict,
        grids: ComputationalGrid
):

    yield {
        name: np.zeros(grids[field_desc["grid"]].shape, dtype=field_desc["dtype"])
        for name, field_desc in temporary_fields.items()
    }

class SmiLagXY:

    def __init__(self, computational_grid: Dict, nitmp: int = 4):

        self.computational_grid = computational_grid
        grid_shape = computational_grid.grids[(I, J, K)].shape
        self.elarche = Elarche(grid_shape, nitmp)

    @cached_property
    def _input_properties(self):
        return {
            "vx": {"grid": (I, J, K), "dtype": float},
            "vy": {"grid": (I, J, K), "dtype": float},
            "tracer": {"grid": (I, J, K), "dtype": float}
        }

    @cached_property
    def _inout_properties(self):
        return {
            "vx_e": {"grid": (I, J, K), "dtype": float},
            "vy_e": {"grid": (I, J, K), "dtype": float},
            "tracer_e": {"grid": (I, J, K), "dtype": float}
        }

    @cached_property
    def _temporaries(self):
        return {
            "lx": {"grid": (I, J, K), "dtype": float},
            "ly": {"grid": (I, J, K), "dtype": float},
            "i_dep": {"grid": (I, J, K), "dtype": int},
            "j_dep": {"grid": (I, J, K), "dtype": int},
            "vx_tmp": {"grid": (I, J, K), "dtype": float},
            "vy_tmp": {"grid": (I, J, K), "dtype": float}
        }

    def __call__(self,
                 state: Dict[str, np.ndarray],
                 grid_spacings: Tuple[float],
                 boundaries: Tuple[int],
                 x_axis: np.ndarray,
                 y_axis: np.ndarray,
                 dth: float
                 ):

        with temporaries(self._temporaries, self.computational_grid.grids) as tempo:

            self.elarche(
                idx=x_axis,
                jdx=y_axis,
                dx=grid_spacings[0],
                dy=grid_spacings[1],
                dth=dth,
                vx_e=state["vx_e"],
                vy_e=state["vy_e"],
                vx=state["vx"],
                vy=state["vy"],
                lx=tempo["lx"],
                ly=tempo["ly"],
                i_dep=tempo["i_dep"],
                j_dep=tempo["j_dep"],
                vx_tmp=tempo["vx_tmp"],
                vy_tmp=tempo["vy_tmp"]
            )

            # Interpolate
            self.elarche.d_interpolate_lin_2d(
            psi=state["tracer_e"],
            psi_dep=state["tracer"],
            lx=tempo["lx"],
            ly=tempo["ly"],
            i_dep=tempo["i_dep"],
            j_dep=tempo["j_dep"],
            **self.elarche.symbol_mapping
        )

