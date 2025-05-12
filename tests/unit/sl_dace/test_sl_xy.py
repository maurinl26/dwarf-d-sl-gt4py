import pytest
from sl_dace.sl_xy import SmiLagXY
import numpy as np
from ifs_physics_common.framework.grid import I, J, K


def test_sl_xy(computational_grid, config):

    grid_shape = computational_grid.grids[(I, J, K)].shape

    grid_spacings = (config.dx, config.dy, config.dz)

    x_axis = np.arange(0, nx)
    y_axis = np.arange(0, ny)
    dth = 0.5


    smi_lag_xy = SmiLagXY(
        computational_grid=computational_grid
    )

    state = {
        "vx": np.ones(grid_shape),
        "vy": np.ones(grid_shape),
        "tracer": np.ones(grid_shape)
    }

    smi_lag_xy(
        state=state,
        grid_spacings=grid_spacings,
        dth=dth,
        x_axis=x_axis,
        y_axis=y_axis
    )

    assert False
