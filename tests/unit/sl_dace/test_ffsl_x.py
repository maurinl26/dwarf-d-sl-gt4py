import pytest
import numpy as np
from sl_dace.ffsl_x import FluxFormSemiLagX

def test_ffsl_x(
    config,
    domain,
    computational_grid,
    gt4py_config
):

    state =  {
        "vx": np.ones(domain),
        "vy": np.ones(domain),
        "rho0": np.ones(domain),
        "rho1": np.zeros(domain)
    }

    geometry = {
        "ds_yz": 1.0,
        "dv": 1.0,
        "dx": 1.0,
        "dt": 1.0
    }

    FluxFormSemiLagX(
        config=config,
        gt4py_config=gt4py_config,
        computational_grid=computational_grid
    )(
        **state,
        **geometry
    )

    assert state["rho1"].mean() == 1.0

