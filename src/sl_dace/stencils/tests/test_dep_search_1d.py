from gt4py.cartesian.gtscript import stencil
import numpy as np
import pytest
import logging

from sl_dace.stencils.dep_search_1d import dep_search_1d
from typing import Union

@pytest.fixture
def dtypes(precision: Union["single", "double"] = 'single'):
    if precision == "single":
        return {
            float: np.float32,
            int: np.int32,
        }
    elif precision == "double":
        return {
            float: np.float64,
            int: np.int64,
        }
    else:
        raise KeyError("precision not in 'single' or 'double'")


@pytest.mark.parametrize("backend", ["numpy", "gt:cpu_ifirst", "dace:cpu"] )
def test_dep_search_1d(backend: str, dtypes: dict):

    logging.info(f"Backend : {backend}")
    for key, type in dtypes.items():
        logging.info(f"{key}, Type : {type}")
    stencil_dep_search_1d = stencil(backend=backend,
                                    definition=dep_search_1d,
                                    name="dep_search_1d",
                                    dtypes=dtypes,
                                    build_info={},
                                    externals={},
                                    rebuild=True,
                                    )

    vx_e = np.ones((50, 50, 10), dtype=dtypes[float])
    vx_tmp = np.ones((50, 50, 10), dtype=dtypes[float])
    i_a = np.ones((50, 50, 10), dtype=dtypes[int])
    i_d = np.zeros((50, 50, 10), dtype=dtypes[int])
    lx = np.zeros((50, 50, 10), dtype=dtypes[float])

    dx = dtypes[float](0.01)
    dth = dtypes[float](0.5)

    stencil_dep_search_1d(
        vx_e=vx_e,
        vx_tmp=vx_tmp,
        i_a=i_a,
        i_d=i_d,
        lx=lx,
        dx=dx,
        dth=dth,
    )

# Test dace functionalities
@pytest.mark.parametrize('backend', ['dace:cpu'])
def test_orchestrate(dtypes: dict, backend: str):
    

    stencil_dep_search_1d = stencil(backend=backend,
                                    definition=dep_search_1d,
                                    name="dep_search_1d",
                                    dtypes=dtypes,
                                    build_info={},
                                    externals={},
                                    rebuild=True,
                                    )

    print(stencil_dep_search_1d._sdfg)

