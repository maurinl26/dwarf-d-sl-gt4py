from gt4py.cartesian.gtscript import stencil
import numpy as np
import pytest
import logging

from sl_dace.stencils.copy import copy, backup
from typing import Tuple


@pytest.mark.parametrize("backend", ["dace:cpu"])
def test_copy(backend: str, dtypes: dict, domain_with_halo: Tuple[int]):

    logging.info(f"Backend : {backend}")
    stencil_copy = stencil(backend=backend,
                                    definition=copy,
                                    name="copy",
                                    dtypes=dtypes,
                                    build_info={},
                                    externals={},
                                    rebuild=True,
                                    )

    vy_tmp = np.zeros((50, 50, 10), dtype=dtypes[float])
    vx_tmp = np.zeros((50, 50, 10), dtype=dtypes[float])

    vy = np.ones((50, 50, 10), dtype=dtypes[float])
    vx = np.ones((50, 50, 10), dtype=dtypes[float])

    stencil_copy(
        vx_tmp=vx_tmp,
        vy_tmp=vy_tmp,
        vx=vx,
        vy=vy,
        domain=domain_with_halo,
        origin=(0, 0, 0)
    )

    assert vx_tmp.mean() == 1
    assert vy_tmp.mean() == 1

@pytest.mark.parametrize("backend", ["dace:cpu"])
def test_backup(backend: str, dtypes: dict, domain: Tuple[int]):
    logging.info(f"Backend : {backend}")
    stencil_backup = stencil(backend=backend,
                           definition=backup,
                           name="backup",
                           dtypes=dtypes,
                           build_info={},
                           externals={},
                           rebuild=True,
                           )

    tracer = np.zeros((50, 50, 10), dtype=dtypes[float])
    tracer_e = np.ones((50, 50, 10), dtype=dtypes[float])

    stencil_backup(
        tracer=tracer,
        tracer_e=tracer_e,
        domain=domain,
        origin=(0, 0, 0)
    )

    assert tracer.mean() == 1


