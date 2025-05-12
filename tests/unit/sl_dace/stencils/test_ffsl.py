import pytest
from gt4py.cartesian.gtscript import stencil
import numpy as np

from tests.conftest import get_cpu_backends
from sl_dace.stencils.ffsl import (
    velocity_on_faces_x,
    split_cfl_x,
    fourth_order_facet_interpolation_x,
    integer_and_fractional_flux_sum
)


@pytest.mark.parametrize("backend", get_cpu_backends())
def test_velocity_on_faces(domain, inner_domain, inner_domain_origin, backend, dtypes):

    velocity_on_faces_x_stencil = stencil(
        name="velocity_on_faces",
        definition=velocity_on_faces_x,
        backend=backend,
        dtypes=dtypes,
        rebuild=True
    )

    # note : the first index (0) for plain cells is not relevant
    # implemented to align plain and half levels
    # n(half_levels) = n(plain_levels) + 1
    state = {
        "vx": np.ones(domain),
        "vy": np.ones(domain),
        "vxh": np.zeros(domain),
        "vyh": np.zeros(domain)
    }

    velocity_on_faces_x_stencil(
        **state,
        domain=inner_domain,
        origin=inner_domain_origin,
    )

    assert state["vxh"][1:-1,:, :].any()== 1

@pytest.mark.parametrize("backend", get_cpu_backends())
def test_split_cfl_x(domain, inner_domain, inner_domain_origin, backend, dtypes):

    split_cfl_x_stencil = stencil(
        name="split_cfl_x",
        definition=split_cfl_x,
        backend=backend,
        dtypes=dtypes,
        rebuild=True
    )

    state = {
        "vxh": 1.5 * np.ones(domain),
        "cxh_int": np.zeros(domain),
        "cxh_frac": np.zeros(domain)
    }

    split_cfl_x_stencil(
        **state,
        dx=1.0,
        dt=1.0,
        domain=inner_domain,
        origin=inner_domain_origin
    )

    assert state["cxh_frac"][1:-1, 1:-1, 1:-1].mean() == -0.5
    assert state["cxh_int"][1:-1, 1:-1, 1:-1].mean() == -1

@pytest.mark.parametrize("backend", get_cpu_backends())
def test_fourth_order_facet_interpolation(domain, inner_domain, inner_domain_origin, backend, dtypes):

    fofi_x_stencil = stencil(
        name="fofi_x",
        definition=fourth_order_facet_interpolation_x,
        backend=backend,
        dtypes=dtypes
    )

    state = {
        "psi": np.ones(domain),
        "psih": np.zeros(domain)
    }

    fofi_x_stencil(
        **state,
        domain=inner_domain,
        origin=inner_domain_origin
    )

    assert state["psih"][1:-1, 1:-1, 1:-1].mean() == 1.0

    fofi_y_stencil = stencil(
        name="fofi_x",
        definition=fourth_order_facet_interpolation_x,
        backend=backend,
        dtypes=dtypes
    )

    state = {
        "psi": np.ones(domain),
        "psih": np.zeros(domain)
    }

    assert state["psih"][1:-1, 1:-1, 1:-1].mean() == 1.0


@pytest.mark.parametrize("backend", get_cpu_backends())
def test_integer_and_fractional_flux_sum(domain, inner_domain, inner_domain_origin, backend, dtypes):

    iffs_stencil = stencil(
        name="iffs",
        definition=integer_and_fractional_flux_sum,
        backend=backend,
        dtypes=dtypes
    )

    state = {
        "fhx": np.zeros(domain),
        "fhx_frac": np.ones(domain),
        "fhx_int": np.ones(domain)
    }

    iffs_stencil(
        **state,
        domain=inner_domain,
        origin=inner_domain_origin
    )

    assert (
            state["fhx"][1:-1, 1:-1, 1:-1]
            - (
                    state["fhx_int"][1:-1, 1:-1, 1:-1]
                    + state["fhx_frac"][1:-1, 1:-1, 1:-1])
           ).mean() == 0.0

