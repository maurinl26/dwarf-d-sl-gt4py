from collections.abc import dict_values
from config import Config
from sl_dace.stencils.ffsl import (
    velocity_on_faces_y,
    split_cfl_y,
    fourth_order_facet_interpolation_y,
    monotonic_limiter_y,
    fractional_flux_y,
    inner_density_update_y,
    integer_and_fractional_flux_sum
)
from sl_dace.ppm import ppm_coefficients_y
from sl_dace.interpolation.flux_integral import (
    integer_flux_integral_y,
    fractional_flux_integral_y
)

from sl_dace.utils.typingx import dace_float
from sl_dace.utils.dims import I, J, K
import numpy as np
import logging

def flux_from_semi_lag_y(
        vy: dace_float[I, J, K],
        rho0: dace_float[I, J, K],
        rho1: dace_float[I, J, K],
        dt: dace_float,
        ds_xz: dace_float,
        dv: dace_float,
        dx: dace_float
):
    chy_int = np.ndarray([I + 1, J, K], dtype=dace_float)
    vyh = np.ndarray([I + 1, J, K], dtype=dace_float)
    chy_frac = np.ndarray([I + 1, J, K], dtype=dace_float)
    rho_hy = np.ndarray([I + 1, J, K], dtype=dace_float)
    fhy = np.ndarray([I + 1, J, K], dtype=dace_float)
    fhy_int = np.ndarray([I + 1, J, K], dtype=dace_float)
    fhy_frac = np.ndarray([I + 1, J, K], dtype=dace_float)
    a0y = np.ndarray([I, J, K], dtype=dace_float)
    a1y = np.ndarray([I, J, K], dtype=dace_float)
    a2y = np.ndarray([I, J, K], dtype=dace_float)

    logging.warning(f"vx, shape {vy.shape}")
    logging.warning(f"vxh, shape {vyh.shape}")

    velocity_on_faces_y(
        vx=vy[:-1, :, :],
        vxh=vyh[:-1, :, :],
    )

    # todo: boundary conditions
    split_cfl_y(
        vxh=vyh,
        cxh_int=chy_int,
        cxh_frac=chy_frac,
        dx=dx,
        dt=dt,
    )

    fourth_order_facet_interpolation_y(
        psihx=rho_hy,
        psi=rho0,
    )

    monotonic_limiter_y(
        psi=rho0,
        psihx=rho_hy,
    )

    ppm_coefficients_y(
        psi=rho0,
        psih=rho_hy,
        a0=a0y,
        a1=a1y,
        a2=a2y,
    )

    integer_flux_integral_y(
        fhx_int=fhy_int,
        chx_int=chy_int,
        rho=rho0,
        ds_yz=ds_xz,
        dv=dv,
        dt=dt,
    )

    fractional_flux_integral_y(
        a0=a0y,
        a1=a1y,
        a2=a2y,
        chx_int=chy_int,
        chx_frac=chy_frac,
        fhx_frac=fhy_frac,
        ds_yz=ds_xz,
        dv=dv,
        dt=dt
    )

    # sum and density update
    integer_and_fractional_flux_sum(
        fhx=fhy,
        fhx_int=fhy_int,
        fhx_frac=fhy_frac,
    )

    inner_density_update_y(
        rho=rho0,
        rho_ix=rho1,
        fhx=fhy,
        ds_xz=ds_xz,
        dv=dv,
        dt=dt,
    )



