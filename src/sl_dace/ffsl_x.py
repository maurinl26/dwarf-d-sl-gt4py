import logging

import dace
import numpy as np
from gt4py.cartesian.gtscript import stencil
from ifs_physics_common.framework.config import GT4PyConfig
from ifs_physics_common.framework.storage import managed_temporary_storage

from config import Config
from sl_dace.interpolation.flux_integral import (
    fractional_flux_integral_x,
    integer_flux_integral_x
)
from sl_dace.stencils.ffsl import (
        fourth_order_facet_interpolation_x,
        inner_density_update_x,
        integer_and_fractional_flux_sum,
        monotonic_limiter_x,
        split_cfl_x,
        velocity_on_faces_x
)
from sl_dace.stencils.ppm import ppm_coefficients_x

from sl_dace.utils.typingx import dtype_int, dtype_float
from sl_dace.utils.dims import I, J, K
from sl_dace.utils.sdfg import build_sdfg


@dace.program
def flux_from_semi_lag_x(self,
             vx: np.ndarray,
             vy: np.ndarray,
             rho0: np.ndarray,
             rho1: np.ndarray,
             dt: float,
             ds_yz: float,
             dv: float,
             dx: float
             ):
    chx_int = np.ndarray([I + 1, J, K], dtype=dtype_float)
    vxh = np.ndarray([I + 1, J, K], dtype=dtype_float)
    chx_frac = np.ndarray([I + 1, J, K], dtype=dtype_float)
    rho_hx = np.ndarray([I + 1, J, K], dtype=dtype_float)
    fhx = np.ndarray([I + 1, J, K], dtype=dtype_float)
    fhx_int = np.ndarray([I + 1, J, K], dtype=dtype_float)
    fhx_frac = np.ndarray([I + 1, J, K], dtype=dtype_float)
    a0x = np.ndarray([I, J, K], dtype=dtype_float)
    a1x = np.ndarray([I, J, K], dtype=dtype_float)
    a2x = np.ndarray([I, J, K], dtype=dtype_float)

    logging.warning(f"vx, shape {vx.shape}")
    logging.warning(f"vxh, shape {vxh.shape}")

    velocity_on_faces_x(
        vx=vx[:-1, :, :],
        vxh=vxh[:-1, :, :],
    )

    # todo: boundary conditions
    split_cfl_x(
        vxh=vxh,
        cxh_int=chx_int,
        cxh_frac=chx_frac,
        dx=dx,
        dt=dt,
    )

    interpolation_domain = (
        self.inner_domain[0] - 2,
        self.inner_domain[1],
        self.inner_domain[2]
    )

    tracer_interpolation_x(
        psihx=rho_hx,
        psi=rho0,
    )

    monotonic_limiter_x(
        psi=rho0,
        psihx=rho_hx,
    )

    ppm_coefficients_x(
        psi=rho0,
        psih=rho_hx,
        a0=a0x,
        a1=a1x,
        a2=a2x,
    )

    integer_flux_integral_x(
        fhx_int=fhx_int,
        chx_int=chx_int,
        rho=rho0,
        ds_yz=ds_yz,
        dv=dv,
        dt=dt,
    )

    fractional_flux_integral_x(
        a0=a0x,
        a1=a1x,
        a2=a2x,
        chx_int=chx_int,
        chx_frac=chx_frac,
        fhx_frac=fhx_frac,
        ds_yz=ds_yz,
        dv=dv,
        dt=dt
    )

    # sum and density update
    flux_sum(
        fhx=fhx,
        fhx_int=fhx_int,
        fhx_frac=fhx_frac,
    )

    inner_density_update_x(
        rho=rho0,
        rho_ix=rho1,
        fhx=fhx,
        ds_yz=ds_yz,
        dv=dv,
        dt=dt,
    )

