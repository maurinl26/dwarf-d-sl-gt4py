import logging
from functools import cached_property, partial
from itertools import repeat

import dace
import numpy as np
from gt4py.cartesian.gtscript import stencil
from ifs_physics_common.framework.config import GT4PyConfig
from ifs_physics_common.framework.storage import managed_temporary_storage

from config import Config
from sl_dace.interpolation.flux_integral import (fractional_flux_integral_x,
                                                 integer_flux_integral_x)
from sl_dace.stencils.ffsl import (fourth_order_facet_interpolation_x,
                                   inner_density_update_x,
                                   integer_and_fractional_flux_sum,
                                   monotonic_limiter_x, split_cfl_x,
                                   velocity_on_faces_x)
from sl_dace.stencils.ppm import ppm_coefficients_x

from sl_dace.utils.typingx import dtype_int, dtype_float
from sl_dace.utils.dims import I, J, K
from sl_dace.utils.sdfg import build_sdfg

class FluxFormSemiLagX:

    def __init__(self, config: Config, gt4py_config: GT4PyConfig):
        self.config = config
        self.gt4py_config = gt4py_config
        self.computational_grid = self.config.computational_grid

        self.domain = self.config.domain
        self.inner_domain = (
            self.domain[0] - 1,
            self.domain[1] - 1,
            self.domain[2]
        )

        logging.warning(f"Inner domain shape {self.inner_domain}")

        nx, ny, nz =  self.config.domain

        # (optim) : contract stencil compilation
        self.d_velocity_on_faces_x = build_sdfg(velocity_on_faces_x)
        self.d_split_cfl_x = build_sdfg(split_cfl_x)
        self.d_tracer_interpolation_x = build_sdfg(fourth_order_facet_interpolation_x)
        self.d_monotonic_limiter_x = build_sdfg(monotonic_limiter_x)
        self.d_ppm_coefficients_x = build_sdfg(ppm_coefficients_x)
        self.d_flux_sum = build_sdfg(integer_and_fractional_flux_sum)
        self.d_inner_density_update_x = build_sdfg(inner_density_update_x)
        self.d_integer_flux_integral_x = build_sdfg(
            integer_flux_integral_x,
        )
        # fractional flux integral
        self.fractional_flux_integral_x = build_sdfg(
            fractional_flux_integral_x,
        )


    # todo: import dims
    @cached_property
    def _input_properties(self):
        return {
            "rho0": {"grid": (I, J, K), "type": "float"},
        }

    @cached_property
    def _output_properties(self):
        return {
            "rho1": {"grid": (I, J, K), "type": "float"}
        }

    @cached_property
    def _temporaries(self):
        return {
            "vxh": {"grid": (I - 1, J, K), "type": "float"},
            "cxh_int": {"grid": (I - 1, J, K), "type": "int"},
            "cxh_frac": {"grid": (I - 1, J, K), "type": "float"},
            "a0x": {"grid": (I, J, K), "type": "float"},
            "a1x": {"grid": (I, J, K), "type": "float"},
            "a2x": {"grid": (I, J, K), "type": "float"},
            "fhx_int": {"grid": (I - 1, J, K), "type": "float"},
            "fhx_frac": {"grid": (I - 1, J, K), "type": "float"},
            "fhx_x": {"grid": (I - 1, J, K), "type": "float"}
        }

    # todo: shift call to dace
    # @dace.method
    def __call__(self,
                 vx: np.ndarray,
                 vy: np.ndarray,
                 rho0: np.ndarray,
                 rho1: np.ndarray,
                 dt: float,
                 ds_yz: float,
                 dv: float,
                 dx: float
                 ):

        with managed_temporary_storage(
            self.computational_grid,
            ((I, J, K), "int"),
            *repeat(((I, J, K), "float"), 6),
            *repeat(((I, J, K), "float"), 3),
            gt4py_config=self.gt4py_config
        ) as (
            chx_int,
            vxh,
            chx_frac,
            rho_hx,
            fhx,
            fhx_int,
            fhx_frac,
            a0x,
            a1x,
            a2x,
        ):

            logging.warning(f"vx, shape {vx.shape}")
            logging.warning(f"vxh, shape {vxh.shape}")

            # Cell faces remapping of v
            self.d_velocity_on_faces_x(
                vx=vx[:-1, :,:],
                vxh=vxh[:-1, :,:],
                domain=self.domain,
                origin=(1, 0, 0)
            )

            # todo: boundary conditions


            self.d_split_cfl_x(
                vxh=vxh,
                cxh_int=chx_int,
                cxh_frac=chx_frac,
                dx=dx,
                dt=dt,
                domain=self.inner_domain,
                origin=(1, 0, 0)
            )

            interpolation_domain = (
                self.inner_domain[0] - 2,
                self.inner_domain[1],
                self.inner_domain[2]
            )

            self.d_tracer_interpolation_x(
                psihx=rho_hx,
                psi=rho0,
                domain=interpolation_domain,
                origin=(2, 0, 0)
            )
            self.d_monotonic_limiter_x(
                psi=rho0,
                psihx=rho_hx,
                domain=self.inner_domain,
                origin=(1, 0, 0)
            )


            self.d_ppm_coefficients_x(
                psi=rho0,
                psih=rho_hx,
                a0=a0x,
                a1=a1x,
                a2=a2x,
                domain=interpolation_domain,
                origin=(2, 0, 0)
            )

            self.d_integer_flux_integral_x(
                fhx_int=fhx_int,
                chx_int=chx_int,
                rho=rho0,
                ds_yz=ds_yz,
                dv=dv,
                dt=dt,
            )
            self.d_fractional_flux_integral_x(
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
            
            #sum and density update
            self.d_flux_sum(
                fhx=fhx,
                fhx_int=fhx_int,
                fhx_frac=fhx_frac,
                domain=self.inner_domain,
                origin=(1, 0, 0)
            )
            self.d_inner_density_update_x(
                rho=rho0,
                rho_ix=rho1,
                fhx=fhx,
                ds_yz=ds_yz,
                dv=dv,
                dt=dt,
                domain=self.inner_domain,
                origin=(1, 0, 0)
            )

