from ifs_physics_common.framework.storage import managed_temporary_storage
from ifs_physics_common.framework.components import DiagnosticComponent
from ifs_physics_common.framework.grid import I, J, K, ComputationalGrid
from ifs_physics_common.framework.stencil import compile_stencil
from gt4py.cartesian.gtscript import stencil
from functools import partial, cached_property
import numpy as np
from config import Config
from sl_dace.stencils.ffsl import (
    velocity_on_faces_x,
    split_cfl_x,
    fourth_order_facet_interpolation_x,
    monotonic_limiter_x,
    fractional_flux_x,
    inner_density_update_x,
    fractional_and_integer_flux_sum
)
from sl_dace.ppm import ppm_coefficients_x
from sl_dace.interpolation.flux_integral import (
    integer_flux_integral_x,
    fractional_flux_integral_x
)


class FluxFormSemiLagX:

    def __init__(self, config: Config):
        self.config = config
        self.backend = self.config.gt4py_backend
        self.domain = self.config.domain

        nx, ny, nz =  self.config.domain

        # (optim) : contract stencil compilation
        self.velocity_on_faces_x = stencil(
            name="velocity_on_faces_x",
            definition=velocity_on_faces_x,
            backend=self.backend
        )
        self.split_cfl_x = stencil(
            name="split_cfl_x",
            definition=split_cfl_x,
            backend=self.backend
        )
        self.tracer_interpolation_x = stencil(
            name="tracer_interpolation_x",
            definition=fourth_order_facet_interpolation_x,
            backend=self.backend
        )
        self.monotonic_limiter_x = stencil(
            name="monotonic_limiter_x",
            definition=monotonic_limiter_x,
            backend=self.backend,
        )
        self.ppm_coefficients_x = stencil(
            name="ppm_coefficients_x",
            definition=ppm_coefficients_x,
            backend=self.backend
        )
        self.flux_sum = stencil(
            name="fractional_and_integer_flux_sum",
            defintion=fractional_and_integer_flux_sum,
            backend=self.backend
        )
        self.inner_density_update_x = stencil(
            name="inner_density_update_x",
            definition=inner_density_update_x,
            backend=self.backend
        )

        # interpolations numpy / dace
        # integer_flux_integral
        self.integer_flux_integral_x = partial(
            integer_flux_integral_x,
            nx=nx,
            ny=ny,
            nz=nz
        )
        # fractional flux integral
        self.fractional_flux_integral_x = partial(
            fractional_flux_integral_x,
            nx=nx,
            ny=ny,
            nz=nz
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
            "vxh": {"grid": (I - 1/2, J, K), "type": "float"},
            "cxh_int": {"grid": (I - 1/2, J, K), "type": "int"},
            "cxh_frac": {"grid": (I - 1/2, J, K), "type": "float"},
            "a0x": {"grid": (I, J, K), "type": "float"},
            "a1x": {"grid": (I, J, K), "type": "float"},
            "a2x": {"grid": (I, J, K), "type": "float"},
            "fhx_int": {"grid": (I - 1/2, J, K), "type": "float"},
            "fhx_frac": {"grid": (I - 1/2, J, K), "type": "float"},
            "fhx_x": {"grid": (I - 1/2, J, K), "type": "float"}
        }

    def __call__(self,
                 vx: np.ndarray,
                 vy: np.ndarray,
                 rho0: np.ndarray,
                 rho1: np.ndarray,
                 dt: float,
                 ds_yz: float,
                 dv: float
                 ):

        # todo : temporary managed storage
        with managed_temporary_storage() as (
            vxh,
            cxh_int,
            cxh_frac,
            rhohx,
            a0x,
            a1x,
            a2x
        ):

            # Cell faces remapping of v
            self.velocity_on_faces_x(
                vx=vx,
                vxh=vxh,
                domain=self.domain,
                origin=(0, 0, 0)
            )
            self.split_cfl_x(
                vxh=vxh,
                cxh_int=cxh_int,
                cxh_frac=cxh_frac,
                dx=dx,
                dt=dt,
                domain=self.domain,
                origin=(0, 0, 0)
            )
            self.tracer_interpolation_x(
                psihx=rhohx,
                psi=rho0,
                domain=self.domain,
                origin=(0, 0, 0)
            )
            self.monotonic_limiter_x(
                psi=rho0,
                psihx=rhohx,
                domain=self.domain,
                origin=(0, 0, 0)
            )
            self.ppm_coefficients_x(
                psi=rho0,
                psih=rhohx,
                a0=a0x,
                a1=a1x,
                a2=a2x,
                domain=self.domain,
                origin=(0, 0, 0)
            )

            # numpy / dace interpolations
            self.integer_flux_integral_x(
            fhx_int,
                chx_int=chx_int,
                fhx_int=rho,
                ds_sz=ds_zy,
                dv=dv,
                dt=dt,
            )
            self.fractional_flux_integral_x(
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
            self.flux_sum(
                fhx=fhx,
                fhx_int=fhx_int,
                fhx_frac=fhx_frac,
                domain=self.domain,
                origin=self.origin
            )
            self.inner_density_update_x(
                rho=rho0,
                rho_ix=rho1,
                fhx_x=fhx_x,
                ds_yz=ds_yz,
                dv=dv,
                dt=dt,
                domain=self.domain,
                origin=(0, 0, 0)
            )