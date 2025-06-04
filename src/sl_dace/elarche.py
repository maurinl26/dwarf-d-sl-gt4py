from typing import Tuple

import dace

from sl_dace.interpolation.interpolation_2d import interpolate_lin_2d
from sl_dace.stencils.copy import copy
from sl_dace.stencils.dep_search_1d import dep_search_1d
from sl_dace.utils.dims import I, J, K
from sl_dace.utils.sdfg import build_sdfg
from sl_dace.utils.typingx import dtype_float, dtype_int


class Elarche:

    def __init__(self, grid_shape: Tuple[int], halo: int = 0, nitmp: int = 4, backend: str = "dace:cpu"):
        self.grid_shape = grid_shape
        self.symbol_mapping = {
                "I": grid_shape[0],
                "J": grid_shape[1],
                "K": grid_shape[2],
                "H": halo
            }
        self.nitmp = nitmp
        self.d_copy = build_sdfg(copy, device="cpu", mode="aot")
        self.d_interpolate_lin_2d = build_sdfg(interpolate_lin_2d, device="cpu", mode="aot")
        self.d_dep_search_1d = build_sdfg(dep_search_1d, device="cpu", mode="aot")


    @dace.method
    def __call__(self,
                 dx: dtype_float,
                 dy: dtype_float,
                 dth: dtype_float,
                 idx: dtype_int[I, J, K],
                 jdx: dtype_int[I, J, K],
                 vx: dtype_float[I, J, K],
                 vy: dtype_float[I, J, K],
                 vx_tmp: dtype_float[I, J, K],
                 vy_tmp: dtype_float[I, J, K],
                 vx_e: dtype_float[I, J, K],
                 vy_e: dtype_float[I, J, K],
                 # Outputs
                 lx: dtype_float[I, J, K],
                 ly: dtype_float[I, J, K],
                 i_dep: dtype_int[I, J, K],
                 j_dep: dtype_int[I, J, K],
        ):

        # Temporaries
        self.d_copy(
            vx=vx,
            vy=vy,
            vx_tmp=vx_tmp,
            vy_tmp=vy_tmp,
            **self.symbol_mapping
        )

        # Array declaration
        for l in range(self.nitmp):
            self.d_dep_search_1d(
                vx_e=vx_e,
                vx_tmp=vx_tmp,
                i_a=idx,
                i_d=i_dep,
                lx=lx,
                dx=dx,
                dth=dth,
                **self.symbol_mapping
            )
            self.d_dep_search_1d(
                vx_e=vy_e,
                vx_tmp=vy_tmp,
                i_a=jdx,
                i_d=jdx,
                lx=ly,
                dx=dy,
                dth=dth,
                **self.symbol_mapping
            )

            # todo: add LipschitzDiag
            ####### Interpolation for fields ########
            self.d_interpolate_lin_2d(
                psi_dep=vx,
                psi=vx_tmp,
                lx=lx,
                ly=ly,
                i_dep=i_dep,
                j_dep=j_dep,
                **self.symbol_mapping
            )

            self.d_interpolate_lin_2d(
                psi_dep=vy,
                psi=vy_tmp,
                lx=lx,
                ly=ly,
                i_dep=i_dep,
                j_dep=j_dep,
                **self.symbol_mapping
            )

        return lx, ly, i_dep, j_dep

    def interpolate_tracer_field(
            self,
            tracer_dep: dtype_float[I, J, K],
            tracer: dtype_float[I, J, K],
            lx: dtype_float[I, J, K],
            ly: dtype_float[I, J, K],
            i_dep: dtype_int[I, J, K],
            j_dep: dtype_int[I, J, K]
    ):

        self.d_interpolate_lin_2d(
            psi_dep=tracer_dep,
            psi=tracer,
            lx=lx,
            ly=ly,
            i_dep=i_dep,
            j_dep=j_dep,
            **self.symbol_mapping
        )


