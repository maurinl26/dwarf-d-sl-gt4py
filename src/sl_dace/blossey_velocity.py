import numpy as np
from gt4py.cartesian.gtscript import stencil

from sl_dace.utils.sdfg import build_sdfg
from sl_dace.stencils.blossey_velocity  import radius, theta, blossey_tracer
from fvms.initialization.velocity import _velocity_blossey


class Blossey:

    def __init__(self):

        self.d_blossey_stf = build_sdfg(blossey_stf)
        self.d_stream_function_xy = build_sdfg(stream_function_xy)
        self.d_blossey_tracer = build_sdfg(blossey_tracer)
        self.d_copy = build_sdfg(copy)

    def blossey_velocity(self,
                         stf: np.ndarray,
                         xcr: np.ndarray,
                         ycr: np.ndarray,
                         mt: float,
                         dx: float,
                         dy: float,
                         ):
        u, v = self.blossey_stf(stf, xcr, ycr, mt)
        self.stream_function_xy(u, v, stf, dx, dy)

    def init_blossey(self, xcr: np.ndarray, ycr: np.ndarray, mt: float, dt: float, dx: float, dy: float, nx: int, ny: int):
        # Vitesses
        self.blossey_velocity(vx, vy, xcr, ycr, mt, dx, dy)
        self.blossey_velocity(vx_p, vy_p, xcr, ycr, mt - dt, dx, dy)
        vx_e, vy_e = np.empty((nx, ny)), np.empty((nx, ny))

        return vx, vy, vx_p, vy_p, vx_e, vy_e

    def tracer_shape(self,
                       xcr: np.ndarray,
                       ycr: np.ndarray,
                       tracer: np.ndarray,
                       tracer_e: np.ndarray):
        # Tracer
        # todo: init tracer shape stencil
        self.blossey_tracer(tracer, xcr, ycr, 0)
        self.copy(tracer, tracer_e)



