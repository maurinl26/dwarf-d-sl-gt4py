import dace
import numpy as np
from sl_dace.utils.typingx import dtype_int, dtype_float
from sl_dace.utils.dims import I, J, K

# 1. remap velocity
def velocity_on_faces_x(
        vx: dtype_float[I, J, K],
        vxh: dtype_float[I, J, K + 1],
):
    """
    compute face velocity given mass point velocity
    with linear interpolation

    faces points are indexed from 0 to nx + 1
    mass points are indexed from 0 to nx

    mass point i is between face point i and face point i + 1

    :param vx: velocity on x (mass point)
    :param vy: velocity on y (mass point)
    :param vxh: velocity on x (face point)
    :param vyh: velocity on y (face point)
    """
    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        vxh[i, j, k] = 0.5 * (vx[i-1, j, k] + vx[i, j, k])

def velocity_on_faces_y(
        vy: dtype_float[I, J, K],
        vyh: dtype_float[I, J, K + 1],
):
    """
        compute face velocity given mass point velocity
        with linear interpolation

        faces points are indexed from 0 to nx + 1
        mass points are indexed from 0 to nx

        mass point i is between face point i and face point i + 1

        :param vy: velocity on y (mass point)
        :param vyh: velocity on y (face point)
        """

    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        vyh[i, j, k] = 0.5 * (vy[i, j-1, k] + vy[i, j-1, k])


# 1.5 convert velocity to integer and fractional cfl
def split_cfl_x(
        vxh: dtype_float[I, J, K + 1],
        cxh_int: dtype_int[I, J, K + 1],
        cxh_frac: dtype_float[I, J, K + 1],
        dx: dtype_float,
        dt: dtype_float
):
    """
    Compute the fractional and the integer parts of the cfl on the faces
    of the domain.

    :param vxh: velocity on faces
    :param cxh_int: integer part of the cfl
    :param cxh_frac: fractional part of the cfl
    :param dx: grid spacing on x
    :param dt: timestep
    """
    cxh = np.ndarray([I, J, K+1], dtype=dtype_float)

    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        cxh[i, j, k] =  - vxh[i, j, k] * dt / dx
        cxh_int[i, j, k] = floor(cxh[i, j, k])
        cxh_frac[i, j, k] = cxh[i, j, k] - floor(cxh[i, j, k])


def split_cfl_y(
        vyh: dtype_float[I, J, K + 1],
        cyh_int: dtype_int[I, J, K + 1],
        cyh_frac: dtype_float[I, J, K + 1],
        dy: dtype_float,
        dt: dtype_float
):
    """

    :param vyh:
    :param cyh_int:
    :param cyh_frac:
    :param dy:
    :param dt:
    """
    cyh = np.ndarray([I, J, K+1], dtype=dtype_float)

    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        cyh[i, j, k] =  - vyh[i, j, k] * dt / dy
        cyh_int[i, j, k] = floor(cyh[i, j, k])
        cyh_frac[i, j, k] = cyh[i, j, k] - floor(cyh[i, j, k])


# 2. Interpolate h (advected tracer)
def fourth_order_facet_interpolation_x(
        psihx: dtype_float[I, J, K + 1],
        psi: dtype_float[I, J, K]
):
    """
    4th order interpolation to compute facet values given mass point values

    :param psih: tracer field at facet points
    :param psi:  tracer field at mass points
    """
    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        psihx[i, j, k] = (psi[i-2, j, k] + 7 * psi[i-1, j, k] + 7 * psi[i, j, k] + psi[i+1, j, k]) / 12


def fourth_order_facet_interpolation_y(
        psihy: dtype_float[I, J, K + 1],
        psi: dtype_float[I, J, K]
):
        for i, j, k in dace.map[0:I, 0:J, 0:K]:
            psihy[i, j, k] = (psi[i, i-2, k] + 7 * psi[i, j-1, k] + 7 * psi[i, j, k] + psi[i, j+1, k]) / 12


# 4. PPM limiter
def monotonic_limiter_x(
        psihx: dtype_float[I, J, K + 1],
        psi: dtype_float[I, J, K]
):
    """
    Monotonic limiter for face values.
    Face value lies between neighbouring mass point values

    :param psih: tracer field at facet point
    :param psi: tracer field at mass point
    """
    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        psihx[i, j, k] = min(
            max(psi[i-1, j, k], psi[i, j, k]),
            max(psihx[i, j, k], min(psi[i-1, j, k], psi[i, j, k]))
        )


def monotonic_limiter_y(
    psihy: dtype_float[I, J, K+1],
    psi: dtype_float[I, J, K]
):
    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        psihy[i, j, k] = min(
            max(psi[i, j-1, k], psi[i, j, k]),
            max(psihy[i,j,k], min(psi[i, j-1, k], psi[i, j, k]))
        )


# note : integer part of the flux is computed in numpy
# sum of integer and fractional part also in numpy

# 4.5 sum of fractional and integer flux
def integer_and_fractional_flux_sum(
        fhx: dtype_float[I, J, K],
        fhx_int: dtype_float[I, J, K + 1],
        fhx_frac: dtype_float[I, J, K + 1]
):
    """
    Sum integer and fractional part of the flux.

    :param fhx: total flux
    :param fhx_int: integer part of the flux
    :param fhx_frac: fractional part of the flux
    """
    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        fhx[i, j, k] = fhx_frac[i, j, k] + fhx_int[i, j, k]


# 5. density update with updated flux
def inner_density_update_x(
        rho: dtype_float[I, J, K],
        rho_ix: dtype_float[I, J, K],
        fhx: dtype_float[I, J, K + 1],
        ds_yz: dtype_float,
        dv: dtype_float,
        dt: dtype_float,
):
    """
    Compute the inner density update on x axis

    :param rho: density to update
    :param fhx: flux at cell face
    :param ds_yz: section of cell face
    :param dv: volume of cell
    :param dt: time step
    """
    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        rho_ix[i, j, k] = (
                rho[i, j, k]
                - dt * (fhx[i, j, k] * ds_yz - fhx[i, j, k] * ds_yz) / dv
        )


def inner_density_update_y(
        rho: dtype_float[I, J, K],
        rho_iy: dtype_float[I, J, K],
        fhy: dtype_float[I, J, K + 1],
        ds_xz: dtype_float,
        dv: dtype_float,
        dt: dtype_float
):
    """
    Compute the inner density update on y axis

    :param rho: density
    :param fhy: flux on face on y
    :param ds_xz: cell's normal section to y axis
    :param dv: volume of cell
    :param dt: time step
    """
    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        rho_iy[i, j, k] = (
                rho[i, j, k]
                - dt * (fhy[i, j, k] * ds_xz - fhy[i, j, k] * ds_xz) / dv
        )


# todo : implement SWIFT outer steps for splitting
def swift_outer_density_update(
        rho_ay: dtype_float[I, J, K],
        sigma_x
):
    ...


