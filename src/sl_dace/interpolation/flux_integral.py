import numpy as np
from itertools import product
from typing import Tuple
from sl_dace.utils.typingx import dtype_float, dtype_int
from sl_dace.utils.dims import I, J, K

@dace.program
def integer_flux_integral_x(
    fhx_int: dtype_float[I + 1, J, K],
    chx_int: dtype_float[I + 1, J, K],
    rho: dtype_float[I, J, K],
    ds_yz: dtype_float,
    dv: dtype_float,
    dt: dtype_float,
):
    """
    Computes the integral of the flux on the path given by c
    the integer cfl on x.

    The path is constituted by upstream cells upto chx_int.

    :param fhx_int: flux at face
    :param chx_int: integer part of the cfl
    :param rho: density to advect
    :param ds_yz: section of the cell
    :param dv:  volume of the cell
    :param dt:  timestep
    :param nx:  number of cells on x
    :param ny:  number of cells on y
    :param nz:  number of cells on z
    """


    for i,j,k in dace.map[0:I, 0:J, 0:K]:

            if chx_int[i, j, k] >= 0:
                fhx_int[i, j, k] = np.sum(
                    rho[i - chx_int[i, j, k]: i, j, k] * dv
                ) / (ds_yz * dt)

            else:
                fhx_int[i, j, k] = np.sum(
                    rho[i + 1: i - chx_int[i, j, k], j, k] * dv
                ) / (ds_yz * dt)

@dace.program
def fractional_flux_integral_x(
    a0: dtype_float[I, J, K],
    a1: dtype_float[I, J, K],
    a2: dtype_float[I, J, K],
    chx_int: dtype_float[I + 1, J, K],
    chx_frac: dtype_float[I + 1, J, K],
    fhx_frac: dtype_float[I + 1, J, K],
    ds_yz: dtype_float,
    dv: dtype_float,
    dt: dtype_float,
):
    """
    Computes the integral of the flux on the path given by c
    the integer cfl on x.

    The path is constituted by upstream cells upto chx_int.

    :param fhx_int: flux at face
    :param chx_int: integer part of the cfl
    :param rho: density to advect
    :param ds_yz: section of the cell
    :param dv:  volume of the cell
    :param dt:  timestep
    :param nx:  number of cells on x
    :param ny:  number of cells on y
    :param nz:  number of cells on z
    """

    for i,j,k in dace.map[0:I, 0:J, 0:K]:

            if chx_int[i, j, k] >= 0:
                ppm_coeffs = (
                    a0[i - chx_int[i, j, k] - 1,j,k],
                    a1[i - chx_int[i, j, k] - 1, j, k],
                    a2[i - chx_int[i, j, k] - 1, j, k]
                )

                fhx_frac[i, j, k] = (
                    ppm_integral(ppm_coeffs, 1.0 - chx_frac[i, j, k], 1.0) * dv
                    / (ds_yz * dt)
                )

            else:
                ppm_coeffs = (
                    a0[i - chx_int[i,j,k] - 1, j, k],
                    a1[i - chx_int[i,j,k] - 1, j, k],
                    a2[i - chx_int[i,j,k] - 1, j, k]
                )

                fhx_frac[i, j, k] = (
                    ppm_integral(ppm_coeffs, 0, abs(chx_frac[i, j, k]), ) * dv
                ) / (ds_yz * dt)

@dace.program
def ppm_integral(
        a: Tuple[dtype_float],
        inf_bound: dtype_float,
        sup_bound: dtype_float
):
    """
    Integral of the parabolic spline of ppm
    given a lower and an upper bound.

    :param a: tuple of parabolic coefficients
    :param x:
    :param inf_bound:
    :param sup_bound:
    :return:
    """
    prim = lambda ppm, x: (
        (ppm[2] * x ** 3 / 3)
        + (ppm[1] * x ** 2 / 2)
        + (ppm[0] * x)
    )

    return (
        prim(a, sup_bound)
        - prim(a, inf_bound)
    )



