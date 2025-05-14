import numpy as np
from itertools import product
from typing import Tuple

def integer_flux_integral_x(
    fhx_int: np.ndarray,
    chx_int: np.ndarray,
    rho: np.ndarray,
    ds_yz: float,
    dv: float,
    dt: float,
    nx: np.ndarray,
    ny: np.ndarray,
    nz: np.ndarray
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
    I = np.arange(0, nx)
    J = np.arange(0, ny)
    K = np.arange(0, nz)

    for i,j in product(I, J):
        for k in K:

            if chx_int[i, j, k] >= 0:
                fhx_int[i, j, k] = np.sum(
                    rho[i - chx_int[i, j, k]: i, j, k] * dv
                ) / (ds_yz * dt)

            else:
                fhx_int[i, j, k] = np.sum(
                    rho[i + 1: i - chx_int[i, j, k], j, k] * dv
                ) / (ds_yz * dt)


def fractional_flux_integral_x(
    a0: np.ndarray,
    a1: np.ndarray,
    a2: np.ndarray,
    chx_int: np.ndarray,
    chx_frac: np.ndarray,
    fhx_frac: np.ndarray,
    ds_yz: float,
    dv: float,
    dt: float,
    nx: np.ndarray,
    ny: np.ndarray,
    nz: np.ndarray
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
    I = np.arange(0, nx)
    J = np.arange(0, ny)
    K = np.arange(0, nz)

    for i,j in product(I, J):
        for k in K:

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


def ppm_integral(
        a: Tuple[float],
        inf_bound: float,
        sup_bound: float):
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



