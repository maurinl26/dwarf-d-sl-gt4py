from typing import Tuple
import numpy as np
import logging

from config import Config
from sl_python.interpolation import interpolate_lin_2d

logging.getLogger(__name__)


def dep_search_1d(
    I: np.ndarray, vx_e: np.ndarray, vx_tmp: np.ndarray, dx: np.ndarray, dth: float
) -> Tuple[np.ndarray]:
    """Compute departure point coordinate (1d)

    Args:
        I (np.ndarray): _description_
        vx_e (np.ndarray): velocity at arrival point (t + dt)
        vx_tmp (np.ndarray): estimate of velocity at departure point
        dx (np.ndarray): grid spacing
        dth (float): half model time step

    Returns:
        Tuple[np.ndarray]:
            i_d: indice of departure point on grid
            lx: adimensionned spacing of departure point from ref grid point
    """

    # Deplacement
    trajx = -dth * (vx_e + vx_tmp) / dx  # dth = dt / 2
    i_d = (I + np.floor(trajx)).astype(int)
    lx = trajx - np.floor(trajx)

    return lx, i_d


# Pointwise
def dep_search_z(
    i_d: int,
    j_d: int,
    k: int,
    vz_e: np.ndarray,
    vz_tmp: np.ndarray,
    dz: np.ndarray,
    dth: float,
) -> Tuple[int, float]:
    """Compute trajectory on an irregular grid.

    Args:
        i_d (int): point de départ i
        j_d (int): point de départ j
        k (int): point de départ k
        vz_e (np.ndarray): velocity at arrival point (t + dt)
        vz_tmp (np.ndarray): estimate of velocity at departure point
        dz (np.ndarray): grid spacing
        dth (float): half model time step

    Returns:
        Tuple[np.ndarray]:
            k: indice of departure point on grid
            lx: spacing of departure point to grid point k
    """

    trajz = -dth * (vz_e + vz_tmp)

    if trajz < 0:
        delta_z = 0
        k_search = k
        while delta_z < trajz:
            k_search -= 1
            delta_z -= dz[i_d, j_d, k_search]

        lambda_z = delta_z - trajz

    else:
        k_search = k
        delta_z = dz[k_search]
        while delta_z > trajz:
            k_search += 1
            delta_z += dz[i_d, j_d, k_search]

        lambda_z = trajz - delta_z

    return k_search, lambda_z


# ELARCHE
def lagrangian_search(
    config: Config,
    vx: np.ndarray,
    vz: np.ndarray,
    vx_e: np.ndarray,
    vz_e: np.ndarray,
    interpolation_function: callable,
    nitmp: int = 4,
) -> Tuple[np.ndarray]:
    """Research departure point for a given grid and velocity field.
    Terminates on nsiter iterations.

    Args:
        x (np.ndarray): grid of arrival points
        v (np.ndarray): velocity fields
        nsiter (int, optional): number of iterations. Defaults to 10.

    Returns:
        np.ndarray: departure point
    """

    vx_tmp = vx.copy()
    vz_tmp = vz.copy()

    # Array declaration
    for l in range(nitmp):
        lx, i_d = dep_search_1d(config.I, vx_e, vx_tmp, config.dx, config.dth)
        lambda_z, k_d = dep_search_z(config.K, vz_e, vz_tmp, config.dz, config.dth)

        ####### Interpolation for fields ########
        vx_tmp = interpolation_function(
            vx, lx, lambda_z, i_d, k_d, config.bcx_kind, config.bcy_kind, config.nx, config.ny
        )

        vz_tmp = interpolation_function(
            vz_tmp, lx, lambda_z, i_d, k_d, config.bcx_kind, config.bcy_kind, config.nx, config.ny
        )

    return lx, lambda_z, i_d, k_d


def sl_init(
    vx_e: np.ndarray,
    vz_e: np.ndarray,
    vx: np.ndarray,
    vz: np.ndarray,
    vx_p: np.ndarray,
    vz_p: np.ndarray,
    lsettls: bool = True,
) -> Tuple[np.ndarray]:
    """Initialize draft velocities with either
    LSETTLS method : 2 fields for velocity (at t and t - dt)
    LNESN (not LSETTLS) method : 1 field for velocity (at t)

    Args:
        vx_e (np.ndarray): outlined velocity at t + dt on x
        vy_e (np.ndarray): outlined velocity at t + dt on y
        vx (np.ndarray): velocity at t on x
        vy (np.ndarray): velocity at t on y
        vx_p (np.ndarray): velocity at t - dt on x
        vy_p (np.ndarray): velcoity at t - dt on y
        lsettls (bool, optional): LSETTLS or LNESC. Defaults to True.

    Returns:
        Tuple[np.ndarray]: velocities at t and t + dt
    """
    # LSETTLS
    if lsettls:
        vx_e = vx.copy()
        vz_e = vz.copy()

        vx = 2 * vx - vx_p
        vz = 2 * vz - vz_p

    # LNESC
    else:
        vx_e = vx.copy()
        vz_e = vz.copy()

    return vx, vz, vx_e, vz_e


def sl_xy(
    config: Config,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    tracer: np.ndarray,
    tracer_e: np.ndarray,
    interpolation_function: callable,
    nitmp: int,
) -> np.ndarray:
    """Performs tracer advection with 2D semi lagrangian.
    1: search for departure point
    2: interpolate tracer field

    Args:
        config (Config): grid configuration
        vx (np.ndarray): velocity on x
        vy (np.ndarray): velocity on y
        vx_e (np.ndarray): velocity on x at t + dt
        vy_e (np.ndarray): velocity on y at t + dt
        tracer (np.ndarray): tracer field
        tracer_e (np.ndarray): tracer at t + dt (ebauche)
        interpolation_function (callable): linear or cubic interpolation
        nitmp (int): number of iterations for departure search

    Returns:
        np.ndarray: tracer outline (ebauche) at t + dt
    """

    # Recherche semi lag
    lx_d, ly_d, i_d, j_d = lagrangian_search(
        config=config,
        vx_e=vx_e,
        vy_e=vy_e,
        vx=vx,
        vy=vy,
        interpolation_function=interpolate_lin_2d,
        nitmp=nitmp,
    )

    # Interpolate
    tracer_e = interpolation_function(
        tracer,
        lx_d,
        ly_d,
        i_d,
        j_d,
        config.bcx_kind,
        config.bcy_kind,
        config.nx,
        config.ny,
    )

    return tracer_e


def backup(
    vx: np.ndarray,
    vz: np.ndarray,
    vx_e: np.ndarray,
    vz_e: np.ndarray,
    tracer: np.ndarray,
    tracer_e: np.ndarray,
) -> Tuple[np.ndarray]:
    """Copy fields for next iteration.
    Ex : vx_e becomes vx at next model time step

    Args:
        vx (np.ndarray): x velocity
        vy (np.ndarray): y velocity
        vx_e (np.ndarray): ebauche vx
        vy_e (np.ndarray): ebauche vy
        tracer (np.ndarray): tracer field
        tracer_e (np.ndarray): ebauche at t + dt for tracer field

    Returns:
        Tuple[np.ndarray]: copy for fields
    """

    # Copie des champs
    tracer = tracer_e.copy()
    vx = vx_e.copy()
    vz = vz_e.copy()

    return vx, vz, tracer
