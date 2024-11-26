from typing import Tuple
import numpy as np
import logging

from config import Config
from tracer_lab_python.diagnostics import diagnostic_lipschitz
from tracer_lab_python.filter import overshoot_filter, undershoot_filter
from tracer_lab_python.interpolation import (
    interpolate_lin_2d,
    max_interpolator_2d,
    min_interpolator_2d,
)
from tracer_lab_python.periodic_filters import (
    periodic_overshoot_filter,
    periodic_undershoot_filter,
)

logging.getLogger(__name__)


def smilag_transport_scheme(
    config: Config,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    tracer: np.ndarray,
    tracer_e: np.ndarray,
    interpolation_function: callable,
    nitmp: int,
    filter: bool = False,
) -> np.ndarray:
    """Performs tracer advection with 2D semi lagrangian.
    1: search for departure point -> larcina
    2: interpolation -> larcinb
    
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

    #############################################
    ######### Departure point search ############
    #############################################
    # TODO : move in larcina
    lx_d, ly_d, i_d, j_d = larcina(
        config=config,
        vx_e=vx_e,
        vy_e=vy_e,
        vx=vx,
        vy=vy,
        interpolation_function=interpolate_lin_2d,
        nitmp=nitmp,
    )

    #############################################
    ######### Interpolation function ############
    #############################################
    # TODO : move in larcinb
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

    ##############################################
    ########## Overshoot filtering ###############
    ##############################################
    # TODO : check filter performances
    if config.filter:
        tracer_sup = max_interpolator_2d(
        tracer, i_d, j_d, config.bcx_kind, config.bcy_kind, config.nx, config.ny
        )

        tracer_inf = min_interpolator_2d(
        tracer, i_d, j_d, config.bcx_kind, config.bcy_kind, config.nx, config.ny
        )

        if config.bcx_kind == 1 and config.bcy_kind == 1:
            tracer_e = periodic_overshoot_filter(
                tracer_e, tracer_sup, config.nx, config.ny
            )
            tracer_e = periodic_undershoot_filter(
                tracer_e, tracer_inf, config.nx, config.ny
            )

        if config.bcx_kind == 0 and config.bcy_kind == 0:
            tracer_e = overshoot_filter(tracer_e, tracer_sup, config.nx, config.ny)
            tracer_e = undershoot_filter(tracer_e, tracer_inf, config.nx, config.ny)

    return tracer_e


def larcina(
    config: Config,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    elarmes_1d: callable,
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
    vy_tmp = vy.copy()
    
    # Spacings and time step
    dx, dy = config.dx, config.dy
    dth = config.dth
    
    # Indexes of arrival points
    i_arr, j_arr = config.I, config.J
    
    # Borders 
    bcx_kind, bcy_kind = config.bcx_kind, config.bcy_kind
    
    # Spacings
    nx, ny = config.nx, config.ny
    
    # Array declaration
    for l in range(nitmp):
        
        ##### ELARCHE
        trajx = elarche_1d(vx_e, vx_tmp, dx, dth)
        trajy = elarche_1d(vy_e, vy_tmp, dy, dth)
        
        ##### ELASCAW
        lx, ix = elascaw_1d(trajx, i_arr)
        ly, jy = elascaw_1d(trajy, j_arr)
        
        ##### ELARMES        
        vx_tmp = elarmes_1d(
            vx, lx, ly, ix, jy, bcx_kind, bcy_kind, nx, ny
        )
        vy_tmp = elarmes_1d(
            vy, lx, ly, ix, jy, bcx_kind, bcy_kind, nx, ny
        )

    return lx, ly, ix, jy

def elarche_1d(
    vhat_arr: np.ndarray,
    vhat_dep: np.ndarray,
    traj: np.ndarray,
    dth: np.float64,
    dx: np.float64
) -> np.ndarray:
    # Deplacement
    traj = -dth * (vhat_arr + vhat_dep) / dx
    return traj

def elascaw_1d(
    traj: np.ndarray,
    i_arr: np.ndarray,
) -> Tuple[np.ndarray]:
    
    ix = (i_arr + np.floor(traj)).astype(int)
    lx = traj - np.floor(traj)
    return ix, lx


def sl_init(
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_p: np.ndarray,
    vy_p: np.ndarray,
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
        vy_e = vy.copy()

        vx = 2 * vx - vx_p
        vy = 2 * vy - vy_p

    # LNESC
    else:
        vx_e = vx.copy()
        vy_e = vy.copy()

    return vx, vy, vx_e, vy_e


def backup(
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
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
    vy = vy_e.copy()

    return vx, vy, tracer
