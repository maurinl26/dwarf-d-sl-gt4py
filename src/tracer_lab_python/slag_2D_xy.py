from typing import Tuple
import numpy as np
import logging

from config import Config
from tracer_lab_python.boundaries import boundaries
from tracer_lab_python.filter import filter_driver
from tracer_lab_python.interpolation import interpolate_lin_2d

logging.getLogger(__name__)


def smilag_transport_scheme(
    config: Config,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    tracer: np.ndarray,
    nitmp: int,
    interpolation_function: callable = interpolate_lin_2d,
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
    dep_weight_x, dep_weight_y, dep_idx_x, dep_idx_y = larcina(
        config=config,
        vx_e=vx_e,
        vy_e=vy_e,
        vx=vx,
        vy=vy,
        elarmes_1d=interpolate_lin_2d,
        nitmp=nitmp,
    )

    #############################################
    ######### Interpolation function ############
    #############################################
    tracer_e = larcinb(
        tracer=tracer,
        weight_x=dep_weight_x,
        weight_y=dep_weight_y,
        dep_idx_x=dep_idx_x,
        dep_idx_y=dep_idx_y,
        config=config
    )

    ##############################################
    ########## Overshoot filtering ###############
    ##############################################
    # TODO : check filter performances
    if config.filter:
        tracer_e = filter_driver(
            tracer, tracer_e, dep_idx_x, dep_idx_y, config.bcx_kind, config.bcy_kind, config.nx, config.ny
        )
  
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
        ix, lx = elascaw_1d(trajx, i_arr)
        jy, ly = elascaw_1d(trajy, j_arr)
                
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
    dth: np.float64,
    dx: np.float64
) -> np.ndarray:
    """Computes the trajectory from departure point to
    arrival point.
    
    The trajectory is basically the CFL field.

    Args:
        vhat_arr (np.ndarray): arrival point
        vhat_dep (np.ndarray): departure point
        traj (np.ndarray): trajectory
        dth (np.float64): half time step
        dx (np.float64): spacing

    Returns:
        np.ndarray: _description_
    """
    
    
    # Deplacement
    traj = -dth * (vhat_arr + vhat_dep) / dx
    return traj

def elascaw_1d(
    traj: np.ndarray,
    i_arr: np.ndarray,
) -> Tuple[np.ndarray]:
    """Returns departure index, and spatial weight given
    a trajectory and a arrival part.
    
    The departure index is the arrival index minus the floor part
    of the trajectory (CFL).
    
    The weight is the decimal part of the CFL.
    
    Args:
        traj (np.ndarray): traject
        i_arr (np.ndarray): arrival point

    Returns:
        Tuple[np.ndarray]: index and weight of departure point
    """
    
    ix = (i_arr + np.floor(traj)).astype(int)
    lx = traj - np.floor(traj)
    return ix, lx

def larcinb(
    tracer: np.ndarray,
    weight_x: np.ndarray,
    weight_y: np.ndarray,
    dep_idx_x: np.ndarray,
    dep_idx_y: np.ndarray,
    config: Config,
) -> np.ndarray:
    """Perform interpolation of a tracer field at departure points.
    

    Args:
        tracer (np.ndarray): tracer field to interpolate
        weight_x (np.ndarray): _description_
        weight_y (np.ndarray): _description_
        dep_idx_x (np.ndarray): index of departure point (on grid)
        dep_idx_y (np.ndarray): index of departure point (on grid)
        bcx_kind (_type_): _description_
        bcy_kind (_type_): _description_
        nx (_type_): _description_
        ny (_type_): _description_
        interpolation_function (callable): interpolation methode

    Returns:
        np.ndarray: _description_
    """
    
    bcx_kind, bcy_kind = config.bcx_kind, config.bcy_kind
    nx, ny = config.nx, config.ny
    
    tracer_e = interpolate_lin_2d(
        psi=tracer,
        lx=weight_x,
        ly=weight_y,
        i_d=dep_idx_x,
        j_d=dep_idx_y,
        bcx_kind=bcx_kind,
        bcy_kind=bcy_kind,
        nx=nx,
        ny=ny
    )
    
    return tracer_e

