from typing import Tuple
import numpy as np
import logging

from gt4py.cartesian.gtscript import stencil, Field, function, floor
from gt4py_config import backend, backend_opts

from config import Config
from tracer_lab_gt4py.interpolation import elarmes_1d_linear
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
    lx_d, ly_d, i_d, j_d = larcina(
        config=config,
        vx_e=vx_e,
        vy_e=vy_e,
        vx=vx,
        vy=vy,
        nitmp=nitmp,
    )

    #############################################
    ######### Interpolation function ############
    #############################################
    # Not in GT4PY
    tracer_e = larcinb(
        tracer,
        lx_d,
        ly_d,
        i_d,
        j_d,
        config.bcx_kind,
        config.bcy_kind,
        config.nx,
        config.ny,
        interpolation_function
    )

    ##############################################
    ########## Overshoot filtering ###############
    ##############################################
    # Not in GT4PY
    if config.filter:
        tracer_e = filter_driver(
            tracer, i_d, j_d, config.bcx_kind, config.bcy_kind, config.nx, config.ny
        )
  
    return tracer_e

##################################################
########### LARCINA GT4PY + DaCe #################
##################################################
def larcina(
    config: Config,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
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
    idx_arr, jdx_arr = config.I, config.J
    
    # Borders 
    bcx_kind, bcy_kind = config.bcx_kind, config.bcy_kind
    
    # Spacings
    nx, ny = config.nx, config.ny
    
    ### tmps 
    # lx, ly
    # idx_dep, jdy_dep
    
    # Array declaration
    # TODO : in GT4PY + DaCe
    for l in range(nitmp):
        
        ##### ELARCHE + ELASCAW
        slag_search(
            weight_x_dep=lx,
            weight_y_dep=ly,
            idx_dep=idx_dep,
            jdy_dep=idy_dep,
            vx_tmp=vx_tmp,
            vy_tmp=vy_tmp,
            vx_e=vx_f,
            vy_e=vy_f,
            idx_arr=idx_arr,
            idx_dep=idx_dep,
            dx=dx,
            dy=dy,
            dth=dth
        )
        
        ##### ELARMES
        elarmes_1d_linear(
            tracer=vx_tmp, 
            tracer_e=vx_e,
            weight_x=weight_x,
            weight_y=weight_y,
            idx_dep=idx_dep,
            idy_dep=idx_dep 
            ) # for vx
        elarmes_1d_linear(
            tracer=vy_tmp, 
            tracer_e=vy_e,
            weight_x=weight_x,
            weight_y=weight_y,
            idx_dep=idx_dep,
            idy_dep=idx_dep 
            ) # for vy 
   
    return lx, ly, ix, jy

@stencil(backend, **backend_opts)
def slag_search(
    weight_x_dep: Field[np.float64], 
    weight_y_dep: Field[np.float64], 
    idx_dep: Field[np.int64],
    idy_dep: Field[np.int64],
    vx_e: Field[np.float64],
    vx_tmp: Field[np.float64],
    vy_e: Field[np.float64],
    vy_tmp: Field[np.float64],
    idx_arr: Field[np.int64],
    idy_arr: Field[np.int64],
    dx: np.float64,
    dy: np.float64,
    dth: np.float64
):
    
    with computation(PARALLEL), interval(...):
        
        ##### ELARCHE
        trajx = elarche_1d(vx_e, vx_tmp, dx, dth)
        trajy = elarche_1d(vy_e, vy_tmp, dy, dth)
        
        ##### ELASCAW
        weight_x_dep = elascaw_1d_weight(trajx)
        idx_dep = elascaw_1d_index(trajx, idx_arr)
        weight_y_dep = elascaw_1d_weight(trajy)
        idy_dep = elascaw_1d_index(trajy, idy_arr)

@function
def elarche_1d(
    vhat_arr: Field[np.float64],
    vhat_dep: Field[np.float64],
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

@function
def elascaw_1d_weight(
    traj: Field[np.float64],
):
    return traj - floor(traj)

@function
def elascaw_1d_index(
    traj: Field[np.float64],
    idx_arr: Field[np.int64],
):
    return idx_arr + floor(traj)


#################################################
#################### LARCINB ####################
#################################################

def larcinb(
    tracer: np.ndarray,
    weight_x: np.ndarray,
    weight_y: np.ndarray,
    dep_idx_x: np.ndarray,
    dep_idx_y: np.ndarray,
    bcx_kind,
    bcy_kind,
    nx, 
    ny,
    interpolation_function: callable
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
    
    tracer_e = interpolation_function(
        tracer,
        weight_x,
        weight_y,
        dep_idx_x,
        dep_idx_y,
        bcx_kind,
        bcy_kind,
        nx,
        ny,
    )
    
    return tracer_e

