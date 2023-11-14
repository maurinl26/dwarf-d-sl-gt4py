from typing import Tuple
import numpy as np
import logging
import sys

from config import Config

logging.getLogger(__name__)


def departure_search(
    xcr: np.ndarray,
    ycr: np.ndarray, 
    I: np.ndarray,
    J: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    vx_tmp: np.ndarray,
    vy_tmp: np.ndarray,
    dth: float,
    dx: float,
    dy: float,
    epsilon: float
) -> Tuple[np.ndarray]:
    """Compute departure points coordinates, with respect 
    to nearest grid points. 

    Args:
        xcr (np.ndarray): x coordinates of gridpoints
        ycr (np.ndarray): y coordinates of gridpoints
        I (np.ndarray): indices of grid points (x direction)
        J (np.ndarray): indices of grid points (y direction)
        vx_e (np.ndarray): velocity on x at t + dt (ebauche)
        vy_e (np.ndarray): velocity on y at t + dt (ebauche)
        vx_tmp (np.ndarray): velocity on x at t for updated departure point
        vy_tmp (np.ndarray): velocity on y at t for updated departure point
        dth (float): half time step
        dx (float): x spacing on grid
        dy (float): y spacing on grid
        epsilon (float): machine minimum

    Returns:
        Tuple[np.ndarray]: 
                lx, ly -> distance from departure point to grid point, 
                id, jd -> indices from  ref grid points (near from departure point)
    """
    
    x_d = xcr - dth * (vx_e + vx_tmp) 
    y_d = ycr - dth * (vy_e + vy_tmp)
        
    lx = (xcr - x_d) / dx
    ly = (ycr - y_d) / dy
        
    i_dep = np.where(abs(lx) < epsilon, 0, np.floor(lx)).astype(int)
    j_dep = np.where(abs(ly) < epsilon, 0, np.floor(ly)).astype(int)
    
    lx = - lx - (i_dep + 1)
    ly = - ly - (j_dep + 1)
     
    i_d = I - i_dep
    j_d = J - j_dep   
    
    return lx, ly, I - i_d, J - j_d
    

# ELARCHE
def lagrangian_search(
    config: Config,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
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
    vy_tmp = vy.copy()

    # Array declaration
    for l in range(nitmp):
        
        lx, ly, i_d, j_d = departure_search(
            config.xcr,
            config.ycr,
            config.I,
            config.J,
            vx_e,
            vy_e,
            vx_tmp,
            vy_tmp,
            config.dth,
            config.dx,
            config.dy,
            np.finfo(float).eps
        )

        ####### Interpolation for fields ########
        vx_tmp = interpolation_function(
            vx,
            lx, 
            ly,
            i_d, 
            j_d,
            config.bcx_kind,
            config.bcy_kind,
            config.nx, 
            config.ny
        )
        
        vy_tmp = interpolation_function(
            vy,
            lx, 
            ly,
            i_d, 
            j_d,
            config.bcx_kind,
            config.bcy_kind,
            config.nx, 
            config.ny
        )

    return lx, ly, i_d, j_d


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
    
    logging.info(f"(sl_xy) Tracer : mean {tracer.mean()}")
    
    # Recherche semi lag
    lx_d, ly_d, i_d, j_d = lagrangian_search(
        config=config,
        vx_e=vx_e,
        vy_e=vy_e,
        vx=vx,
        vy=vy,
        interpolation_function=interpolation_function,
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
        config.ny
    )
    
    logging.info(f"(sl_xy) Tracer (t + dt) : mean {tracer_e.mean()}")
    
    return tracer_e


def backup(
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    tracer: np.ndarray,
    tracer_e: np.ndarray,
) -> Tuple[np.ndarray]:
    """Remap ebauche fields at t into fields at tau = t + dt

    Args:
        vx (np.ndarray): x velocity
        vy (np.ndarray): y velocity
        vx_e (np.ndarray): ebauche vx
        vy_e (np.ndarray): ebauche vy
        tracer (np.ndarray): tracer field
        tracer_e (np.ndarray): ebauche at t + dt for tracer field

    Returns:
        _type_: _description_
    """

    # Copie des champs
    tracer = tracer_e.copy()
    vx = vx_e.copy()
    vy = vy_e.copy()

    return vx, vy, tracer
