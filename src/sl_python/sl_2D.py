from typing import Tuple
import numpy as np
import itertools

from config import Config
from sl_python.interpolation import interpolate_lin_2d

def boundaries(
    bc_kind: int,
    indices: np.ndarray,
    n: int,
):
    """Apply boundary conditions on field.

    Args:
        bc_kind (int): _description_
        field (np.ndarray): _description_
        min (np.float): _description_
        max (np.float): _description_
        nx (int): _description_
    """
    left_exceed = indices >= n
    right_exceed = indices < 0

    # Periodic boundaries
    if bool(bc_kind):
        indices %= n

    # Fixed boundaries
    else:
        indices = (
            indices * (1 - right_exceed) * (1 - left_exceed)
            + 0 * left_exceed
            + (n - 1) * right_exceed
        )

    return indices


def departure_search(
    xcr: np.ndarray,
    ycr: np.ndarray, 
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    vx_tmp: np.ndarray,
    vy_tmp: np.ndarray,
    dth: float,
    dx: float,
    dy: float,
    bcx_kind: int,
    bcy_kind: int,
    nx: int, 
    ny: int,
    epsilon: float
):
    
    x_d = xcr - dth * (vx_e + vx_tmp) 
    y_d = ycr - dth * (vy_e + vy_tmp)
        
    lx = (xcr - x_d) / dx
    ly = (ycr - y_d) / dy
        
    i_d = np.where(abs(lx) < epsilon, 0, np.floor(lx))
    j_d = np.where(abs(ly) < epsilon, 0, np.floor(ly))
    
    lx = -lx - (i_d + 1)
    ly = -ly - (j_d + 1)
    
    i_d = boundaries(bcx_kind, i_d, nx)
    j_d = boundaries(bcy_kind, j_d, ny)
    
    return lx, ly, i_d, j_d
    

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
    for l in range(0, nitmp):
        
        lx, ly, i_d, j_d = departure_search(
            config.xcr,
            config.ycr,
            vx_e,
            vy_e,
            vx_tmp,
            vy_tmp,
            config.dth,
            config.bcx_kind,
            config.bcy_kind,
            config.nx,
            config.ny,
            np.finfo(float).eps
        )

        ####### Interpolation for fields ########
        vx_tmp = interpolate_lin_2d(
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
        
        vy_tmp = interpolate_lin_2d(
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

def sl_xy_tracer(
    config: Config,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    tracer: np.ndarray,
    tracer_e: np.ndarray,
    interpolation_function: callable,
    nsiter: int,
):
    """Interpolate only tracers

    Args:
        config (Config): _description_
        vx (np.ndarray): _description_
        vy (np.ndarray): _description_
        vx_e (np.ndarray): _description_
        vy_e (np.ndarray): _description_
        tracer (np.ndarray): _description_
        tracer_e (np.ndarray): _description_
        interpolation_function (callable): _description_
        nsiter (int): _description_
    """
    
    # Recherche semi lag
    lx_d, ly_d, i_d, j_d = lagrangian_search(
        config=config,
        vx_e=vx_e,
        vy_e=vy_e,
        vx=vx,
        vy=vy,
        interpolation_function=interpolation_function,
        nsiter=nsiter,
    )

    # Interpolate
    tracer_e = interpolate_lin_2d(
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
    
    return tracer_e


def sl_xy(
    config: Config,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    tracer: np.ndarray,
    tracer_e: np.ndarray,
    interpolation_function: callable,
    nsiter: int,
):
    # Recherche semi lag
    lx_d, ly_d, i_d, j_d = lagrangian_search(
        config=config,
        vx_e=vx_e,
        vy_e=vy_e,
        vx=vx,
        vy=vy,
        interpolation_function=interpolation_function,
        nsiter=nsiter,
    )
    
    # Interpolation
    tracer_e = interpolate_lin_2d(
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
    vx_e = interpolate_lin_2d(
        vx,
        lx_d, 
        ly_d,
        i_d, 
        j_d,
        config.bcx_kind,
        config.bcy_kind,
        config.nx, 
        config.ny
    )
    vy_e = interpolate_lin_2d(
        vy,
        lx_d, 
        ly_d,
        i_d, 
        j_d,
        config.bcx_kind,
        config.bcy_kind,
        config.nx, 
        config.ny
    )

    return vx_e, vy_e, tracer_e


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
