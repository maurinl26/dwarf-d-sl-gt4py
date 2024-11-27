from typing import Tuple
import numpy as np
import logging

from gt4py.cartesian.gtscript import stencil, Field, function, floor
from tracer_lab_gt4py.gt4py_config import backend, backend_opts, origin
from gt4py.storage import from_array, empty

from utils.config import Config
from tracer_lab_python.filter import filter_driver
from dace import program

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
    if filter:
        tracer_e = filter_driver(
            tracer, i_d, j_d, config.bcx_kind, config.bcy_kind, config.nx, config.ny
        )
  
    return tracer_e

##################################################
########### LARCINA GT4PY + DaCe #################
##################################################
@program
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

    # vx should be in reading
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
    nx, ny, nz = config.nx, config.ny, config.nz
    
    ############################################
    ######## numpy to gt4py storage ############
    ############################################
    # In arrays
    vx_e_in = from_array(vx_e, np.float64, backend=backend, aligned_index=origin)
    vy_e_in = from_array(vy_e, np.float64, backend=backend, aligned_index=origin)
        
    # Temporary arrays
    weight_x_dep_out = empty((nx, ny, nz), np.float64, backend=backend, aligned_index=origin)
    weight_y_dep_out = empty((nx, ny, nz), np.float64, backend=backend, aligned_index=origin)
    idx_x_dep_out = empty((nx, ny, nz), np.int64, backend=backend, aligned_index=origin)
    idx_y_dep_out = empty((nx, ny, nz), np.int64, backend=backend, aligned_index=origin)
    idx_x_arr_in = empty((nx, ny, nz), np.int64, backend=backend, aligned_index=origin)
    idx_y_arr_in = empty((nx, ny, nz), np.int64, backend=backend, aligned_index=origin)
    
    
    # Array declaration
    # TODO : in GT4PY + DaCe
    for l in range(nitmp):
        
        # Updated fields
        vx_tmp_in = from_array(vx_tmp, np.float64, backend=backend, aligned_index=origin)
        vy_tmp_in = from_array(vy_tmp, np.float64, backend=backend, aligned_index=origin)
        
        ###################################
        ##### ELARCHE + ELASCAW GT4Py #####
        ###################################
        slag_search(
            weight_x_dep=weight_x_dep_out,
            weight_y_dep=weight_y_dep_out,
            idx_dep=idx_x_dep_out,
            jdy_dep=idx_y_dep_out,
            vx_tmp=vx_tmp_in,
            vy_tmp=vy_tmp_in,
            vx_e=vx_e_in,
            vy_e=vy_e_in,
            idx_arr=idx_x_arr_in,
            idy_arr=idx_y_arr_in,
            dx=dx,
            dy=dy,
            dth=dth
        )
        
        
        ######################################
        ########## GT4Py to Numpy ############
        ######################################
        lx_dep = np.asarray(weight_x_dep_out)
        ly_dep = np.asarray(weight_y_dep_out)
        idx_dep = np.asarray(idx_x_dep_out)
        idy_dep = np.asarray(idx_y_dep_out)
        
        ######################################
        ############ ELARMES #################
        ######################################
        # TODO: elarmes in gt4py
        vx_tmp = elarmes_1d(
            vx, lx_dep, ly_dep, idx_dep, idy_dep, bcx_kind, bcy_kind, nx, ny
        )
        vy_tmp = elarmes_1d(
            vy, lx_dep, ly_dep, idx_dep, idy_dep, bcx_kind, bcy_kind, nx, ny
        )
   
    return lx_dep, ly_dep, idx_dep, idy_dep

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


#################################################
############## ELARMES DaCe #####################
#################################################
# TODO : elarmes in DaCe
def elarmes_1d(
    psi: np.ndarray,
    lx: np.ndarray,
    ly: np.ndarray,
    i_d: np.ndarray,
    j_d: np.ndarray,
    bcx_kind: int,
    bcy_kind: int,
    nx: int,
    ny: int,
):
    """Perform a 1d linear interpolation

    Args:
        lx (np.ndarray): _description_
        ly (np.ndarray): _description_
        psi (np.ndarray): _description_
        i_d (np.ndarray): _description_
        j_d (np.ndarray): _description_
        bcx_kind (int): _description_
        bcy_kind (int): _description_
        nx (int): _description_
        ny (int): _description_
    """
    # Polynômes de lagrange
    p0 = lambda l: 1 - l
    p1 = lambda l: l

    # Interp selon x -> field_hat_x
    px = np.array([p0(lx), p1(lx)])
    py = np.array([p0(ly), p1(ly)])

    # 1. Construire les tableaux d'indices i_d0, i_d1 / j_d0, j_d1
    # Non periodique
    id_0 = boundaries(i_d, nx, bcx_kind)
    id_p1 = boundaries(i_d + 1, nx, bcx_kind)
    
    jd_0 = boundaries(j_d, ny, bcy_kind)
    jd_p1 = boundaries(j_d + 1, ny, bcy_kind)

    # Lookup
    psi_d_i = np.zeros((4, nx, ny))
    for i in range(nx):
        for j in range(ny):
            psi_d_i[0, i, j] = psi[id_0[i, j], jd_0[i, j]]
            psi_d_i[1, i, j] = psi[id_p1[i, j], jd_0[i, j]]

            psi_d_i[2, i, j] = psi[id_0[i, j], jd_p1[i, j]]
            psi_d_i[3, i, j] = psi[id_p1[i, j], jd_p1[i, j]]

    psi_d_j = np.zeros((2, nx, ny))
    psi_d_j[0] = px[0] * psi_d_i[0] + px[1] * psi_d_i[1]
    psi_d_j[1] = px[0] * psi_d_i[2] + px[1] * psi_d_i[3]

    psi_d = py[0] * psi_d_j[0] + py[1] * psi_d_j[1]

    return psi_d
    
    


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

