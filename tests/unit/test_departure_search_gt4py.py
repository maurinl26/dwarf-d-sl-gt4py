import logging
import numpy as np
import sys
import yaml
import gt4py
from typing import Tuple

sys.path.append("/home/maurinl/sl_gt4py/src")
sys.path.append("/home/maurinl/sl_gt5py/tests")
print(sys.path)

from tracer_lab_gt4py.gt4py_config import dtype, dtype_int, backend, origin
from sandbox.sl_gt4py.departure_search import dep_search_1d
from utils.config import Config


logging.basicConfig(stream=sys.stdout, level=logging.INFO)


# Taken from Christian Kühnlein
def blossey_stf(xcr: np.ndarray, ycr: np.ndarray, mt: float) -> np.ndarray:
    """Compute blossey stream function as described by
    Durran & Blossey (Selective monotonicity preservation in scalar advection)

    Args:
        t (float): _description_
        T (float): _description_
        r (np.ndarray): _description_

    Returns:
        _type_: _description_
    """
    rsq = (xcr - 0.5) ** 2.0 + (ycr - 0.5) ** 2.0

    term1 = 0.5 * rsq
    term2 = np.log(1.0 - 16.0 * rsq + 256.0 * rsq**2.0) / 96.0
    term3 = -np.log(1.0 + 16.0 * rsq) / 48.0
    term4 = -np.sqrt(3.0) * np.arctan((-1.0 + 32.0 * rsq) / np.sqrt(3.0)) / 48.0

    stf = (
        4.0
        * np.pi
        * (term1 + np.cos(2.0 * np.pi * mt) * (term1 + term2 + term3 + term4))
    )

    return stf


# Taken from FVM Slim (Christian)
def stream_function_xy(stf: np.ndarray, dx, dy) -> Tuple[np.ndarray]:
    """Computes velocity components based on stream function

    Args:
        stf (np.ndarray): stream function

    Returns:
        Tuple[np.ndarray]: u (velocity on x axis), v (velocity on y axis)
    """
    
    # TODO : implementation du gradient a la main
    dstf_dx = np.gradient(stf, dx, axis=0)
    dstf_dy = np.gradient(stf, dy, axis=1)
    u = dstf_dy
    v = - dstf_dx

    return u, v


def blossey_velocity(xcr: np.ndarray, ycr: np.ndarray, mt, dx, dy) -> Tuple[np.ndarray]:
    stf = blossey_stf(xcr, ycr, mt)
    u, v = stream_function_xy(stf, dx, dy)

    return u, v

def init_blossey(xcr: np.ndarray, ycr: np.ndarray, mt, dt, dx, dy, nx: int, ny: int):
    # Vitesses
    vx, vy = blossey_velocity(xcr, ycr, mt, dx, dy)
    vx_p, vy_p = blossey_velocity(xcr, ycr, mt - dt, dx, dy), vy.copy()
    vx_e, vy_e = np.empty((nx, ny)), np.empty((nx, ny))

    return vx, vy, vx_p, vy_p, vx_e, vy_e

if __name__ == "__main__":
    # Shift in config file
    
    config_file = "./config/durran_blossey.yml"
    with open(config_file, 'r') as file:
        conf_dict = yaml.safe_load(file)
        
        config = Config(**conf_dict)
    
    # LSETTLS  
    lsettls = True
    nz = 1
    
    # Pour info
    T = config.model_endtime - config.model_starttime
    t = config.model_starttime
    nstep = np.ceil((config.model_endtime - config.model_starttime) / config.dt)

    logging.info(f"Time step dt : {config.dt:.06f} s")
    logging.info(f"N steps : {nstep:.06f} s")
    
    # Indices 
    aligned_index = (0, 0, 0)
    grid_indices = gt4py.storage.zeros((config.nx, config.ny, nz), dtype_int, backend=backend, aligned_index=aligned_index)
    departure_indices = gt4py.storage.zeros((config.nx, config.ny, nz), dtype, backend=backend, aligned_index=aligned_index)
    grid_indices[:, :] = config.I[:, :, np.newaxis]
    departure_indices[:, :] = config.I[:, :, np.newaxis]
    
    # Velocities
    vx = gt4py.storage.zeros((config.nx, config.ny, nz), dtype, backend=backend, aligned_index=aligned_index)
    vx_e = gt4py.storage.zeros((config.nx, config.nx, nz), dtype, backend=backend, aligned_index=aligned_index)
    
    # Init velocities (numpy -> gt4py)
    init_state = init_blossey(
         config.xcr, config.ycr, t, config.dt, config.dx, config.dy, config.nx, config.ny
    )
    vx[:, :] = init_state[0][:, :, np.newaxis]
    vx_e[:, :] = init_state[4][:, :, np.newaxis]
    
    # lambda
    lx = gt4py.storage.zeros((config.nx, config.ny, nz), dtype, backend=backend, aligned_index=aligned_index)
    lx[:, :] = 0
    
    # Temporary fields 
    vx_tmp = gt4py.storage.zeros((config.nx, config.ny, nz), dtype, backend=backend, aligned_index=aligned_index)

    dep_search_1d(
        grid_indices, 
        vx_e,
        vx_tmp,
        lx, 
        departure_indices, 
        config.dx,
        config.dth
    )