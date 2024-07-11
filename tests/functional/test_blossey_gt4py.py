import logging
import time
from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt
import sys
import yaml
import gt4py


from sl_gt4py.gt4py_config import dtype, backend, origin, backend_opts
from sl_python_numba.interpolation import interpolate_cub_2d
from sl_python.sl_2D import sl_xy, sl_init
from config import Config
from .test_blossey import plot_blossey, plot_tracer_against_reference

logging.basicConfig(stream=sys.stdout, level=logging.INFO)


# Pointwise functions
def r(x, y):
    return np.sqrt((x - 0.5) ** 2 + (y - 0.5) ** 2)


def r_tilde(x, y):
    return 5 * np.sqrt((x - 0.3) ** 2 + (y - 0.5) ** 2)


def theta(x, y):
    return np.arctan((y - 0.5) / (x - 0.5))


def u(r, t, T):
    return (4 * np.pi * r / T) * (
        1 + np.cos(2 * np.pi * t / T) * (1 - (4 * r) ** 6) / (1 + (4 * r) ** 6)
    )


def cfl_1d(u: float, dx: float, dt: float):
    return u * dt / dx


def tracer_shape(x: float, y: float, tracer0: float):
    r_t = r_tilde(x, y)
    return np.where(r_t < 1, tracer0 + ((1 + np.cos(np.pi * r_t)) / 2) ** 2, tracer0)


# Field functions
def polar_coordinates(xcr: np.ndarray, ycr: np.ndarray):
    R = r(xcr, ycr)
    Theta = theta(xcr, ycr)
    return R, Theta


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


def blossey_tracer(xcr, ycr):
    # Tracer
    tracer = tracer_shape(xcr, ycr, 0)
    tracer_e = tracer.copy()
    return tracer, tracer_e


def backup(tracer, tracer_e):
    tracer = tracer_e.copy()
    return tracer


# Driver
def sl_driver(
    config: Config,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    vx_p: np.ndarray,
    vy_p: np.ndarray,
    tracer: np.ndarray,
    tracer_e: np.ndarray,
    lsettls: bool,
    model_starttime: float,
    model_endtime: float,
):
    tracer_ref = tracer.copy()

    # Advection
    # TODO: shift sl init to gt4py stencil
    vx_e, vy_e, vx, vy = sl_init(
        vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, lsettls=lsettls
    )

    t = model_starttime
    jstep = 0
    while t < model_endtime:
        
        jstep += 1
        t += config.dt
        logging.info(f"Step : {jstep}")
        logging.info(f"Time : {100*t/(model_endtime - model_starttime):.02f}%")

        # TODO: shift blossey init to gt4py stencil
        vx, vy = blossey_velocity(config.xcr, config.ycr, t, config.dx, config.dy)
        vx_e, vy_e = blossey_velocity(
            config.xcr, config.ycr, t + config.dt, config.dx, config.dy
        )

        # Estimations
        # TODO : shift first part of sl_xy into gt4py
        tracer_e = sl_xy(
            config=config,
            vx=vx,
            vy=vy,
            vx_e=vx_e,
            vy_e=vy_e,
            tracer=tracer,
            tracer_e=tracer_e,
            interpolation_function=interpolate_cub_2d,
            nitmp=4,
        )

        # TODO : backup in GT4Py
        tracer = backup(tracer=tracer, tracer_e=tracer_e)

        # Diagnostics and outputs
        # TODO : cfl1d in GT4Py
        courant_xmax = np.max(cfl_1d(vx_e, config.dx, config.dt))
        courant_ymax = np.max(cfl_1d(vy_e, config.dy, config.dt))

        logging.info(f"Maximum courant number : {max(courant_xmax, courant_ymax):.02f}")

if __name__ == "__main__":
    # Shift in config file
    
    config_file = "./config/durran_blossey.yml"
    with open(config_file, 'r') as file:
        conf_dict = yaml.safe_load(file)
        config = Config(**conf_dict)
    
    # LSETTLS  
    lsettls = True
    
    # Pour info
    T = config.model_endtime - config.model_starttime
    t = config.model_starttime
    nstep = np.ceil((config.model_endtime - config.model_starttime) / config.dt)

    logging.info(f"Time step dt : {config.dt:.06f} s")
    logging.info(f"N steps : {nstep:.06f} s")
    
    # Velocities
    vx = gt4py.storage((config.nx, config.ny), dtype, backend=backend, aligned_index=origin)
    vy = gt4py.storage((config.nx, config.ny), dtype, backend=backend, aligned_index=origin)
    vx_p = gt4py.storage((config.nx, config.ny), dtype, backend=backend, aligned_index=origin)
    vy_p = gt4py.storage((config.nx, config.ny), dtype, backend=backend, aligned_index=origin)
    vx_e = gt4py.storage((config.nx, config.ny), dtype, backend=backend, aligned_index=origin)
    vy_e = gt4py.storage((config.nx, config.ny), dtype, backend=backend, aligned_index=origin)
    
    # Init velocities (numpy -> gt4py)
    init_state = init_blossey(
         config.xcr, config.ycr, t, config.dt, config.dx, config.dy, config.nx, config.ny
    )
    vx[:, :] = init_state[0]
    vy[:, :] = init_state[1]
    vx_p[:, :] = init_state[2]
    vy_p[:, :] = init_state[3]
    vx_e[:, :] = init_state[4]
    vy_e[:, :] = init_state[5]    
    
    # Tracer
    tracer = gt4py.storage((config.nx, config.ny), dtype, backend=backend, aligned_index=origin)
    tracer_e = gt4py.storage((config.nx, config.ny), dtype, backend=backend, aligned_index=origin)
    
    tracer_state = blossey_tracer(config.xcr, config.ycr)
    tracer[:, :] = tracer_state[0]
    tracer_e[:, :] = tracer_state[1]

    # Advection encapsulation
    start_time = time.time()
    sl_driver(
        config,
        vx,
        vy,
        vx_e,
        vy_e,
        vx_p,
        vy_p,
        tracer,
        tracer_e,
        lsettls,
        config.model_starttime,
        config.model_endtime,
    )
    duration = time.time() - start_time
    logging.info(f"Duration : {duration} s")
    
    # Pickup fields
    out_tracer = np.asarray(tracer)
    out_vx = np.asarray(vx)
    out_vy = np.asarray(vy)
    
    e_inf = np.max(np.abs(tracer - tracer_state[0]))
    e_2 = np.sqrt((1 / (config.nx * config.ny)) * np.sum((tracer - tracer_state[0]) ** 2))

    plot_blossey(config.xcr, config.ycr, out_vx, out_vy, out_tracer, 1)
    plot_tracer_against_reference(config.xcr, config.yxr, out_tracer, tracer_state[0], e_2, e_inf)
