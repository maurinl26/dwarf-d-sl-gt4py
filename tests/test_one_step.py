import logging
import time
from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("/home/maurinl/sl_gt4py/src")
print(sys.path)

from sl_python.interpolation import interpolate_cub_2d, interpolate_lin_2d
from sl_python.sl_2D import sl_init, sl_xy
from config import Config
from utils.cfl import cfl_1d
logging.basicConfig(stream=sys.stdout, level=logging.INFO)


# Pointwise functions
def r_tilde(x, y):
    return 5 * np.sqrt((x - 0.5) ** 2 + (y - 0.5) ** 2)


def tracer_shape(x: float, y: float, tracer0: float):
    r_t = r_tilde(x, y)
    return np.where(r_t < 1, tracer0 + ((1 + np.cos(np.pi * r_t)) / 2) ** 2, tracer0)


def init_uniform(U: float, V: float, nx: int, ny: int):
    # Vitesses
    vx, vy = U * np.ones((nx, ny)), V * np.ones((nx, ny))
    vx_p, vy_p = vx.copy(), vy.copy()
    vx_e, vy_e = np.empty((nx, ny)), np.empty((nx, ny))

    return vx, vy, vx_p, vy_p, vx_e, vy_e


def blossey_tracer(xcr, ycr):
    # Tracer
    tracer = tracer_shape(xcr, ycr, 0)
    tracer_e = np.zeros(xcr.shape)
    return tracer, tracer_e


def backup(tracer, tracer_e):
    tracer = tracer_e.copy()
    return tracer


# Driver
def one_step_driver(
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
    interp_function: callable = interpolate_lin_2d
):
    # Advection
    vx_e, vy_e, vx, vy = sl_init(
        vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, lsettls=lsettls
    )
        
    logging.info(f"Backup")
    logging.info(f"Tracer : min {tracer.min():.02f}, max {tracer.max():.02f}, mean {tracer.mean():.02f}")
    logging.info(f"Tracer e : min {tracer_e.min():.02f}, max {tracer_e.max():.02f}, mean {tracer_e.mean():.02f}")

    # Estimations
    tracer_e = sl_xy(
        config=config,
        vx=vx,
        vy=vy,
        vx_e=vx_e,
        vy_e=vy_e,
        tracer=tracer,
        tracer_e=tracer_e,
        interpolation_function=interp_function,
        nitmp=4,
        )

    # Diagnostics and outputs
    courant_xmax = np.max(cfl_1d(vx_e, config.dx, config.dt))
    courant_ymax = np.max(cfl_1d(vy_e, config.dy, config.dt))

    logging.info(f"Maximum courant number : {max(courant_xmax, courant_ymax):.02f}")
    logging.info(f"Tracer (t + dt) : min {tracer_e.min():.02f}, max {tracer_e.max():.02f}, mean {tracer_e.mean():.02f}")

    plot_blossey(config.xcr, config.ycr, vx, vy, tracer_e, f"./figures/one_step/one_step_end.pdf")
    
    return tracer_e


def plot_blossey(
    xcr: np.ndarray,
    ycr: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    tracer: np.ndarray,
    output_file: str,
):
    # Shape
    fig, ax = plt.subplots()

    # Vent
    ax.quiver(
        xcr[::2, ::2], ycr[::2, ::2], vx[::2, ::2], vy[::2, ::2],
        color="C0",
        angles="xy",
        scale_units="xy",
        scale=5,
        width=0.002,
    )
    ax.set(xlim=(0, 1), ylim=(0, 1))

    levels = [0.05 + i * 0.1 for i in range(0, 10)]
    ax.contour(xcr, ycr, tracer, colors="black", vmin=0.05, vmax=0.95, levels=levels)

    plt.savefig(output_file)


if __name__ == "__main__":
    # Shift in config file
    T = 1
    t = 0

    # TODO : send config to config field (.yml)
    model_starttime = 0
    model_endtime = 1
    nstep = 1
    nitmp = 4
    dt = (model_endtime - model_starttime) / nstep
    xmin, xmax = 0, 1
    ymin, ymax = 0, 1
    nx, ny = 50, 50
    U, V = 0.2, 0.2
    lsettls = True
    bcx_kind, bcy_kind = 1, 1


    # TODO : initialize config
    config = Config(dt, xmin, xmax, nx, ymin, ymax, ny, bcx_kind, bcy_kind)
    
    logging.info("Config")
    logging.info(f"time step dt : {config.dt} s")
    logging.info(f"dx : {config.dx}, dy : {config.dy}")
    logging.info(f"Uniform velocity U = {U}, V = {V}")

    vx, vy, vx_p, vy_p, vx_e, vy_e = init_uniform(
        U, V, nx, ny
    )
    tracer, tracer_e = blossey_tracer(config.xcr, config.ycr)
    tracer_ref = tracer.copy()
    
    logging.info(f"Tracer : min {tracer.min():.02f}, max {tracer.max():.02f}, mean {tracer.mean():.02f}")
    plot_blossey(
        config.xcr, config.ycr, vx, vy, tracer, f"./figures/one_step/one_step_start.pdf"
    )
    
    #### Interp Lineaire

    # Advection encapsulation
    start_time = time.time()
    tracer  = one_step_driver(config, vx, vy, vx_e, vy_e, vx_p, vy_p, tracer, tracer_e, lsettls, interpolate_lin_2d)
    duration = time.time() - start_time
    logging.info(f"Duration : {duration:.02f} s")
    
    # Errors 
    e_inf = np.max(np.abs(tracer - tracer_ref))
    e_2 = np.sqrt((1 / (config.nx * config.ny)) * np.sum((tracer - tracer_ref) ** 2))

    logging.info(f"Error E_inf : {e_inf}")
    logging.info(f"Error E_2 : {e_2}")
    
    #### Interp Cubique
    
    # Advection encapsulation
    
    tracer, tracer_e = blossey_tracer(config.xcr, config.ycr)
    tracer_ref = tracer.copy()
    
    start_time = time.time()
    tracer  = one_step_driver(config, vx, vy, vx_e, vy_e, vx_p, vy_p, tracer, tracer_e, lsettls, interpolate_cub_2d)
    duration = time.time() - start_time
    logging.info(f"Duration : {duration:.02f} s")
    
    # Errors 
    e_inf = np.max(np.abs(tracer - tracer_ref))
    e_2 = np.sqrt((1 / (config.nx * config.ny)) * np.sum((tracer - tracer_ref) ** 2))

    logging.info(f"Error E_inf : {e_inf}")
    logging.info(f"Error E_2 : {e_2}")