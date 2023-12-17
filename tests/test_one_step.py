import logging
import time
from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt
import sys

import yaml

from sl_python.plot import plot_blossey


from sl_python.interpolation import interpolate_cub_2d, interpolate_lin_2d
from sl_python.sl_2D import sl_init, sl_xy
from test_blossey import tracer_shape, blossey_tracer, r_tilde
from config import Config
from utils.cfl import cfl_1d

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logging.getLogger(__name__)

def init_uniform(U: float, V: float, nx: int, ny: int):
    # Vitesses
    vx, vy = U * np.ones((nx, ny)), V * np.ones((nx, ny))
    vx_p, vy_p = vx.copy(), vy.copy()
    vx_e, vy_e = np.empty((nx, ny)), np.empty((nx, ny))

    return vx, vy, vx_p, vy_p, vx_e, vy_e


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

if __name__ == "__main__":
    
    
    config_file = "./config/uniform_advection.yml"
    with open(config_file, 'r') as file:
        conf_dict = yaml.safe_load(file)
        config = Config(**conf_dict)

    # TODO : send config to config field (.yml)
    U, V = 0.2, 0.2
    lsettls = True

    logging.info("Config")
    logging.info(f"time step dt : {config.dt} s")
    logging.info(f"dx : {config.dx}, dy : {config.dy}")
    logging.info(f"Uniform velocity U = {U}, V = {V}")

    vx, vy, vx_p, vy_p, vx_e, vy_e = init_uniform(
        U, V, config.nx, config.ny
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