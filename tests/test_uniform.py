import logging
import time
import numpy as np
import matplotlib.pyplot as plt
import sys

import yaml

from sl_python.plot import plot_blossey
from tests.test_blossey import backup, blossey_tracer

sys.path.append("/home/maurinl/sl_gt4py/src")
print(sys.path)

from sl_python.interpolation import interpolate_lin_2d
from sl_python.sl_2D import sl_init, sl_xy
from config import Config
from utils.cfl import cfl_1d
logging.basicConfig(stream=sys.stdout, level=logging.INFO)

def init_uniform(U: float, V: float, nx: int, ny: int):
    # Vitesses
    vx, vy = U * np.ones((nx, ny)), V * np.ones((nx, ny))
    vx_p, vy_p = vx.copy(), vy.copy()
    vx_e, vy_e = np.empty((nx, ny)), np.empty((nx, ny))

    return vx, vy, vx_p, vy_p, vx_e, vy_e


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
    nitmp: int,
    model_endtime: float,
    model_starttime: float
):
    tracer_ref = tracer.copy()

    # Advection
    vx_e, vy_e, vx, vy = sl_init(
        vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, lsettls=lsettls
    )

    t = model_starttime
    jstep = 0
    while t < model_endtime:
        jstep += 1
        t += config.dt
        logging.info(f"Step : {jstep}")
        logging.info(f"Time : {100*t/model_endtime:.02f}%")

        # Estimations
        tracer_e = sl_xy(
            config=config,
            vx=vx,
            vy=vy,
            vx_e=vx_e,
            vy_e=vy_e,
            tracer=tracer,
            tracer_e=tracer_e,
            interpolation_function=interpolate_lin_2d,
            nitmp=4,
        )
        
        tracer = backup(tracer=tracer, tracer_e=tracer_e)

        # Diagnostics and outputs
        courant_xmax = np.max(cfl_1d(vx_e, config.dx, config.dt))
        courant_ymax = np.max(cfl_1d(vy_e, config.dy, config.dt))

        logging.info(f"Maximum courant number : {max(courant_xmax, courant_ymax):.02f}")
        logging.info(f"Tracer : min {tracer.min():.02f}, max {tracer.max():.02f}, mean {tracer.mean():.02f}")
        
        if t >= (T / 4) and t < (T / 4) + config.dth:
            plot_blossey(config.xcr, config.ycr, vx, vy, tracer, t,  f"./figures/uniform_advection/uniform_advection_{t:.02f}.pdf")
            
        if t >= (T / 2) and t < (T / 2) + config.dth:
            plot_blossey(config.xcr, config.ycr, vx, vy, tracer, t,  f"./figures/uniform_advection/uniform_advection_{t:.02f}.pdf")
            
        if t >= (3 * T / 4) and t < (3 * T / 4) + config.dth:
            plot_blossey(config.xcr, config.ycr, vx, vy, tracer, t,  f"./figures/uniform_advection/uniform_advection_{t:.02f}.pdf")

    e_inf = np.max(np.abs(tracer - tracer_ref))
    e_2 = np.sqrt((1 / (config.nx * config.ny)) * np.sum((tracer - tracer_ref) ** 2))

    logging.info(f"Error E_inf : {e_inf}")
    logging.info(f"Error E_2 : {e_2}")

    plot_blossey(config.xcr, config.ycr, vx, vy, tracer, t, f"./figures/uniform_advection/uniform_advection_{t:.02f}.pdf")


if __name__ == "__main__":
    
    config_file = "./config/uniform_advection.yml"
    with open(config_file, 'r') as file:
        conf_dict = yaml.safe_load(file)
        config = Config(**conf_dict)
        
        
    # Shift in config file
    T = 1
    t = 0
    U, V = 1, 1
    lsettls = True
    nitmp=4
    
    logging.info("Config")
    logging.info(f"time step dt : {config.dt} s")
    logging.info(f"dx : {config.dx}, dy : {config.dy}")
    logging.info("Uniform velocity U = {U}, V = {V}")

    vx, vy, vx_p, vy_p, vx_e, vy_e = init_uniform(
        U, V, config.nx, config.ny
    )
    tracer, tracer_e = blossey_tracer(config.xcr, config.ycr)
    logging.info(f"Tracer : min {tracer.min():.02f}, max {tracer.max():.02f}, mean {tracer.mean():.02f}")
    plot_blossey(
        config.xcr, config.ycr, vx, vy, tracer, 0, f"./figures/uniform_advection/uniform_advection_{t:.02f}.pdf"
    )

    # Advection encapsulation
    start_time = time.time()
    sl_driver(config, vx, vy, vx_e, vy_e, vx_p, vy_p, tracer, tracer_e, lsettls, nitmp, config.model_endtime, config.model_starttime)
    duration = time.time() - start_time
    logging.info(f"Duration : {duration} s")
