import logging
import numpy as np
import sys

from sl_python.blossey import blossey_tracer, blossey_velocity, init_blossey, tracer_shape
from utils.cfl import cfl_1d
from sl_python.plot import plot_blossey, plot_tracer_against_reference
from sl_python.interpolation import interpolate_cub_2d
from sl_python.sl_2D import sl_xy, sl_init
from config import Config
from utils.cfl import cfl_1d

logging.basicConfig(stream=sys.stdout, level=logging.INFO)


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
    T: float = 1
):
    tracer_ref = tracer.copy()

    # Advection
    vx_e, vy_e, vx, vy = sl_init(
        vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, lsettls=lsettls
    )
    
    cfl_max = 0

    t = model_starttime
    jstep = 0
    while t < model_endtime:
        
        jstep += 1
        t += config.dt
        logging.info(f"Step : {jstep}")
        logging.info(f"Time : {100*t/(model_endtime - model_starttime):.02f}%")

        vx, vy = blossey_velocity(config.xcr, config.ycr, t, config.dx, config.dy)
        vx_e, vy_e = blossey_velocity(
            config.xcr, config.ycr, t + config.dt, config.dx, config.dy
        )

        # Estimations
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

        tracer = backup(tracer=tracer, tracer_e=tracer_e)

        # Diagnostics and outputs
        courant_xmax = np.max(cfl_1d(vx_e, config.dx, config.dt))
        courant_ymax = np.max(cfl_1d(vy_e, config.dy, config.dt))
        
        if courant_xmax > cfl_max:
            cfl_max = courant_xmax
        if courant_ymax > cfl_max:
            cfl_max = courant_ymax

        logging.info(f"Maximum courant number : {max(courant_xmax, courant_ymax):.02f}")
        
        if t >= (T / 4) and t < (T / 4) + config.dt:
            plot_blossey(config.xcr, config.ycr, vx, vy, tracer, f"./figures/blossey/blossey_{t:.03f}.pdf")

        if t >= (T / 2) and t < (T / 2) + config.dt:
            plot_blossey(config.xcr, config.ycr, vx, vy, tracer, f"./figures/blossey/blossey_{t:.03f}.pdf")
            
        if t >= (0.75 * T) and t < (0.75 * T) + config.dt:
            plot_blossey(config.xcr, config.ycr, vx, vy, tracer, f"./figures/blossey/blossey_{t:.03f}.pdf")

    e_inf = np.max(np.abs(tracer - tracer_ref))
    e_2 = np.sqrt((1 / (config.nx * config.ny)) * np.sum((tracer - tracer_ref) ** 2))

    logging.info(f"Error E_inf : {e_inf}")
    logging.info(f"Error E_2 : {e_2}")

    plot_blossey(config.xcr, config.ycr, vx, vy, tracer, f"./figures/blossey/blossey_{t:.03f}.pdf")
    
    plot_tracer_against_reference(config.xcr, config.ycr, tracer, tracer_ref, e_2, e_inf, f"./figures/blossey/blossey_ref.pdf", cfl_max, config.dx)

    
