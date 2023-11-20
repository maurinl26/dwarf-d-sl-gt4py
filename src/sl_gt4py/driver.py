import logging
import numpy as np

from config import Config
from sl_gt4py.sl_2D import sl_init, sl_xy
from sl_python.sl_2D import backup
from tests.test_blossey import blossey_velocity
from utils.cfl import cfl_1d

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
            nitmp=4,
        )

        # TODO : backup in GT4Py
        tracer = backup(tracer=tracer, tracer_e=tracer_e)

        # Diagnostics and outputs
        # TODO : cfl1d in GT4Py
        courant_xmax = np.max(cfl_1d(vx_e, config.dx, config.dt))
        courant_ymax = np.max(cfl_1d(vy_e, config.dy, config.dt))

        logging.info(f"Maximum courant number : {max(courant_xmax, courant_ymax):.02f}")