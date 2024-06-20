import logging
import numpy as np
import sys

from sl_python.plot import plot_blossey
from sl_python.interpolation import interpolate_lin_2d
from sl_python.sl_2D import sl_init, sl_xy
from config import Config
from utils.cfl import cfl_1d

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logging.getLogger(__name__)

def init_uniform(U: float, V: float, nx: int, ny: int):
    """Initialize horizontal velocities with uniform values (U, V)

    Args:
        U (float): _description_
        V (float): _description_
        nx (int): _description_
        ny (int): _description_

    Returns:
        _type_: _description_
    """
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

