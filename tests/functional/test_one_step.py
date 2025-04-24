import logging
import time
import numpy as np
import sys

from sl_python.plot import plot_blossey
from sl_python.interpolation.interpolation import interpolate_cub_2d, interpolate_lin_2d
from sl_python.sl_2D import sl_init, sl_xy
from test_blossey import blossey_tracer, set_experiment
from config import Config
from sl_python.cfl import cfl_1d

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logging.getLogger(__name__)

def init_uniform(U: float, V: float, nx: int, ny: int):
    # velocity
    vx, vy = U * np.ones((nx, ny)), V * np.ones((nx, ny))
    vx_p, vy_p = vx.copy(), vy.copy()
    vx_e, vy_e = np.empty((nx, ny)), np.empty((nx, ny))

    return vx, vy, vx_p, vy_p, vx_e, vy_e

class OneStepTest:

    def __init__(self, config: Config,
                 interp_function: callable = interpolate_lin_2d):
        self.config = config
        self.interp_function = interp_function

    def __call__(self,
                 vx: np.ndarray,
                 vy: np.ndarray,
                 vx_e: np.ndarray,
                 vy_e: np.ndarray,
                 vx_p: np.ndarray,
                 vy_p: np.ndarray,
                 tracer: np.ndarray,
                 tracer_e: np.ndarray,

                 ):

        tracer_ref = tracer.copy()

        # Advection
        vx_e, vy_e, vx, vy = sl_init(
            vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, lsettls=self.config.lsettls
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
            interpolation_function=self.interp_function,
            nitmp=4,
        )

        # Diagnostics and outputs
        cfl_max = max(
            np.max(cfl_1d(vx_e, config.dx, config.dt)),
            np.max(cfl_1d(vy_e, config.dy, config.dt))
        )

        logging.info(f"Maximum courant number : {cfl_max:.02f}")
        logging.info(
            f"Tracer (t + dt) : min {tracer_e.min():.02f}, max {tracer_e.max():.02f}, mean {tracer_e.mean():.02f}")

        plot_blossey(config.xcr, config.ycr, vx, vy, tracer_e, f"./figures/one_step/one_step_end.pdf")

        # Errors
        e_inf = np.max(np.abs(tracer - tracer_ref))
        e_2 = np.sqrt((1 / (config.nx * config.ny)) * np.sum((tracer - tracer_ref) ** 2))

        logging.info(f"Error E_inf : {e_inf}")
        logging.info(f"Error E_2 : {e_2}")

        return tracer_e


def init_one_step_state(config: Config):

        U, V = 0.2, 0.2

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

        return {
            "vx": vx,
            "vy": vy,
            "vx_p": vx_p,
            "vy_p": vy_p,
            "vx_e": vx_e,
            "vy_e": vy_e,
            "tracer": tracer,
            "tracer_e": tracer_e,
            "tracer_ref": tracer_ref
        }

if __name__ == "__main__":

    config, run_figures_path = set_experiment("uniform_advection.yml")
    #### Interp Lineaire
    state = init_one_step_state(config)

    # Advection encapsulation
    start_time = time.time()
    OneStepTest(config, interpolate_lin_2d)(**state)
    duration = time.time() - start_time
    logging.info(f"Duration : {duration:.02f} s")

    #### Interp Cubique
    state = init_one_step_state(config)

    start_time = time.time()
    tracer  = OneStepTest(config, interpolate_cub_2d)(**state)
    duration = time.time() - start_time
    logging.info(f"Duration : {duration:.02f} s")
