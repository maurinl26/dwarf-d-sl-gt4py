import logging
import time
import numpy as np
import sys

from sl_python.plot import plot_blossey
from test_blossey import backup, blossey_tracer
from sl_python.interpolation.interpolation import interpolate_lin_2d
from sl_python.sl_2D import sl_init, sl_xy
from config import Config
from sl_python.cfl import cfl_1d
from test_blossey import set_experiment

logging.basicConfig(stream=sys.stdout, level=logging.INFO)

def init_uniform(U: float, V: float, nx: int, ny: int):
    # Vitesses
    vx, vy = U * np.ones((nx, ny)), V * np.ones((nx, ny))
    vx_p, vy_p = vx.copy(), vy.copy()
    vx_e, vy_e = np.empty((nx, ny)), np.empty((nx, ny))

    return vx, vy, vx_p, vy_p, vx_e, vy_e

def init_uniform_state(config: Config):
    U, V = 1, 1

    logging.info("Config")
    logging.info(f"time step dt : {config.dt} s")
    logging.info(f"dx : {config.dx}, dy : {config.dy}")
    logging.info(f"Uniform velocity U = {U}, V = {V}")

    vx, vy, vx_p, vy_p, vx_e, vy_e = init_uniform(
        U, V, config.nx, config.ny
    )
    tracer, tracer_e = blossey_tracer(config.xcr, config.ycr)
    logging.info(f"Tracer : min {tracer.min():.02f}, max {tracer.max():.02f}, mean {tracer.mean():.02f}")
    plot_blossey(
        config.xcr, config.ycr, vx, vy, tracer,
        f"./figures/uniform_advection/uniform_advection_{config.model_starttime:.02f}.pdf"
    )

    return {
        "vx": vx,
        "vy": vy,
        "vx_p": vx_p,
        "vy_p": vy_p,
        "vx_e": vx_e,
        "vy_e": vy_e,
        "tracer": tracer,
        "tracer_e": tracer_e
    }

class UniformTestDriver:

    def __init__(self, config: Config):
        self.config = config

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

        t = self.config.model_starttime
        t_end = self.config.model_endtime
        jstep = 0
        while t < self.config.model_endtime:
            jstep += 1
            t += self.config.dt
            logging.info(f"Step : {jstep}")
            logging.info(f"Time : {100 * t / self.config.model_endtime:.02f}%")

            # Estimations
            tracer_e = sl_xy(
                config=self.config,
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
            cfl_max = max(
                np.max(cfl_1d(vx_e, self.config.dx, self.config.dt)),
                np.max(cfl_1d(vy_e, self.config.dy, self.config.dt))
            )

            logging.info(f"Maximum courant number : {cfl_max:.02f}")
            logging.info(f"Tracer : min {tracer.min():.02f}, max {tracer.max():.02f}, mean {tracer.mean():.02f}")

            if (
                    0.25 * t_end <= t < 0.25 * t_end + self.config.dt
                    or 0.5 * t_end <= t < 0.5 * t_end + self.config.dt
                    or 0.75 * t_end <= t < 0.75 * t_end + self.config.dt
            ):
                plot_blossey(self.config.xcr, self.config.ycr, vx, vy, tracer,
                             f"./figures/uniform_advection/uniform_advection_{t:.02f}.pdf")

        e_inf = np.max(np.abs(tracer - tracer_ref))
        e_2 = np.sqrt((1 / (self.config.nx * self.config.ny)) * np.sum((tracer - tracer_ref) ** 2))

        logging.info(f"Error E_inf : {e_inf}")
        logging.info(f"Error E_2 : {e_2}")

        plot_blossey(self.config.xcr, self.config.ycr, vx, vy, tracer,
                     f"./figures/uniform_advection/uniform_advection_{t:.02f}.pdf")


if __name__ == "__main__":

    config, run_figures_path = set_experiment("uniform_advection.yml")

    state = init_uniform_state(config)

    # Advection encapsulation
    start_time = time.time()
    UniformTestDriver(config)(**state)
    duration = time.time() - start_time
    logging.info(f"Duration : {duration} s")
