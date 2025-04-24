import logging
import time
import numpy as np
import sys
import yaml
from pathlib import Path
import os
from typing import Union

from sl_python.blossey import (
    blossey_tracer,
    blossey_velocity,
    init_blossey,
    tracer_shape,
)
from sl_python.cfl import cfl_1d
from sl_python.plot import plot_blossey, plot_tracer_against_reference
from sl_python.interpolation.interpolation import interpolate_cub_2d
from sl_python.sl_2D import sl_xy, sl_init
from sl_dace.sl_xy import SmiLagXY
from config import Config

logging.basicConfig(stream=sys.stdout, level=logging.INFO)


def backup(tracer, tracer_e):
    tracer = tracer_e.copy()
    return tracer


class SLDriver:
    def __init__(
        self,
        figure_dir: Path,
        config: Config,
        exec_mode: Union["python", "dace"] = "python",
    ):
        self.figure_dir = figure_dir
        self.config = config

        match exec_mode:
            case "python":
                self.sl_xy = sl_xy
            case "dace":
                self.sl_xy = SmiLagXY(grid=config.domain)

    def __call__(
        self,
        tracer: np.ndarray,
        tracer_e: np.ndarray,
        vx_e: np.ndarray,
        vy_e: np.ndarray,
        vx: np.ndarray,
        vy: np.ndarray,
        vx_p: np.ndarray,
        vy_p: np.ndarray,
    ):
        tracer_ref = tracer.copy()

        # Advection
        vx_e, vy_e, vx, vy = sl_init(
            vx_e=vx_e,
            vy_e=vy_e,
            vx=vx,
            vy=vy,
            vx_p=vx_p,
            vy_p=vy_p,
            lsettls=self.config.lsettls,
        )

        cfl_max = 0
        t = self.config.model_starttime
        T = self.config.model_endtime
        jstep = 0

        # todo: all components in timestep should be shifted to dace or compiled stencil
        while t < self.config.model_endtime:
            jstep += 1
            t += self.config.dt
            logging.info(f"Step : {jstep}")
            logging.info(
                f"Time : {100 * t / (self.config.model_endtime - self.config.model_starttime):.02f}%"
            )

            ############### Velocity ################
            vx, vy = blossey_velocity(
                self.config.xcr, self.config.ycr, t, self.config.dx, self.config.dy
            )
            vx_e, vy_e = blossey_velocity(
                self.config.xcr,
                self.config.ycr,
                t + self.config.dt,
                self.config.dx,
                self.config.dy,
            )

            ##################### SL scheme run ######################
            # todo: same signatures for sl_xy and smilag_xy
            tracer_e = self.sl_xy(
                config=self.config,
                vx=vx,
                vy=vy,
                vx_e=vx_e,
                vy_e=vy_e,
                tracer=tracer,
                tracer_e=tracer_e,
                nitmp=4,
            )

            tracer = backup(tracer=tracer, tracer_e=tracer_e)

            # Diagnostics and outputs
            cfl_max = max(
                np.max(cfl_1d(vx_e, self.config.dx, self.config.dt)),
                np.max(cfl_1d(vy_e, self.config.dy, self.config.dt)),
            )
            logging.info(f"Maximum courant number : {cfl_max:.02f}")

            # Plots
            if (
                (0.25 * T) <= t < (0.25 * T) + self.config.dt
                or (0.5 * T) <= t < (0.5 * T) + self.config.dt
                or (0.75 * T) <= t < (0.75 * T) + self.config.dt
                or t >= self.config.model_endtime - self.config.dt
            ):
                plot_blossey(
                    self.config.xcr,
                    self.config.ycr,
                    vx,
                    vy,
                    tracer,
                    Path(self.figure_dir, f"blossey_{t:.03f}.pdf"),
                )

        e_inf = np.max(np.abs(tracer - tracer_ref))
        e_2 = np.sqrt(
            (1 / (self.config.nx * self.config.ny)) * np.sum((tracer - tracer_ref) ** 2)
        )

        logging.info(f"Error E_inf : {e_inf}")
        logging.info(f"Error E_2 : {e_2}")

        plot_tracer_against_reference(
            self.config.xcr,
            self.config.ycr,
            tracer,
            tracer_ref,
            e_2,
            e_inf,
            Path(self.figure_dir, "blossey_ref.pdf"),
            cfl_max,
            self.config.dx,
        )


def init_state_blossey(config: Config, figures_dir: Path) -> dict:
    vx, vy, vx_p, vy_p, vx_e, vy_e = init_blossey(
        config.xcr,
        config.ycr,
        config.model_starttime,
        config.dt,
        config.dx,
        config.dy,
        config.nx,
        config.ny,
    )
    tracer, tracer_e = blossey_tracer(config.xcr, config.ycr)
    plot_blossey(
        config.xcr, config.ycr, vx, vy, tracer, Path(figures_dir, "blossey_0.pdf")
    )

    return {
        "vx": vx,
        "vy": vy,
        "vx_e": vx_e,
        "vy_e": vy_e,
        "vx_p": vx_p,
        "vy_p": vy_p,
        "tracer": tracer,
        "tracer_e": tracer_e,
    }


def set_experiment(config_file_name: str) -> Config:
    current_path = Path(__file__)
    root_dir = current_path.parent.parent.parent

    # set paths
    figures_path = Path(root_dir, "figures")
    run_figures_path = Path(figures_path, "durran_blossey")
    if not figures_path.exists():
        os.mkdir(figures_path)
    if not run_figures_path.exists():
        os.mkdir(run_figures_path)

    # load config file
    config_file = Path(root_dir, "config", config_file_name)
    with open(config_file, "r") as file:
        conf_dict = yaml.safe_load(file)
        return Config(**conf_dict), run_figures_path


if __name__ == "__main__":
    # Shift in config file
    config, run_figures_path = set_experiment("durran_blossey.yml")
    state = init_state_blossey(config, run_figures_path)
    # Advection encapsulation
    start_time = time.time()
    SLDriver(figure_dir=run_figures_path, config=config, exec_mode="dace")(**state)
    duration = time.time() - start_time
    logging.info(f"Duration : {duration} s")
