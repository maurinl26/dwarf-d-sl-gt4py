import logging
import time
import numpy as np
import sys
import yaml
from pathlib import Path
import os
import dace

from sl_dace.blossey import (
    blossey_tracer,
    blossey_velocity,
    init_blossey,
    tracer_shape,
)
from sl_dace.cfl import cfl_1d
from sl_dace.plot import plot_blossey, plot_tracer_against_reference
from sl_dace.interpolation.interpolation import interpolate_cub_2d
from sl_dace.sl_xy import sl_xy
from sl_dace.sl_init import sl_init
from config import Config

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
    figure_dir: Path,
):
    tracer_ref = tracer.copy()

    # Advection
    vx_e, vy_e, vx, vy = sl_init(
        vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, lsettls=config.lsettls
    )

    cfl_max = 0

    t = config.model_starttime
    T = config.model_endtime - config.model_starttime
    jstep = 0
    while t < config.model_endtime:
        jstep += 1
        t += config.dt
        logging.info(f"Step : {jstep}")
        logging.info(
            f"Time : {100 * t / (config.model_endtime - config.model_starttime):.02f}%"
        )

        ###################### Velocity precription ################
        vx, vy = blossey_velocity(config.xcr, config.ycr, t, config.dx, config.dy)
        vx_e, vy_e = blossey_velocity(
            config.xcr, config.ycr, t + config.dt, config.dx, config.dy
        )

        ##################### SL scheme run ######################
        tracer_e = sl_xy(
            bcx_kind=config.bcx_kind,
            bcy_kind=config.bcy_kind,
            nx=config.nx,
            ny=config.ny,
            dx=config.dx,
            dy=config.dy,
            dth=config.dth,
            I=config.I,
            J=config.J,
            vx=vx,
            vy=vy,
            vx_e=vx_e,
            vy_e=vy_e,
            tracer=tracer,
            tracer_e=tracer_e,
            interpolation_function=interpolate_cub_2d,
            nitmp=4,
        )

        # backup
        tracer = backup(tracer=tracer, tracer_e=tracer_e)

        #################### Diagnostics and outputs ####################
        courant_xmax = np.max(cfl_1d(vx_e, config.dx, config.dt))
        courant_ymax = np.max(cfl_1d(vy_e, config.dy, config.dt))

        if courant_xmax > cfl_max:
            cfl_max = courant_xmax
        if courant_ymax > cfl_max:
            cfl_max = courant_ymax

        logging.info(f"Maximum courant number : {max(courant_xmax, courant_ymax):.02f}")

        if (
            (T / 4) <= t < (T / 4) + config.dt
            or (T / 2) <= t < (T / 2) + config.dt
            or (0.75 * T) <= t < (0.75 * T) + config.dt
        ):
            plot_blossey(
                config.xcr,
                config.ycr,
                vx,
                vy,
                tracer,
                Path(figure_dir, f"blossey_{t:.03f}.pdf"),
            )

    e_inf = np.max(np.abs(tracer - tracer_ref))
    e_2 = np.sqrt((1 / (config.nx * config.ny)) * np.sum((tracer - tracer_ref) ** 2))

    logging.info(f"Error E_inf : {e_inf}")
    logging.info(f"Error E_2 : {e_2}")

    plot_blossey(
        config.xcr, config.ycr, vx, vy, tracer, Path(figure_dir, f"blossey_{t:.03f}.pdf")
    )

    plot_tracer_against_reference(
        config.xcr,
        config.ycr,
        tracer,
        tracer_ref,
        e_2,
        e_inf,
        Path(figure_dir, "blossey_ref.pdf"),
        cfl_max,
        config.dx,
    )

def run(config: Config):

        # Pour info
        T = config.model_endtime - config.model_starttime
        t = config.model_starttime
        nstep = np.ceil((config.model_endtime - config.model_starttime) / config.dt)

        logging.info(f"Time step dt : {config.dt:.06f} s")
        logging.info(f"N steps : {nstep:.06f} s")

        vx, vy, vx_p, vy_p, vx_e, vy_e = init_blossey(
            config.xcr, config.ycr, t, config.dt, config.dx, config.dy, config.nx, config.ny
        )
        tracer, tracer_e = blossey_tracer(config.xcr, config.ycr)
        plot_blossey(
            config.xcr, config.ycr, vx, vy, tracer, Path(run_figures_path, "blossey_0.pdf")
        )

        state = {
            "vx": vx,
            "vy": vy,
            "vx_e": vx_e,
            "vy_e": vy_e,
            "vx_p": vx_p,
            "vy_p": vy_p,
            "tracer": tracer,
            "tracer_e": tracer_e,
        }

        # Advection encapsulation
        start_time = time.time()
        sl_driver(
            config=config,
            figure_dir=run_figures_path,
            **state,
        )
        duration = time.time() - start_time
        logging.info(f"Duration : {duration} s")


if __name__ == "__main__":
    # Shift in config file

    current_path = Path(__file__)
    root_dir = current_path.parent.parent.parent
    config_file = Path(root_dir, "config", "durran_blossey.yml")
    with open(config_file, "r") as file:
        conf_dict = yaml.safe_load(file)
        config = Config(**conf_dict)

    figures_path = Path(root_dir, "figures")
    run_figures_path = Path(figures_path, "durran_blossey")
    if not figures_path.exists():
        os.mkdir(figures_path)
    if not run_figures_path.exists():
        os.mkdir(run_figures_path)

    # Dace compilation
    #logging.info("Compiling to dace")
    #sl_xy_dace = dace.program(sl_xy)
    #sl_xy_sdfg = sl_xy_dace.to_sdfg()


    logging.info("Run sl advection")
    run(config)

