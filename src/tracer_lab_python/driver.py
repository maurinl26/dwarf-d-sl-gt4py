import logging

import numpy as np

from config import Config
from tracer_lab_python.slag_2D_xy import smilag_transport_scheme
from tracer_lab_python.smilag_init import slag_init
from utils.blossey import blossey_tracer, blossey_velocity
from utils.budget import budget
from utils.plot import plot_blossey, plot_tracer_against_reference

logging.getLogger(__name__)

def slag_driver(
    config: Config,
    lsettls: bool = False,
):
    """Driver for blossey advection test

    Args:
        config (Config): configuration file
        lsettls (bool): velocity initialization
    """
    
    t = config.model_starttime
    
    tracer = blossey_tracer(config.xcr, config.ycr)
    vx, vy = blossey_velocity(config.xcr, config.ycr, 0, config.dx, config.dy)
    
    # TODO set initial shape
    tracer_ref = tracer.copy()
    # TODO save initial shape

    # Advection
    vx_e, vy_e = slag_init(
        vx=vx, vy=vy, lsettls=lsettls
    )
    
    for nstep in range(config.ntimestep):
        
        t += config.dt
        logging.info(f"Step : {nstep}")
        logging.info(f"Time : {100*t/(config.model_endtime - config.model_starttime):.02f}%, {t} s")
        
        vx, vy, vx_e, vy_e = set_velocity(vx, vx_e, vy, vy_e, config, t)

        # Estimations
        tracer_e = smilag_transport_scheme(
            config=config,
            vx=vx,
            vy=vy,
            vx_e=vx_e,
            vy_e=vy_e,
            tracer=tracer,
            nitmp=4,
        )
        
        # Backup is swapp in fortran
        tracer = swapp(tracer=tracer, tracer_e=tracer_e)
        
        if nstep % config.nfreqoutput == 0:
            plot_blossey(config.xcr, config.ycr, vx, vy, tracer, f"./figures/tracer_{t:.03f}.pdf")

    process_output(config, tracer, tracer_ref, vx, vy, t)
   

def set_velocity(
    vx: np.ndarray,
    vx_e: np.ndarray,
    vy: np.ndarray,
    vy_e: np.ndarray,
    config: Config,
    t: float,
):
    """Set velocity for tracers advection

    Args:
        vx (np.ndarray): _description_
        vx_e (np.ndarray): _description_
        vy (np.ndarray): _description_
        vy_e (np.ndarray): _description_
        config (Config): _description_
        t (float): _description_

    Returns:
        _type_: _description_
    """
    
    xcr, ycr = config.xcr, config.ycr
    dx, dy = config.dx, config.dy
    dt = config.dt
    
    vx, vy = blossey_velocity(xcr, ycr, t, dx, dy)
    vx_e, vy_e = blossey_velocity(
            xcr, ycr, t + dt, dx, dy
        )
    
    return vx, vy, vx_e, vy_e

def process_output(
    config: Config,
    tracer: np.ndarray,
    tracer_ref: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    t: np.ndarray
):
    """Compute error and plot advected field against reference.

    Args:
        config (Config): config
        tracer (np.ndarray): advected tracer field
        tracer_ref (np.ndarray): reference tracer field
        vx (np.ndarray): contravariant velocity on x
        vy (np.ndarray): contravariant velocity on y
        t (np.ndarray): time of simulation (s)
    """
    
    e_inf = np.max(np.abs(tracer - tracer_ref))
    e_2 = np.sqrt((1 / (config.nx * config.ny)) * np.sum((tracer - tracer_ref) ** 2))
    tracer_ratio = budget(tracer, tracer_ref)

    logging.info(f"Ratio : {tracer_ratio * 100:.03f} %")
    logging.info(f"Error E_inf : {e_inf}")
    logging.info(f"Error E_2 : {e_2}")
    
    plot_blossey(config.xcr, config.ycr, vx, vy, tracer, f"./figures/blossey_{t:.03f}.pdf")
    # plot_tracer_against_reference(config.xcr, config.ycr, tracer, tracer_ref, e_2, e_inf, f"./figures/blossey_ref.pdf", cfl_max, config.dx)

def swapp(tracer, tracer_e):
    """Swapp tracers in memory

    Args:
        tracer (float): field at t
        tracer_e (float): field at t + dt

    Returns:
        float: _description_
    """
    tracer = tracer_e.copy()
    return tracer


if __name__ == "__main__":
    
    config = Config(
        dt=1,
        xmin=0,
        xmax=1,
        ymin=0,
        ymax=1,
        nx=200,
        ny=200,
        bcx_kind=1,
        bcy_kind=1,
        ntimestep=250,
        nfreqoutput=25
    )
    
    slag_driver(config, lsettls=False)
    
    