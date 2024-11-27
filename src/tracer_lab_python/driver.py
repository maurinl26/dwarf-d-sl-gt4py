import logging
import sys
import numpy as np

from utils.config import Config
from tracer_lab_python.slag_2D_xy import smilag_transport_scheme
from tracer_lab_python.smilag_init import slag_init
from utils.blossey import blossey_tracer, blossey_velocity
from utils.budget import budget
from utils.plot import plot_blossey

logging.basicConfig(level=logging.INFO, stream=sys.stdout)
logging.getLogger()


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
    
    tracer_dep = blossey_tracer(config.xcr, config.ycr)
    vx0, vy0 = blossey_velocity(config.xcr, config.ycr, 0, config.dx, config.dy)
    plot_blossey(config.xcr, config.ycr, vx0, vy0, tracer_dep, f"./figures/tracer_0.pdf")

    
    # TODO set initial shape
    tracer_ref = tracer_dep.copy()
    # TODO save initial shape

    # Advection
    vx1, vy1 = slag_init(
        vx=vx0, vy=vy0, lsettls=lsettls
    )
    
    logging.info(f"Config : dt {config.dt}, ntimestep {config.ntimestep} ")
    
    for nstep in range(config.ntimestep):
        
        t += config.dt
        vx1, vy1 = set_velocity(config, t)

        logging.info(f"Step : {nstep}")
        logging.info(f"Time : {100*t/(config.model_endtime - config.model_starttime):.02f}%, {t:.03f} s")
    
        # Estimations
        tracer_arr = smilag_transport_scheme(
            config=config,
            vx=vx0,
            vy=vy0,
            vx_e=vx1,
            vy_e=vy1,
            tracer=tracer_dep,
            nitmp=4,
        )
 
        tracer_dep = swapp(tracer_arr=tracer_arr)
        
        if nstep % config.nfreqoutput == 0:
            plot_blossey(config.xcr, config.ycr, vx0, vy0, tracer_arr, f"./figures/tracer_{t:.3f}.pdf")

    process_output(config, tracer_arr, tracer_ref, vx0, vy0, t)
   

def set_velocity(
    config: Config,
    t: float,
):
    """Set velocity for tracers advection

    Args:
        vx (np.ndarray): _description_
        vx1 (np.ndarray): _description_
        vy (np.ndarray): _description_
        vy1 (np.ndarray): _description_
        config (Config): _description_
        t (float): _description_

    Returns:
        _type_: _description_
    """
    
    xcr, ycr = config.xcr, config.ycr
    dx, dy = config.dx, config.dy
    dt = config.dt
    
    if config.isetup == "blossey":
    
        vx1, vy1 = blossey_velocity(
            xcr, ycr, t, dx, dy
            )
    
    return vx1, vy1

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

def swapp(tracer_arr: np.ndarray):
    """Swapp tracers in memory

    Args:
        tracer (float): field at t
        tracer_e (float): field at t + dt

    Returns:
        float: _description_
    """
    tracer_dep = tracer_arr.copy()
    return tracer_dep


if __name__ == "__main__":
    
    config = Config(
        dt=0.004,
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
    
    