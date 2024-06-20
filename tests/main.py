import logging
import time
import numpy as np
import yaml
from config import Config
from sl_python.driver import sl_driver
from sl_python.blossey import blossey_tracer, init_blossey
from sl_python.plot import plot_blossey
import typer

from test_one_step import one_step_driver
from test_uniform import init_uniform 

app = typer.Typer()

@app.command()
def run_blossey(config_file: str):
    """Run blossey experiment

    Args:
        config_file (str): _description_
    """
    
    # config_file = "./config/durran_blossey.yml"
    with open(config_file, 'r') as file:
        conf_dict = yaml.safe_load(file)
        
        config = Config(**conf_dict)
    
    # LSETTLS  
    lsettls = True
    
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
    plot_blossey(config.xcr, config.ycr, vx, vy, tracer, "./figures/blossey/blossey_0.pdf")

    # Advection encapsulation
    start_time = time.time()
    sl_driver(
        config,
        vx,
        vy,
        vx_e,
        vy_e,
        vx_p,
        vy_p,
        tracer,
        tracer_e,
        lsettls,
        config.model_starttime,
        config.model_endtime,
    )
    duration = time.time() - start_time
    logging.info(f"Duration : {duration} s")
    
@app.command()
def run_uniform(config_file: str):
    """Run 2D advection with uniform prescribed verlocity

    Args:
        config_file (str): . config file yml
    """
    
    # config_file = "./config/uniform_advection.yml"
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
    
    
@app.command()
def run_one_step(config_file: str):
    """Run one step for uniform advection.

    Args:
        config_file (str): .yml configuration file
    """
    
    # config_file = "./config/uniform_advection.yml"
    with open(config_file, 'r') as file:
        conf_dict = yaml.safe_load(file)
        config = Config(**conf_dict)

    # TODO : send config to config field (.yml)
    U, V = 0.2, 0.2
    lsettls = True

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
    
    #### Interp Lineaire

    # Advection encapsulation
    start_time = time.time()
    tracer  = one_step_driver(config, vx, vy, vx_e, vy_e, vx_p, vy_p, tracer, tracer_e, lsettls, interpolate_lin_2d)
    duration = time.time() - start_time
    logging.info(f"Duration : {duration:.02f} s")
    
    # Errors 
    e_inf = np.max(np.abs(tracer - tracer_ref))
    e_2 = np.sqrt((1 / (config.nx * config.ny)) * np.sum((tracer - tracer_ref) ** 2))

    logging.info(f"Error E_inf : {e_inf}")
    logging.info(f"Error E_2 : {e_2}")
    
    #### Interp Cubique
    
    # Advection encapsulation
    tracer, tracer_e = blossey_tracer(config.xcr, config.ycr)
    tracer_ref = tracer.copy()
    
    start_time = time.time()
    tracer  = one_step_driver(config, vx, vy, vx_e, vy_e, vx_p, vy_p, tracer, tracer_e, lsettls, interpolate_lin_2d)
    duration = time.time() - start_time
    logging.info(f"Duration : {duration:.02f} s")
    
    # Errors 
    e_inf = np.max(np.abs(tracer - tracer_ref))
    e_2 = np.sqrt((1 / (config.nx * config.ny)) * np.sum((tracer - tracer_ref) ** 2))

    logging.info(f"Error E_inf : {e_inf}")
    logging.info(f"Error E_2 : {e_2}")

if __name__ == "__main__":
    app()