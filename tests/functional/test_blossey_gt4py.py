import logging
import time
import numpy as np
import sys
import yaml
import gt4py

from config import Config
from sl_gt4py.build import dtype, backend
from sl_gt4py.blossey_velocity import init_blossey
from sl_python.plot import plot_blossey, plot_tracer_against_reference

logging.basicConfig(stream=sys.stdout, level=logging.INFO)


if __name__ == "__main__":
    # Shift in config file
    
    config_file = "./config/durran_blossey.yml"
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

    fields = [
        "vx",
        "vy",
        "vx_p",
        "vy_p",
        "vx_e",
        "vy_e",
        "tracer",
        "tracer_e"
    ]
    
    aligned_index = (0, 0)
    
    # Velocities
    vx = gt4py.storage.ones((config.nx, config.ny), dtype, backend=backend, aligned_index=aligned_index)
    vy = gt4py.storage.ones((config.nx, config.ny), dtype, backend=backend, aligned_index=aligned_index)
    vx_p = gt4py.storage.ones((config.nx, config.ny), dtype, backend=backend, aligned_index=aligned_index)
    vy_p = gt4py.storage.ones((config.nx, config.ny), dtype, backend=backend, aligned_index=aligned_index)
    vx_e = gt4py.storage.ones((config.nx, config.ny), dtype, backend=backend, aligned_index=aligned_index)
    vy_e = gt4py.storage.ones((config.nx, config.ny), dtype, backend=backend, aligned_index=aligned_index)
    # Tracer
    tracer = gt4py.storage.ones((config.nx, config.ny), dtype, backend=backend, aligned_index=aligned_index)
    tracer_e = gt4py.storage.ones((config.nx, config.ny), dtype, backend=backend, aligned_index=aligned_index)
    
    # Init velocities (numpy -> gt4py)
    init_state = init_blossey(
         config.xcr, config.ycr, t, config.dt, config.dx, config.dy, config.nx, config.ny
    )
    vx[:, :] = init_state[0]
    vy[:, :] = init_state[1]
    vx_p[:, :] = init_state[2]
    vy_p[:, :] = init_state[3]
    vx_e[:, :] = init_state[4]
    vy_e[:, :] = init_state[5]    
    
    
    tracer_state = blossey_tracer(config.xcr, config.ycr)
    tracer[:, :] = tracer_state[0]
    tracer_e[:, :] = tracer_state[1]

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
    
    # Pickup fields
    out_tracer = np.asarray(tracer)
    out_vx = np.asarray(vx)
    out_vy = np.asarray(vy)
    
    e_inf = np.max(np.abs(tracer - tracer_state[0]))
    e_2 = np.sqrt((1 / (config.nx * config.ny)) * np.sum((tracer - tracer_state[0]) ** 2))

    plot_blossey(config.xcr, config.ycr, out_vx, out_vy, out_tracer, 1)
    plot_tracer_against_reference(config.xcr, config.yxr, out_tracer, tracer_state[0], e_2, e_inf)
