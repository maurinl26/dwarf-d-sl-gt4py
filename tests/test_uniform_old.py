import sys
import logging
from matplotlib import pyplot as plt
import numpy as np

sys.path.append("/home/maurinl/FVM_GT4Py_slim/src")
sys.path.append("/home/maurinl/sl_gt4py/src")

from sl_python.sl_2D import sl_init, sl_xy, backup
from sl_python.plot import plot_2D_scalar_field
from sl_python.interpolation import interpolate_cubic_2d, cubic_interpolation
from config import Config

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

def run_sl_2D(lsettls: bool, nsiter: int, nitmp: int):
    lnesc = not lsettls

    dt = 1
    
    # Grid
    xmin, xmax = 0, 29
    ymin, ymax = 0, 19
    nx, ny = 30, 20

    # Boundaries
    # 1 : PERIODIC
    # 0 : FIXED
    bcx_kind = 1
    bcy_kind = 1
    
    config = Config(
        dt, xmin, xmax, nx, ymin, ymax, ny, bcx_kind, bcy_kind
    )
        
    # Initialisation simple
    # Vent uniforme
    U, V = 1.5, 0

    ############## Declaration des champs #############
    vx, vy = (U * np.ones((nx, ny)), V * np.ones((nx, ny)))
    
    logging.debug(f"Shape vx : {vx.shape}, Shape vy : {vy.shape}")

    # Champs de vent à t - dt
    vx_p, vy_p = vx.copy(), vy.copy()

    # Champs de vent à t + dt
    vx_e, vy_e = np.zeros((nx, ny)), np.zeros((nx, ny))
    
    # Tracer Gaussien
    T = 20
    rsq = (config.xcr - 6) ** 2 + (config.ycr - 7) ** 2

    tracer = T * np.exp(-(rsq / (2 * 2)))
    tracer_e = tracer.copy()

    logging.info("")    
    ###### Plot Initial Fields #####
    plot_2D_scalar_field(config.xcr, config.ycr, tracer, 25)
    logging.debug(f"Step : Initialisation, Time : {0} s")
    logging.debug(
        f"Tracer : min {np.min(tracer)}, max {np.max(tracer)}, mean {np.mean(tracer)}"
    )
    imax, jmax = np.unravel_index(np.argmax(tracer, axis=None), tracer.shape)
    logging.debug(f"Max of tracer field,  X : {config.xcr[imax, jmax]:.2f} Y : {config.ycr[imax, jmax]:.2f}")

    ########## Advection ###########
    ######### Premier pas ##########

    ######### jstep = 0 ##########

    # Initialisation vitesses
    vx_e, vy_e, vx, vy = sl_init(
        vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, lsettls=lsettls
    )

    ######### jstep > 0 ##########
    for jstep in range(nitmp):
        # Copie des champs
        vx, vy, tracer = backup(
            vx=vx, vy=vy, vx_e=vx_e, vy_e=vy_e, tracer=tracer, tracer_e=tracer_e
        )

        # Estimations
        vx_e, vy_e, tracer_e = sl_xy(
            config=config,
            vx=vx,
            vy=vy,
            vx_e=vx_e,
            vy_e=vy_e,
            tracer=tracer,
            tracer_e=tracer_e,
            interpolation_function=interpolate_cubic_2d,
            nsiter=nsiter
        )
        
        plot_2D_scalar_field(config.xcr, config.ycr, tracer, 25)

        logging.debug(f"Step : {jstep}, Time : {jstep * dt} s")
        logging.debug(
            f"Tracer : min {np.min(tracer)}, max {np.max(tracer)}, mean {np.mean(tracer)}"
        )

        imax, jmax = np.unravel_index(np.argmax(tracer, axis=None), tracer.shape)
        logging.debug(f"Max of tracer field,  X : {config.xcr[imax, jmax]:.2f} Y : {config.ycr[imax, jmax]:.2f}")


    #### Print ####
    plot_2D_scalar_field(config.xcr, config.ycr, tracer, 25)

    logging.debug(
        f"Tracer : min {np.min(tracer)}, max {np.max(tracer)}, mean {np.mean(tracer)}"
    )
    imax, jmax = np.unravel_index(np.argmax(tracer, axis=None), tracer.shape)
    logging.debug(f"Max of tracer field,  X : {config.xcr[imax, jmax]:.2f} Y : {config.ycr[imax, jmax]:.2f}")


if __name__ == "__main__":
    
    LSETTLS = True
    NITMP = 21
    NSITER = 5

    run_sl_2D(LSETTLS, NSITER, NITMP)


