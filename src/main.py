""" Taken from FVM_GT4Py_slim blossey_xy """
import sys
import logging
import time
from matplotlib import pyplot as plt
import numpy as np


sys.path.append("/home/maurinl/FVM_GT4Py_slim/src")
sys.path.append("/home/maurinl/SL_GT4Py/SL8GT4Py/src")

from fvms.model.config import Config
from fvms.model.fields import FieldContainer
from fvms.utils.storage import to_numpy

from sl_python.run_model_2D import sl_init, sl_xy, backup
from sl_python.plot import plot_2D_scalar_field

def model_driver(config: Config, fields: FieldContainer):
    """Driver for semi lagrangian advection.
    
    Args:
        config (Config): _description_
        fields (FieldContainer): _description_

    Returns:
        _type_: _description_
    """
    
    dt = config.dt
    NITMP = config.outer_nstep
    
    # Initiate 
    tracer = to_numpy(fields.tracer["covering"])
    tracer_e = to_numpy(fields.tracer["covering"])
    
    vx, vy = to_numpy(fields.vel[0]["covering"]), to_numpy(fields.vel[1]["covering"])
    vx_e, vy_e = to_numpy(fields.vel[0]["covering"]), to_numpy(fields.vel[1]["covering"])
    vx_p, vy_p = to_numpy(fields.vel_bck[0]["covering"]), to_numpy(fields.vel_bck[1]["covering"])
    
    # Coordinates
    xcr = config.coordinates.xcr["covering"]
    ycr = config.coordinates.ycr["covering"]
    
    # Indices 
    I = np.arange(0, nx).astype(np.int8)
    J = np.arange(0, ny).astype(np.int8)
    
    # 
    vx_e, vy_e, vx, vy = sl_init(
        vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, LSETTLS=True
    )

    ######### j_iter > 0 ##########
    for jstep in range(1, NITMP):
        
        # Copie des champs
        tracer, vy, tracer = backup(
            vx=vx,
            vy=vy,
            vx_e=vx_e,
            vy_e=vy_e,
            tracer=tracer,
            tracer_e=tracer_e
        )


        # Estimations
        vx_e, vy_e, tracer_e = sl_xy(
            I=I,
            J=J,
            x=xcr,
            y=ycr,
            vx=vx,
            vy=vy,
            vx_e=vx_e,
            vy_e=vy_e,
            tracer=tracer,
            tracer_e=tracer_e,
            dt=dt,
            dx=dx, 
            dy=dy,
            nx=nx,
            ny=ny,
            bcx_kind=bcx_kind,
            bcy_kind=bcy_kind,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax
        )
    
    return vx_e, vy_e, tracer_e

if __name__ == "__main__":
    
    # Init option
    LSETTLS = True
    LNESC = not LSETTLS  # Pour info

    # Iterations
    NSITER = 5  # Semi lagrangian step
    NITMP = 5
    dt = 1

    # Grid
    xmin, xmax = 0, 100
    ymin, ymax = 0, 100
    nx, ny = 100, 100

    # spacings
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny

    # Boundaries
    # 1 : PERIODIC
    # 0 : FIXED
    bcx_kind = 0
    bcy_kind = 0

    # Spacing
    xc = np.linspace(xmin, xmax, nx)
    yc = np.linspace(ymin, ymax, ny)
    xcr, ycr = np.meshgrid(xc, yc)

    # Horizontal indexes
    i_indices = np.arange(0, nx).astype(np.int8)
    j_indices = np.arange(0, ny).astype(np.int8)
    I, J = np.meshgrid(i_indices, j_indices)

    # Initialisation simple
    # Vent uniforme
    U, V = 20, 20

    ############## Declaration des champs #############
    vx, vy = (U * np.ones((nx, ny)), V * np.ones((nx, ny)))

    # Champs de vent à t - dt
    vx_p, vy_p = vx.copy(), vy.copy()

    # Champs de vent à t + dt
    vx_e, vy_e = (np.zeros((nx, ny)), np.zeros((nx, ny)))

    # Tracer Gaussien
    T = 20
    rsq = (xcr - 25) ** 2 + (ycr - 25) ** 2

    tracer = T * np.exp(-(rsq / (2 * 10)))
    tracer_e = tracer.copy()
    
    
    ###### Plot Initial Fields #####
    fig = plot_2D_scalar_field(xcr, ycr, tracer, 25)
    plt.show()
    print(f"Step : {0}, Time : {0} s")
    print(f"Tracer : min {np.min(tracer)}, max {np.max(tracer)}, mean {np.mean(tracer)}")
    imax, jmax = np.unravel_index(np.argmax(tracer, axis=None), tracer.shape)
    print(f"Max of tracer field,  X : {xc[imax]:.2f} Y : {yc[jmax]:.2f}")


    ########### Advection ##########
    ######### Premier pas ##########
    
    ######### jstep = 0 ##########

    # Initialisation vitesses
    vx_e, vy_e, vx, vy = sl_init(
        vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, LSETTLS=True
    )

    ######### jstep > 0 ##########
    for jstep in range(0, NITMP):
        
        # Copie des champs
        vx, vy, tracer = backup(
            vx=vx, vy=vy, vx_e=vx_e, vy_e=vy_e, tracer=tracer, tracer_e=tracer_e
        )

        # Estimations
        vx_e, vy_e, tracer_e = sl_xy(
            I=I,
            J=J,
            x=xcr,
            y=ycr,
            vx=vx,
            vy=vy,
            vx_e=vx_e,
            vy_e=vy_e,
            tracer=tracer,
            tracer_e=tracer_e,
            dt=dt,
            dx=dx,
            dy=dy,
            nx=nx,
            ny=ny,
            bcx_kind=bcy_kind,
            bcy_kind=bcy_kind,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax
        )
        
        print(f"Step : {jstep}, Time : {jstep * dt} s")
        print(f"Tracer : min {np.min(tracer)}, max {np.max(tracer)}, mean {np.mean(tracer)}")
        
        imax, jmax = np.unravel_index(np.argmax(tracer, axis=None), tracer.shape)
        print(f"Max of tracer field,  X : {xc[imax]:.2f} Y : {yc[jmax]:.2f}")
        

    #### Print ####
    fig = plot_2D_scalar_field(xcr, ycr, tracer, 25)
    plt.show()
    
    print(f"Tracer : min {np.min(tracer)}, max {np.max(tracer)}, mean {np.mean(tracer)}")
    imax, jmax = np.unravel_index(np.argmax(tracer, axis=None), tracer.shape)
    print(f"Max of tracer field,  X : {xc[imax]:.2f} Y : {yc[jmax]:.2f}")
    
    
    

    
    