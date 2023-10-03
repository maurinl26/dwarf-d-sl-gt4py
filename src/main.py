""" Taken from FVM_GT4Py_slim blossey_xy """
import sys

from matplotlib import pyplot as plt
import numpy as np


sys.path.append("/home/maurinl/FVM_GT4Py_slim/src")
sys.path.append("/home/maurinl/SL_GT4Py/SL8GT4Py/src")

from sl_python.sl_2D import sl_init, sl_xy, backup
from sl_python.plot import plot_2D_scalar_field


def main(lsettls: bool, nsiter: int, nitmp: int):
    lnesc = not lsettls
    dt = 1

    # Grid
    xmin, xmax = 0, 123
    ymin, ymax = 0, 123
    nx, ny = 123, 123

    # spacings
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny

    # Boundaries
    # 1 : PERIODIC
    # 0 : FIXED
    bcx_kind = 1
    bcy_kind = 1

    # Spacing
    xc = np.linspace(xmin, xmax, nx)
    yc = np.linspace(ymin, ymax, ny)
    xcr, ycr = np.meshgrid(xc, yc)

    # Horizontal indexes
    i_indices = np.arange(0, nx).astype(np.int64)
    j_indices = np.arange(0, ny).astype(np.int64)
    I, J = np.meshgrid(i_indices, j_indices)

    # Initialisation simple
    # Vent uniforme
    U, V = 20.5, 20.5

    ############## Declaration des champs #############
    vx, vy = (U * np.ones((nx, ny)), V * np.ones((nx, ny)))

    # Champs de vent à t - dt
    vx_p, vy_p = vx.copy(), vy.copy()

    # Champs de vent à t + dt
    vx_e, vy_e = (np.zeros((nx, ny)), np.zeros((nx, ny)))

    # Tracer Gaussien
    T = 20
    rsq = (xcr - 20) ** 2 + (ycr - 20) ** 2

    tracer = T * np.exp(-(rsq / (2 * 10)))
    tracer_e = tracer.copy()

    ###### Plot Initial Fields #####
    fig = plot_2D_scalar_field(xcr, ycr, tracer, 25)
    plt.show()
    print(f"Step : {0}, Time : {0} s")
    print(
        f"Tracer : min {np.min(tracer)}, max {np.max(tracer)}, mean {np.mean(tracer)}"
    )
    imax, jmax = np.unravel_index(np.argmax(tracer, axis=None), tracer.shape)
    print(f"Max of tracer field,  X : {xc[imax]:.2f} Y : {yc[jmax]:.2f}")

    ########### Advection ##########
    ######### Premier pas ##########

    ######### jstep = 0 ##########

    # Initialisation vitesses
    vx_e, vy_e, vx, vy = sl_init(
        vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, lsettls=lsettls
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
            bcx_kind=bcx_kind,
            bcy_kind=bcy_kind,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
        )

        print(f"Step : {jstep}, Time : {jstep * dt} s")
        print(
            f"Tracer : min {np.min(tracer)}, max {np.max(tracer)}, mean {np.mean(tracer)}"
        )

        imax, jmax = np.unravel_index(np.argmax(tracer, axis=None), tracer.shape)
        print(f"Max of tracer field,  X : {xc[imax]:.2f} Y : {yc[jmax]:.2f}")

    #### Print ####
    fig = plot_2D_scalar_field(xcr, ycr, tracer, 25)
    plt.show()

    print(
        f"Tracer : min {np.min(tracer)}, max {np.max(tracer)}, mean {np.mean(tracer)}"
    )
    imax, jmax = np.unravel_index(np.argmax(tracer, axis=None), tracer.shape)
    print(f"Max of tracer field,  X : {xc[imax]:.2f} Y : {yc[jmax]:.2f}")


if __name__ == "__main__":
    LSETTLS = False
    NITMP = 7
    NSITER = 5

    main(LSETTLS, NSITER, NITMP)
