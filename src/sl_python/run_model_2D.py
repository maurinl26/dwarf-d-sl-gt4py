from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt

import sys

sys.path.append("/home/maurinl/FVM_GT4Py_slim/src")
sys.path.append("/home/maurinl/sl_gt4py/SL_GT4Py/src")

from sl_python.plot import plot_2D_scalar_field


def interpolate_linear_2d(
    lx: np.float64, ly: np.float64, ii: int, jj: int, field: np.ndarray
):
    """Interpolate sequentially on x axis and on y axis

    Args:
        l (np.float64): _description_
        ii (int): _description_
        field (np.ndarray): _description_

    Returns:
        _type_: _description_
    """
    p0 = lambda l: 1 - l
    p1 = lambda l: l

    px = np.array([p0(lx), p1(lx)])
    py = np.array([p0(ly), p1(ly)])

    return field[ii : ii + 2, jj : jj + 2] * px * py.T


def interpolate_cubic_2d(
    lx: np.float64, ly: np.float64, ii: int, jj: int, field: np.ndarray
) -> np.float64:
    """Interpolate sequentially on x axis and on y axis.

    Args:
        lx (np.float64): _description_
        ly (np.float64): _description_
        ii (int): _descriptionsaliur
        jj (int): _description_
        field (np.ndarray): _description_

    Returns:
        _type_: _description_
    """

    # Padding on interpolation field
    # Fixed
    if bcx_kind == 0:
        padded_field = np.pad(field, (1, 2), "edge")
        ii += 1
        jj += 1

    # Periodic
    else:
        padded_field = np.pad(field, (1, 2), "wrap")
        ii += 1
        jj += 1

    # Polynomes de lagrange d'ordre 3
    p_1 = lambda l: 0.5 * l * (l - 1) * (2 - l) / 3
    p0 = lambda l: (1 - l**2) * (1 - l / 2)
    p1 = lambda l: 0.5 * l * (l + 1) * (l - 2)
    p2 = lambda l: 0.5 * l * (l**2 - 1) / 3

    px = np.array([p_1(lx), p0(lx), p1(lx), p2(lx)])
    py = np.array([p_1(ly), p0(ly), p1(ly), p2(ly)])

    psi_hat = np.dot(np.matmul(padded_field[ii - 1 : ii + 3, jj - 1 : jj + 3], px), py)

    return psi_hat


def boundaries(
    bc_kind: int,
    points: np.ndarray,
    indices: np.ndarray,
    l: np.ndarray,
    xmin: np.float64,
    xmax: np.float64,
    nx: int,
):
    """Apply boundary conditions on field.

    Args:
        bc_kind (int): _description_
        field (np.ndarray): _description_
        min (np.float): _description_
        max (np.float): _description_
        nx (int): _description_
    """
    left_exceed = points > xmax
    right_exceed = points < xmin

    left_recovering = xmin + points - xmax
    right_recovering = xmax + points - xmin

    # Periodic boundaries
    if bool(bc_kind):
        points = (
            points * (1 - right_exceed) * (1 - left_exceed)
            + right_exceed * left_recovering
            + left_exceed * right_recovering
        )

        indices = (
            indices * (1 - right_exceed) * (1 - left_exceed)
            + right_exceed * np.floor(left_recovering / dx)
            + left_exceed * np.floor(right_recovering / dx)
        )

        l = (
            l * (1 - right_exceed) * (1 - left_exceed)
            + right_exceed * (left_recovering / dx - np.floor(left_recovering / dx))
            + left_exceed * (right_recovering / dx)
            - np.floor(right_recovering / dx)
        )

    # Fixed boundaries
    else:
        points = (
            points * (1 - right_exceed) * (1 - left_exceed)
            + xmin * left_exceed
            + xmax * right_exceed
        )

        indices = (
            indices * (1 - right_exceed) * (1 - left_exceed)
            + 0 * left_exceed
            + (nx - 1) * right_exceed
        )

        l = l * (1 - right_exceed) * (1 - left_exceed)

    return indices.astype(np.int8), l


def lagrangian_search(
    I: np.ndarray,
    J: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    dt: np.float64,
    dx: np.float64,
    dy: np.float64,
    nx: int,
    ny: int,
    bcx_kind: int,
    bcy_kind: int,
    xmin: np.float64,
    xmax: np.float64,
    ymin: np.float64,
    ymax: np.float64,
    nsiter: int = 10,
) -> Tuple[np.ndarray]:
    """Research departure point for a given grid and velocity field.
    Terminates on nsiter iterations.

    Args:
        x (np.ndarray): grid of arrival points
        v (np.ndarray): velocity fields
        nsiter (int, optional): number of iterations. Defaults to 10.

    Returns:
        np.ndarray: departure point
    """

    # Array declaration
    ################ l > 1 ##################
    for l in range(0, nsiter):
        if l == 0:
            disp_x = dt * vx
            disp_y = dt * vy

        else:
            disp_x = 0.5 * dt * (vx_e + vx)
            disp_y = 0.5 * dt * (vy_e + vy)

        x_dep = x - disp_x
        y_dep = y - disp_y

        lx = disp_x / dx - np.floor(disp_x / dx)
        ly = disp_y / dy - np.floor(disp_y / dy)

        i_d = I - np.floor(disp_x / dx)
        j_d = J - np.floor(disp_y / dy)

        i_d, lx = boundaries(bcx_kind, x_dep, i_d, lx, xmin, xmax, nx)
        j_d, ly = boundaries(bcy_kind, y_dep, j_d, ly, ymin, ymax, ny)

        ####### Interpolation for fields ########
        for i in range(nx):
            for j in range(ny):
                # interpolation en r_d(l) -> i_d(l), j_d(l)
                vx_e[i, j] = interpolate_cubic_2d(
                    lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vx
                )
                vy_e[i, j] = interpolate_cubic_2d(
                    lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vy
                )

                # interpolation en r_d(l) -> i_d(l), j_d(l)
                vx[i, j] = interpolate_cubic_2d(
                    lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vx_p
                )
                vy[i, j] = interpolate_cubic_2d(
                    lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vy_p
                )

    return lx, ly, i_d.astype(np.int8), j_d.astype(np.int8), vx_e, vy_e


def sl_init(
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_p: np.ndarray,
    vy_p: np.ndarray,
    LSETTLS: bool = True,
) -> Tuple[np.ndarray]:
    # LSETTLS
    if LSETTLS:
        vx_e = vx.copy()
        vy_e = vy.copy()

        vx = 2 * vx - vx_p
        vy = 2 * vy - vy_p

    # LNESC
    else:
        vx_e = vx.copy()
        vy_e = vy.copy()

    return vx, vy, vx_e, vy_e


def sl_xy(
    I,
    J,
    x,
    y,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    tracer: np.ndarray,
    tracer_e: np.ndarray,
    dt: np.float64,
    dx: np.float64,
    dy: np.float64,
    nx: int,
    ny: int,
    bcx_kind: int,
    bcy_kind: int,
    xmin: np.float64,
    xmax: np.float64,
    ymin: np.float64,
    ymax: np.float64
):
    # Recherche semi lag
    lx_d, ly_d, i_d, j_d, vx_e, vy_e = lagrangian_search(
        x=x,
        y=y,
        I=I,
        J=J,
        vx_e=vx_e,
        vy_e=vy_e,
        vx=vx,
        vy=vy,
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

    # Interpolate
    for i in range(nx):
        for j in range(ny):
            # Interpolate tracer in T(r, t) = T(r_d, t - dt)

            tracer_e[i, j] = interpolate_cubic_2d(
                lx=lx_d[i, j], ly=ly_d[i, j], ii=i_d[i, j], jj=j_d[i, j], field=tracer
            )

    return vx_e, vy_e, tracer_e


def backup(
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    tracer: np.ndarray,
    tracer_e: np.ndarray,
) -> Tuple[np.ndarray]:
    """Remap ebauche fields at t into fields at tau = t + dt

    Args:
        vx (np.ndarray): x velocity
        vy (np.ndarray): y velocity
        vx_e (np.ndarray): ebauche vx
        vy_e (np.ndarray): ebauche vy
        tracer (np.ndarray): tracer field
        tracer_e (np.ndarray): ebauche at t + dt for tracer field

    Returns:
        _type_: _description_
    """

    # Copie des champs
    tracer = tracer_e.copy()
    vx = vx_e.copy()
    vy = vy_e.copy()

    return vx, vy, tracer


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
    U, V = 10, 10

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
            bcx_kind=bcx_kind,
            bcy_kind=bcy_kind,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax
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
