from typing import Tuple
import numpy as np


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
    lx: np.float64, ly: np.float64, ii: int, jj: int, field: np.ndarray, bc_kind: int
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
    if bc_kind == 0:
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
    dx: float
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
            traj_x = 0.5 * dt * (vx_e + vx)
            traj_y = 0.5 * dt * (vy_e + vy)

        else:
            traj_x = 0.5 * dt * (vx_e + vx)
            traj_y = 0.5 * dt * (vy_e + vy)

        x_dep = x - traj_x
        y_dep = y - traj_y

        lx = traj_x / dx - np.floor(traj_x / dx)
        ly = traj_y / dy - np.floor(traj_y / dy)

        i_d = I - np.floor(traj_x / dx)
        j_d = J - np.floor(traj_y / dy)

        i_d, lx = boundaries(bcx_kind, x_dep, i_d, lx, xmin, xmax, nx, dx)
        j_d, ly = boundaries(bcy_kind, y_dep, j_d, ly, ymin, ymax, ny, dy)

        ####### Interpolation for fields ########
        vx_a = np.zeros((nx, ny))
        vy_a = np.zeros((nx, ny))
        for i in range(nx):
            for j in range(ny):
            
                # interpolation en r_d(l) -> i_d(l), j_d(l)
                vx_a[i, j] = interpolate_cubic_2d(
                    lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vx, bcx_kind
                )
                vy_a[i, j] = interpolate_cubic_2d(
                    lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vy, bcx_kind
                )

    return lx, ly, i_d.astype(np.int8), j_d.astype(np.int8)

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
    lx_d, ly_d, i_d, j_d = lagrangian_search(
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
                lx=lx_d[i, j], ly=ly_d[i, j], ii=i_d[i, j], jj=j_d[i, j], field=tracer, bc_kind=bcx_kind
            )
            
            # Ebauche vitesse
            vx_e[i, j] = interpolate_cubic_2d(
                lx=lx_d[i, j], ly=ly_d[i, j], ii=i_d[i, j], jj=j_d[i, j], field=vx, bc_kind=bcx_kind
            )
            vy_e[i, j] = interpolate_cubic_2d(
                lx=lx_d[i, j], ly=ly_d[i, j], ii=i_d[i, j], jj=j_d[i, j], field=vy, bc_kind=bcx_kind
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

