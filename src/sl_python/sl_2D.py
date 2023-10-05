from typing import Tuple
import numpy as np
import itertools

from config import Config

def boundaries(
    bc_kind: int,
    indices: np.ndarray,
    n: int,
):
    """Apply boundary conditions on field.

    Args:
        bc_kind (int): _description_
        field (np.ndarray): _description_
        min (np.float): _description_
        max (np.float): _description_
        nx (int): _description_
    """
    left_exceed = indices >= n
    right_exceed = indices < 0

    # Periodic boundaries
    if bool(bc_kind):
        indices %= n

    # Fixed boundaries
    else:
        indices = (
            indices * (1 - right_exceed) * (1 - left_exceed)
            + 0 * left_exceed
            + (n - 1) * right_exceed
        )

    return indices


def lagrangian_search(
    config: Config,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    interpolation_function: callable,
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

    vx_tmp = vx.copy()
    vy_tmp = vy.copy()

    # Array declaration
    for l in range(0, nsiter):
        traj_x = config.dth * (vx_e + vx_tmp)
        traj_y = config.dth * (vy_e + vy_tmp)

        lx = traj_x / config.dx - np.floor(traj_x / config.dx)
        ly = traj_y / config.dy - np.floor(traj_y / config.dy)

        i_d = np.floor(config.I - traj_x / config.dx).astype(np.int64)
        j_d = np.floor(config.J - traj_y / config.dy).astype(np.int64)

        i_d = boundaries(config.bcx_kind, i_d, config.nx)
        j_d = boundaries(config.bcy_kind, j_d, config.ny)

        ####### Interpolation for fields ########
        for i, j in itertools.product(range(config.nx), range(config.ny)):
            # interpolation en r_d(l) -> i_d(l), j_d(l)
            vx_tmp[i, j] = interpolation_function(
                lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vx, config.bcx_kind
            )
            vy_tmp[i, j] = interpolation_function(
                lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vy, config.bcy_kind
            )

    return lx, ly, i_d, j_d


def sl_init(
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_p: np.ndarray,
    vy_p: np.ndarray,
    lsettls: bool = True,
) -> Tuple[np.ndarray]:
    # LSETTLS
    if lsettls:
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
    config: Config,
    vx: np.ndarray,
    vy: np.ndarray,
    vx_e: np.ndarray,
    vy_e: np.ndarray,
    tracer: np.ndarray,
    tracer_e: np.ndarray,
    interpolation_function: callable,
    nsiter: int,
):
    # Recherche semi lag
    lx_d, ly_d, i_d, j_d = lagrangian_search(
        config=config,
        vx_e=vx_e,
        vy_e=vy_e,
        vx=vx,
        vy=vy,
        interpolation_function=interpolation_function,
        nsiter=nsiter,
    )

    # Interpolate
    for i, j in itertools.product(range(config.nx), range(config.ny)):
        # Interpolate tracer in T(r, t) = T(r_d, t - dt)
        tracer_e[i, j] = interpolation_function(
            lx=lx_d[i, j],
            ly=ly_d[i, j],
            ii=i_d[i, j],
            jj=j_d[i, j],
            field=tracer,
            bc_kind=config.bcx_kind,
        )

        # Ebauche vitesse
        vx_e[i, j] = interpolation_function(
            lx=lx_d[i, j],
            ly=ly_d[i, j],
            ii=i_d[i, j],
            jj=j_d[i, j],
            field=vx,
            bc_kind=config.bcx_kind,
        )
        vy_e[i, j] = interpolation_function(
            lx=lx_d[i, j],
            ly=ly_d[i, j],
            ii=i_d[i, j],
            jj=j_d[i, j],
            field=vy,
            bc_kind=config.bcy_kind,
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
