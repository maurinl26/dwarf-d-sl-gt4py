from typing import Tuple
import numpy as np
import itertools


def interpolate_linear_2d(
    lx: np.float64, ly: np.float64, ii: np.int64, jj: np.int64, field: np.ndarray, bc_kind: int
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
    
    if bc_kind == 0:
        padded_field = np.pad(field, (1, 2), "edge")
        
    # Periodic
    else:
        padded_field = np.pad(field, (1, 2), "wrap")

    return np.dot(np.matmul(padded_field[ii + 1 : ii + 3, jj + 1 : jj + 3], px), py)



def interpolate_cubic_2d(
    lx: np.float64, ly: np.float64, ii: np.int64, jj: np.int64, field: np.ndarray, bc_kind: int
) -> np.float64:
    """Interpolate sequentially on x axis and on y axis.

    Args:
        lx (np.float64): _description_
        ly (np.float64): _description_
        ii (int): _description_
        jj (int): _description_
        field (np.ndarray): _description_

    Returns:
        _type_: _description_
    """
    
    # Padding on interpolation field
    # Fixed
    if bc_kind == 0:
        padded_field = np.pad(field, (1, 2), "edge")
        
    # Periodic
    else:
        padded_field = np.pad(field, (1, 2), "wrap")

    # Polynomes de lagrange d'ordre 3
    p_1 = lambda l: (1/6) * l * (l - 1) * (2 - l) 
    p0 = lambda l: (1 - l**2) * (1 - l / 2)
    p1 = lambda l: (1/2) * l * (l + 1) * (2 - l)
    p2 = lambda l: (1/6) * l * (l**2 - 1) 

    px = np.array([p_1(lx), p0(lx), p1(lx), p2(lx)])
    py = np.array([p_1(ly), p0(ly), p1(ly), p2(ly)])
    
    assert(round(sum(px)) == 1)
    assert(round(sum(py)) == 1)

    return np.dot(np.matmul(padded_field[ii : ii + 4, jj : jj + 4], px), py)


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
    left_exceed = (indices >= n)
    right_exceed = (indices < 0)

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
    I: np.ndarray,
    J: np.ndarray,
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
        
        traj_x = (dt/2) * (vx_e + vx_tmp)
        traj_y = (dt/2) * (vy_e + vy_tmp)       
                
        lx = traj_x / dx - np.floor(traj_x / dx) 
        ly = traj_y / dy - np.floor(traj_y / dy)
        
        i_d = np.floor(I - traj_x / dx).astype(np.int64)
        j_d = np.floor(J - traj_y / dy).astype(np.int64)
                    
        i_d = boundaries(bcx_kind, i_d, nx)
        j_d = boundaries(bcy_kind, j_d, ny)

        ####### Interpolation for fields ########
        for i, j in itertools.product(range(nx), range(ny)):
            
                # interpolation en r_d(l) -> i_d(l), j_d(l)
                vx_tmp[i, j] = interpolate_cubic_2d(
                    lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vx, bcx_kind
                )
                vy_tmp[i, j] = interpolate_cubic_2d(
                    lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vy, bcy_kind
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
    I: np.ndarray,
    J: np.ndarray,
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
    nsiter: int
):
    
    # Recherche semi lag
    lx_d, ly_d, i_d, j_d = lagrangian_search(
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
        nsiter=nsiter
    )
    
    # Interpolate
    for i, j in itertools.product(range(nx), range(ny)):
            
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

