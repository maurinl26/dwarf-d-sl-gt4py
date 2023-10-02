""" Run semi lagrangian advection in 2D xy or yz with padded array to process boundaries or interpolations"""

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

    # Polynomes de lagrange d'ordre 3
    p_1 = lambda l: 0.5 * l * (l - 1) * (2 - l) / 3
    p0 = lambda l: (1 - l**2) * (1 - l / 2)
    p1 = lambda l: 0.5 * l * (l + 1) * (l - 2)
    p2 = lambda l: 0.5 * l * (l**2 - 1) / 3

    px = np.array([p_1(lx), p0(lx), p1(lx), p2(lx)])
    py = np.array([p_1(ly), p0(ly), p1(ly), p2(ly)])      
    
    psi_hat = np.dot(np.matmul(field[ii - 1 : ii + 3, jj - 1 : jj + 3], px), py)
    
    return psi_hat

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
    nx: int,
    ny: int,
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

        lx = disp_x / dx - np.floor(disp_x / dx)
        ly = disp_y / dy - np.floor(disp_y / dy)

        i_d = I#_pad - np.floor(disp_x / dx).astype(np.int64)
        j_d = J#_pad - np.floor(disp_y / dy).astype(np.int64)

        ####### Interpolation for fields ########
        for i in range(nx):
            for j in range(ny):
                                 
                # interpolation en r_d(l) -> i_d(l), j_d(l)
                vx_e[i, j] = interpolate_cubic_2d(
                    lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vx#_padded
                )
                vy_e[i, j] = interpolate_cubic_2d(
                    lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vy#_padded
                )

                # interpolation en r_d(l) -> i_d(l), j_d(l)
                vx[i, j] = interpolate_cubic_2d(
                    lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vx_p#_padded
                )
                vy[i, j] = interpolate_cubic_2d(
                    lx[i, j], ly[i, j], i_d[i, j], j_d[i, j], vy_p#_padded
                )

    return lx, ly, i_d, j_d, vx_e, vy_e


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
    nx: int,
    ny: int
):
    # Recherche semi lag
    lx_d, ly_d, i_d, j_d, vx_e, vy_e = lagrangian_search(
        x=x, y=y, I=I, J=J, vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, dt=dt, nx=nx, ny=ny
    )

    # Interpolate
    for i in range(nx):
        for j in range(ny):
            # Interpolate tracer in T(r, t) = T(r_d, t - dt)
            
            tracer_e[i, j] = interpolate_cubic_2d(
                lx=lx_d[i, j], ly=ly_d[i, j], ii=i_d[i, j], jj=j_d[i, j], field=tracer.#_padded
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

# For fields with dx = dy, bcx_kind = bcy_kind
def padding(field: np.ndarray, I, J, dx, dy, dt, Vmax, bcx_kind):
    
    pad_width_x = int(np.floor(Vmax * dt / dx) + 1)
    pad_width_y = int(np.floor(Vmax * dt / dy) + 1)
    
    I_pad = I + pad_width_x
    J_pad = J + pad_width_y
        
    padded_field = np.pad(field, ((pad_width_x, pad_width_x),(pad_width_y, pad_width_y)), "edge") 

    return I_pad, J_pad, padded_field

if __name__ == "__main__":
    # Init option
    LSETTLS = True
    LNESC = not LSETTLS  # Pour info

    # Iterations
    NSITER = 5  # Semi lagrangian step
    NITMP = 20
    dt = 1
    
    # 
    V_max = 20

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
    x, y = np.meshgrid(xc, yc)

    # Horizontal indexes
    i_indices = np.arange(0, nx).astype(np.int8)
    j_indices = np.arange(0, ny).astype(np.int8)
    I, J = np.meshgrid(i_indices, j_indices)

    # Initialisation simple
    # Vent uniforme
    U, V = 20, 0

    ############## Declaration des champs #############
    vx, vy = (U * np.ones((nx, ny)), V * np.ones((nx, ny)))

    # Champs de vent à t - dt
    vx_p, vy_p = vx.copy(), vy.copy()

    # Champs de vent à t + dt
    vx_e, vy_e = (np.zeros((nx, ny)), np.zeros((nx, ny)))

    # Tracer Blossey / Duran
    T = 300
    tracer = T * np.ones((nx, ny))
    tracer_e = np.zeros((nx, ny))

    ############# Advection
    ########## Premier pas ###########
    # jstep = 0

    # Initialisation vitesses
    vx_e, vy_e, vx, vy = sl_init(
        vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, LSETTLS=True
    )

    ######### j_iter > 0 ##########
    for jstep in range(0, NITMP):
        
        
        # Copie des champs
        vx, vy, tracer = backup(
            vx=vx,
            vy=vy,
            vx_e=vx_e,
            vy_e=vy_e,
            tracer=tracer,
            tracer_e=tracer_e
        )
        
        # Padding 
        _, _, x_padded = padding(x, I, J, dx, dy, dt, V_max, bcx_kind)
        _, _, y_padded = padding(x, I, J, dx, dy, dt, V_max, bcx_kind)
        _, _, vx_padded = padding(vx, I, J, dx, dy, dt, V_max, bcx_kind)
        _, _, vy_padded = padding(vy, I, J, dx, dy, dt, V_max, bcx_kind)
        
        I_pad, J_pad, tracer_padded = padding(tracer, I, J, dx, dy, dt, V_max, bcx_kind)

        # Estimations
        vx_e, vy_e, tracer_e = sl_xy(
            I=I_pad,
            J=J_pad,
            x=x,
            y=y,
            vx=vx,
            vy=vy,
            vx_e=vx_e,
            vy_e=vy_e,
            tracer=tracer_padded,
            tracer_e=tracer_e,
            dt=dt,
            nx=nx,
            ny=ny
        )
        
        print(f"vx_e, shape : {vx_e.shape}")
        print(f"vy_e, shape : {vy_e.shape}")
        print(f"tracer_e, shape : {tracer_e.shape}")