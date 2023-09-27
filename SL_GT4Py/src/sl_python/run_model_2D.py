from typing import Tuple
import numpy as np


def interpolate_linear_2d(lx: np.float64, ly: np.float64, ii: int, jj: int, field: np.ndarray):
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
    
    return field[ii: ii + 2, jj: jj + 2] * px * py.T


def interpolate_cubic_2d(lx: np.float64, ly: np.float64, ii: int, jj: int, field: np.ndarray):
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

    # Polynomes de lagrange d'ordre 3
    p_1 = lambda l: 0.5 * l * (l - 1) * (2 - l) / 3
    p0 = lambda l: (1 - l**2) * (1 - l / 2)
    p1 = lambda l: 0.5 * l * (l + 1) * (l - 2)
    p2 = lambda l: 0.5 * l * (l**2 - 1) / 3
    
    px = np.array([p_1(lx), p0(lx), p1(lx), p2(lx)])
    py = np.array([p_1(ly), p0(ly), p1(ly), p2(ly)])
    
    return field[ii - 1: ii + 3, jj - 1: jj + 3] * px * py.T


def lagrangian_search(
    x: np.ndarray,
    y: np.ndarray,
    I: np.ndarray,
    J: np.ndarray,
    v_e: Tuple[np.ndarray],
    v: Tuple[np.ndarray],
    dt: np.float64,
    nx: int,
    ny: int,
    nsiter: int,
    LNESC: bool,
    LSETTLS: bool
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

    vx_e, vy_e = v_e[0], v_e[1]
    vx, vy = v[0], v[1]
    
    # Array declaration
    vx_d = np.zeros((nx, ny))
    vy_d = np.zeros((nx, ny))
    
    vx_d_p = np.zeros((nx, ny))+4
    vy_d_p = np.zeros((nx, ny))

    ################ l = 1 ##############

    # trajectories
    disp_x = dt * vx
    disp_y = dt * vy

    #  r_d
    x_dep = x - disp_x
    y_dep = y - disp_y

    # Indices interpolation
    ii = I - np.floor(disp_x / dx)
    jj = J - np.floor(disp_y / dy)

    # Lambda 
    lx = dx * (disp_x / dx - np.floor(disp_x / dx))
    ly = dy * (disp_y / dy - np.floor(disp_y / dy))
    
    
    for i in range(nx):
        for j in range(ny):
            vx_d[i, j] = interpolate_linear_2d(lx=lx, ly=ly, ii=ii, jj=jj, field=vx)
            vy_d[i, j] = interpolate_linear_2d(lx=lx, ly=ly, ii=ii, jj=jj, field=vy)
            
            vx_d_p[i, j] = interpolate_linear_2d(lx=lx, ly=ly, ii=ii, jj=jj, field=vx_p)
            vy_d_p[i, j] = interpolate_linear_2d(lx=lx, ly=ly, ii=ii, jj=jj, field=vy_p)
    
    ################ l > 1 ##################
    for l in range(1, nsiter):
            
        disp_x = 0.5 * dt * (vx_e + vx_d)
        x_dep = x - disp_x 
        
        disp_y = 0.5 * dt * (vy_e + vy_d)
        y_dep = y - disp_y 
        
        lx = disp_x / dx - np.floor(disp_x / dx)
        ly = disp_y / dy - np.floor(disp_y / dy)
        
        ii = I - np.floor(disp_x / dx)
        jj = J - np.floor(disp_y / dy)

        left_exceed = (x_dep > xmax).astype(np.int8)
        right_exceed = (x_dep < xmin).astype(np.int8)
        
        top_exceed = (y_dep > ymax).astype(np.int8)
        bottom_exceed = (y_dep < ymin).astype(np.int8)
        
        # Periodic x boundaries
        if bool(bcx_kind):
            x_dep = (
                x_dep * (1 - right_exceed)* (1 - left_exceed)
                +  (xmin + x_dep - xmax) * right_exceed
                +  (xmax + x_dep - xmin) * left_exceed
                )
            
            ii = (
                ii * (1 - right_exceed)* (1 - left_exceed)
                +  right_exceed * np.floor((xmin + x_dep - xmax) / dx) 
                +  left_exceed * np.floor((xmax + x_dep - xmin) / dx)
                )
            
            lx = (
                lx * (1 - right_exceed)* (1 - left_exceed)
                +  right_exceed * ((xmin + x_dep - xmax) / dx - np.floor((xmin + x_dep - xmax) / dx))
                +  left_exceed * ((xmax + x_dep - xmin) / dx) - np.floor((xmax + x_dep - xmin) / dx)
            )
           
        # Fixed x boundaries
        else:
            x_dep = (
                x_dep * (1 - right_exceed) * (1 - left_exceed)
                +  xmin * left_exceed
                +  xmax  * right_exceed
                )    
            
            ii = (
                ii * (1 - right_exceed) * (1 - left_exceed)
                +  0 * left_exceed
                +  (nx - 1)  * right_exceed
                ) 
            
            lx = lx * (1 - right_exceed) * (1 - left_exceed)      
        
        # Periodic y boundaries
        if bool(bcy_kind):
            y_dep = (
                y_dep * (1 - top_exceed)* (1 - bottom_exceed)
                +  (ymin + y_dep - ymax) * top_exceed
                +  (ymax + y_dep - ymin) * bottom_exceed
                )
            
            jj = (
                jj * (1 - top_exceed)* (1 - bottom_exceed)
                +  top_exceed * np.floor((ymin + y_dep - ymax) / dy) 
                +  bottom_exceed * np.floor((ymax + y_dep - ymin) / dy)
                )
            
            ly = (
                ly * (1 - top_exceed)* (1 - bottom_exceed)
                +  top_exceed * ((ymin + y_dep - xmax) / dy - np.floor((ymin + y_dep - ymax) / dy))
                +  bottom_exceed * ((ymax + y_dep - ymin) / dy) - np.floor((ymax + y_dep - ymin) / dy)
                )            
           
        # Fixed y boundaries 
        else:
            y_dep = (
                y_dep * (1 - top_exceed) * (1 - bottom_exceed)
                +  ymin * top_exceed
                +  ymax * bottom_exceed
                )  
            
            jj = (
                jj * (1 - top_exceed) * (1 - bottom_exceed)
                +  0 * bottom_exceed
                +  (ny - 1)  * top_exceed
                ) 
            
            ly = ly * (1 - top_exceed) * (1 - bottom_exceed)  
                    
        
        ####### Interpolation for fields ########
        for i in range(nx):
            for j in range(ny):
                vx_e[i, j] = interpolate_cubic_2d(lx[i, j], ly[i, j], ii[i, j], jj[i, j], vx_e)
                vy_e[i, j] = interpolate_cubic_2d(lx[i, j], ly[i, j], ii[i, j], jj[i, j], vy_e)
                
                vx[i, j] = interpolate_cubic_2d(lx[i, j], ly[i, j], ii[i, j], jj[i, j], vx_d)
                vy[i, j] = interpolate_cubic_2d(lx[i, j], ly[i, j], ii[i, j], jj[i, j], vy_d)
                
    return x_dep, y_dep, lx, ly, ii, jj    
    
    

if __name__ == "__main__":
    
    # Init option
    LSETTLS = True
    LNESC = False
    
    # Grid
    xmin, xmax = 0, 100
    ymin, ymax = 0, 100
    zmin, zmax = 0, 100
    nx, ny, nz = 100, 100, 100

    # spacings
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    dz = (zmax - zmin) / nz

    # Boundaries
    # 1 : PERIODIC
    # 0 : FIXED
    bcx_kind = 1
    bcy_kind = 1
    bcz_kind = 0

    # Spacing
    xc = np.linspace(xmin, xmax, nx)
    yc = np.linspace(ymin, ymax, ny)
    zc = np.linspace(zmin, zmax, nz)
    x, y = np.meshgrid(xc, yc)
    
    # Horizontal indexes
    i_indices = np.arange(0, nx)
    j_indices = np.arange(0, ny)

    I, J = np.meshgrid((i_indices, j_indices))
    
    # Initialisation simple
    # Vent uniforme
    U, V = 20, 0, 0

    # Iterations
    NSITER = 5  # Semi lagrangian step
    NITMP = 20
    dt = 10

    # TODO: Revoir l'initilisation
    vx, vy = U * np.ones((nx, ny)), V * np.ones((nx, ny))
    vx_p, vy_p = vx.copy(), vy.copy()
    
    # Champs de vent à t + dt
    vx_e, vy_e = np.zeros((nx, ny))
    

    psi = np.zeros((nx, ny, nz))
    psib = psi.copy()
    
    ########## Premier pas ###########
    # j_iter = 0
    if LNESC:
        vx_e = vx.copy()
        
    elif LSETTLS:
        vx_e = vx.copy()
        vx = 2 * vx - vx_p    
        
    ######### j_iter > 0 ##########   
    for j_iter in range(1, NITMP):
        
        # Backup psi(t - dt)
        psi_b = psi.copy()
    
        # Recherche semi lag
        x_dep, y_dep, lx, ly, ii, jj = lagrangian_search(
            x=x, y=y, I=I, J=J, v_e=(vx_e, vy_e), v=(vx, vy)
        )

        for i in range(nx):
            for j in range(ny):
                psi[i, j] = interpolate_cubic_2d(lx=lx, ly=ly, ii=ii, jj=jj, field=psi_b)
                
                # Interpolation vitesse ?
