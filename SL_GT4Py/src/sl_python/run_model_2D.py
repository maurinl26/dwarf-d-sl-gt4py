from typing import Tuple
import numpy as np


def interpolate_linear_1d(l: np.float64, ii: int, psi: np.ndarray):
    p0 = lambda x: x
    p1 = lambda x: x
    
    
    return None


def interpolate_cubic_2d(lx: np.float64, ly: np.float64, ii: int, jj: int, field: np.ndarray):
    # Polynomes de lagrange d'ordre 3
    p_1 = lambda l: 0.5 * l * (l - 1) * (2 - l) / 3
    p0 = lambda l: (1 - l**2) * (1 - l / 2)
    p1 = lambda l: 0.5 * l * (l + 1) * (l - 2)
    p2 = lambda l: 0.5 * l * (l**2 - 1) / 3
    
    px = np.array([p_1(lx), p0(lx), p1(lx), p2(lx)])
    py = np.array([p_1(ly), p0(ly), p1(ly), p2(ly)])
    
    return field[ii - 1: ii + 3, jj - 1: jj + 3] * px * py.T

    # Greedy interpolation
    # Interp sur X
    interp_x = (
        p_1(lx) * field[ii - 1, jj - 1 : jj + 3]
        + p0(lx) * field[ii, jj - 1 : jj + 3]
        + p1(lx) * field[ii + 1, jj - 1 : jj + 3]
        + p2(lx) * field[ii + 2, jj - 1 : jj + 3]
    )
        
    # Interp sur Y
    interp_y = (
        p_1(ly) * interp_x[jj - 1]
        + p0(ly) * interp_x[jj]
        + p1(ly) * interp_x[jj + 1]
        + p2(ly) * interp_x[jj + 2]
    )
    
    return interp_y


def interpolate_cubic(distance: Tuple[float], indices: Tuple[int], psi: np.ndarray):
    return None


def lagrangian_research(
    grid: Tuple[np.ndarray],
    x: np.ndarray,
    y: np.ndarray,
    v: Tuple[np.ndarray],
    v_p: Tuple[np.ndarray],
    dt: np.float64,
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

    vx, vy, vz = v[0], v[1], v[2]

    # Horizontal indexes
    i_indices = np.arange(0, nx)
    j_indices = np.arange(0, ny)

    I, J = np.meshgrid((i_indices, j_indices))

    ################ l = 1 ##############

    # trajectories
    disp_x = dt * vx
    disp_y = dt * vy

    x_dep = x - disp_x
    y_dep = y - disp_y

    # Indices interpolation
    ii = I - np.floor(dt * v[0] / dx)
    jj = J - np.floor(dt * v[1] / dy)
    # kk = k - np.floor(dt * v[2] / dz)

    # Lambda interpolation
    lx = dx * (disp_x / dx - np.floor(disp_x / dx))
    ly = dy * (disp_y / dy - np.floor(disp_x / dy))

    # X_dep[2] = z - dt * v[2]

    lx = dx * (dt * v[0] / dx - np.floor(dt * v[0] / dx))
    ly = dy * (dt * v[0] / dy - np.floor(dt * v[0] / dy))
    # lz = dz * (dt * v[0] / dz - np.floor(dt * v[0] / dz))

    psi_b = psi.copy()
    for i in range(nx):
        for j in range(ny):
            psi_b = 0
            #

    v_d = interpolate_linear(distance=(lx, ly), indices=(ii, jj), field=v)
    v_d_p = interpolate_linear(l, i, v_p)

    ################ l > 1 ##################
    for l in range(2, nsiter):
        x_dep = x - 0.5 * dt * (v + 2 * v_d - v_d_p)
        i, j, lx, ly = distance_to_grid(x_dep, grid)
        x_dep_b = x_dep.copy()

    return x_dep


if __name__ == "__main__":
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
    # To be used in algorithms
    bcx_kind = 1
    bcy_kind = 1
    bcz_kind = 0

    # Spacing
    xc = np.linspace(xmin, xmax, nx)
    yc = np.linspace(ymin, ymax, ny)
    zc = np.linspace(zmin, zmax, nz)
    grid = np.meshgrid(xc, yc)

    # Iterations
    NSITER = 5  # Semi lagrangian step
    NITMP = 20
    dt = 10

    vx, vy, vz = np.zeros((nx, ny, nz)), np.zeros((nx, ny, nz)), np.zeros((nx, ny, nz))
    vx_p, vy_p, vz_p = vx.copy(), vy.copy(), vz.copy()

    psi = np.zeros((nx, ny, nz))
    psib = psi.copy()

    # Initialization
    ## TODO : work on init cases
    for j_iter in range(0, NITMP):
        # Loops over I, J, K to be set

        # Backup
        psi_b = psi.copy()

        # Recherche semi lag
        X_dep = lagrangian_research(
            grid=grid, x=xc, y=yc, z=zc, v=(vx, vy, vz), v_p=(vx_p, vy_p, vz_p)
        )

        X_dep[0] = x - dt * v[0]
        X_dep[1] = y - dt * v[1]
        X_dep[2] = z - dt * v[2]

        ii = i - np.floor(dt * v[0] / dx)
        jj = j - np.floor(dt * v[1] / dy)
        kk = k - np.floor(dt * v[2] / dz)

        lx = dx * (dt * v[0] / dx - np.floor(dt * v[0] / dx))
        ly = dy * (dt * v[0] / dy - np.floor(dt * v[0] / dy))
        lz = dz * (dt * v[0] / dz - np.floor(dt * v[0] / dz))

        # Update champs psi
        psi = interpolate_cubic(indices=(ii, jj, kk), distances=(lx, ly, lz), field=psi)
