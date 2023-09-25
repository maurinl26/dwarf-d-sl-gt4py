from typing import Tuple
import numpy as np


def interpolate_linear(l: float, i: int, psi: np.ndarray):
    return None


def interpolate_cubic(l: float, i: int, psi: np.ndarray):
    return None


def distance_to_grid(x: Tuple[np.ndarray], grid: Tuple[np.ndarray]):
    i, j = 0, 0
    lx, ly = 0, 0

    return i, j, lx, ly


def lagrangian_research(
    grid: Tuple[np.ndarray],
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    v: Tuple[np.ndarray],
    v_p: Tuple[np.ndarray],
    dt: np.float64,
    nsiter = 10
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

    # l = 1
    X_dep[0] = x - dt * v[0]
    X_dep[1] = y - dt * v[1]
    X_dep[3] = z - dt * v[3]

    i, j, k, lx, ly, lz = distance_to_grid(x_dep, grid)

    v_d = interpolate_linear(l, i, v)
    v_d_p = interpolate_linear(l, i, v_p)

    # l > 1
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

    # Spacing
    xc = np.linspace(xmin, xmax, nx)
    yc = np.linspace(ymin, ymax, ny)
    grid = np.meshgrid(xc, yc)

    # Iterations
    NSITER = 5  # Semi lagrangian step
    NITMP = 20
    dt = 10

    vx, vy, vz = np.zeros((nx, ny, nz)), np.zeros((nx, ny, nz)), np.zeros((nx, ny, nz))

    psi = np.zeros((nx, ny, nz))
    psib = np.zeros((nx, ny, nz))
    
    # Initialization
    ## TODO : work on init cases

    for j_iter in range(0, NITMP):
        
        # Backup
        psi_b = psi.copy()

        # Recherche semi lag
        X_dep = lagrangian_research()
        ii, jj, kk, lx, ly, lz = distance_to_grid(X_dep, grid)

        # Update champs psi
        psi_x = interpolate_cubic(lx, ii, psi)
        psi_y = interpolate_cubic(ly, jj, psi_x)
        psi = interpolate_cubic(lz, kk, psi_y)
        
        