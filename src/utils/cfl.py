import numpy as np


def cfl_1d(
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    dx: float,
    dy: float,
    dz: float,
    dt: float
):
    """Maximum 1D CFL

    Args:
        u (np.ndarray): horizontal velocity on x
        v (np.ndarray): horizontal velocity on y 
        w (np.ndarray): vertical velocity on z
        dx (float): x spacing
        dy (float): y spacing
        dz (float): z spacing
        dt (float): time step

    Returns:
        float: _description_
    """
    return max(np.max((u * dt) / dx, (v * dt) / dy, (w * dt) / dz))
    