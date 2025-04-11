import numpy as np
import matplotlib.pyplot as plt

@dace.program
def diagnostic_lipschitz(
        u: np.ndarray,
        v: np.ndarray,
        du_dx: np.ndarray,
        dv_dx: np.ndarray,
        du_dy: np.ndarray,
        dv_dy: np.ndarray,
        dx: float,
        dy: float,
        dth: float):
    """Diagnostic for Semi-Lagrangian Research stability

    Args:
        u (np.ndarray): velocity on x 
        v (np.ndarray): velocity on y
        dx (float): spacing on x
        dy (float): spacing on y
        dth (float): half time step

    Returns:
        float: lipschitz condition on stability
    """

    du_dx = c2_x(u, du_dx, dx, origin=(1, 0, 0))
    dv_dx = c2_x(v, dv_dx, dx, origin=(1, 0, 0))

    du_dy = c2_y(u, du_dy, dy, origin=(0, 1, 0))
    dv_dy = c2_y(v, dv_dy, dy, origin=(0, 1, 0))

    
    return dth * np.max(np.maximum(du_dx, du_dy), np.maximum(dv_dx, dv_dy))

