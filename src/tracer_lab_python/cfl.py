def cfl_1d(u: float, dx: float, dt: float):
    """Computes directional CFD

    Args:
        u (float): contravariant velocity
        dx (float): spacing on axis
        dt (float): time step

    Returns:
        float: Courant number for a given direction
    """
    return u * dt / dx