import numpy as np
from typing import Optional, Tuple

def slag_init(
    vx: np.ndarray,
    vy: np.ndarray,
    vx_p: Optional[np.ndarray] = None,
    vy_p: Optional[np.ndarray] = None,
    lsettls: bool = False,
) -> Tuple[np.ndarray]:
    """Initialize draft velocities with either
    LSETTLS method : 2 fields for velocity (at t and t - dt)
    LNESN (not LSETTLS) method : 1 field for velocity (at t)

    Args:
        vx_e (np.ndarray): outlined velocity at t + dt on x
        vy_e (np.ndarray): outlined velocity at t + dt on y
        vx (np.ndarray): velocity at t on x
        vy (np.ndarray): velocity at t on y
        vx_p (np.ndarray): velocity at t - dt on x
        vy_p (np.ndarray): velcoity at t - dt on y
        lsettls (bool, optional): LSETTLS or LNESC. Defaults to True.

    Returns:
        Tuple[np.ndarray]: velocities at t and t + dt
    """
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

    return vx_e, vy_e