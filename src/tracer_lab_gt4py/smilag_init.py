import logging

from gt4py.cartesian.gtscript import Field, stencil, computation, PARALLEL, interval

from tracer_lab_gt4py.gt4py_config import backend, backend_opts, dtype_float

logging.getLogger(__name__)

@stencil(backend=backend, **backend_opts)
def sl_init(
    vx_e: Field[dtype_float],
    vy_e: Field[dtype_float],
    vx: Field[dtype_float],
    vy: Field[dtype_float],
    vx_p: Field[dtype_float],
    vy_p: Field[dtype_float],
    lsettls: bool = True,
):
    """Initialize draft velocities with either 
    LSETTLS method : 2 fields for velocity (at t and t - dt)
    LNESN (not LSETTLS) method : 1 field for velocity (at t)

    Args:
        vx_e (Field[dtype]): outlined velocity at t + dt on x
        vy_e (Field[dtype]): outlined velocity at t + dt on y
        vx (Field[dtype]): velocity at t on x
        vy (Field[dtype]): velocity at t on y
        vx_p (Field[dtype]): velocity at t - dt on x
        vy_p (Field[dtype]): velcoity at t - dt on y 
        lsettls (bool, optional): LSETTLS or LNESC. Defaults to True.

    Returns:
        Tuple[Field[dtype]]: velocities at t and t + dt
    """
    
    # LSETTLS
    with computation(PARALLEL), interval(...):
        if lsettls:
            vx_e[0, 0, 0] = vx[0, 0, 0]
            vy_e[0, 0, 0] = vy[0, 0, 0]

            vx[0, 0, 0] = 2 * vx[0, 0, 0] - vx_p[0, 0, 0]
            vy[0, 0, 0] = 2 * vy[0, 0, 0] - vy_p[0, 0, 0]
            
        else:
                
            vx_e[0, 0, 0] = vx[0, 0, 0]
            vy_e[0, 0, 0] = vy[0, 0, 0]
                
