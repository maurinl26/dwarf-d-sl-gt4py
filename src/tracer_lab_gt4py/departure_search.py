from gt4py.cartesian.gtscript import stencil, floor, PARALLEL, interval, computation
from gt4py.cartesian import gtscript

from sl_gt4py.gt4py_config import backend, backend_opts, dtype, dtype_int

@stencil(backend=backend, **backend_opts)
def dep_search_1d(
    grid_indices: gtscript.Field[dtype_int],
    vx_e: gtscript.Field[dtype],
    vx_tmp: gtscript.Field[dtype],
    lx: gtscript.Field[dtype],
    departure_indices: gtscript.Field[dtype_int],
    dx: dtype,
    dth: dtype,
):
    """Compute departure point coordinate (1d)
    
    Args:
        I (gtscript.Field[dtype]): _description_
        vx_e (gtscript.Field[dtype]): velocity at arrival point (t + dt)
        vx_tmp (gtscript.Field[dtype]): estimate of velocity at departure point
        dx (gtscript.Field[dtype]): grid spacing
        dth (float): half model time step

    Returns:
        Tuple[gtscript.Field[dtype]]: 
            i_d: indice of departure point on grid 
            lx: adimensionned spacing of departure point from ref grid point 
    """
    
    with computation(PARALLEL), interval(...):
        # Deplacement
        trajx = - dth * (vx_e[0, 0, 0] + vx_tmp[0, 0, 0]) / dx # dth = dt / 2
        departure_indices[0, 0, 0] = (grid_indices + floor(trajx))
        lx[0, 0, 0] = trajx - floor(trajx)