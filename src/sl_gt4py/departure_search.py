from gt4py.cartesian.gtscript import stencil, floor, Field

from sl_gt4py.gt4py_config import backend, backend_opts, dtype, dtype_int

@stencil(backend=backend, **backend_opts)
def dep_search_1d(
    grid_indices: Field[dtype_int],
    vx_e: Field[dtype],
    vx_tmp: Field[dtype],
    lx: Field[dtype],
    departure_indices: Field[dtype_int],
    dx: dtype,
    dth: dtype,
):
    """Compute departure point coordinate (1d)
    
    Args:
        I (Field[dtype]): _description_
        vx_e (Field[dtype]): velocity at arrival point (t + dt)
        vx_tmp (Field[dtype]): estimate of velocity at departure point
        dx (Field[dtype]): grid spacing
        dth (float): half model time step

    Returns:
        Tuple[Field[dtype]]: 
            i_d: indice of departure point on grid 
            lx: adimensionned spacing of departure point from ref grid point 
    """
    
    with computation(PARALLEL), interval(...):
        # Deplacement
        trajx = - dth * (vx_e[0, 0, 0] + vx_tmp[0, 0, 0]) / dx # dth = dt / 2
        departure_indices[0, 0, 0] = (grid_indices + floor(trajx))
        lx[0, 0, 0] = trajx - floor(trajx)