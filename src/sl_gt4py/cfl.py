from gt4py.cartesian import gtscript
from sl_gt4py.gt4py_config import dtype

@gtscript.stencil()
def cfl_1d(u: gtscript.Field[dtype], cfl_u: gtscript.Field[dtype], dx: dtype):
    
    with computation(PARALLEL), interval(...):
        cfl_u = u / dx