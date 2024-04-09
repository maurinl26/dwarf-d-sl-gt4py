<<<<<<< HEAD
import gt4py.cartesian as gtscript
=======
from gt4py.cartesian.gtscript import Field, stencil, function
>>>>>>> 024babf82bebe5e64c3cff23bac714c34320f09a

from sl_gt4py.gt4py_config import dtype, backend
from fvms.geometry.coordinates import Grid

@function(backend)
def copya2b(
    a: Field[dtype],
    b: Field[dtype]
):
    b[0, 0, 0] = a[0, 0, 0]
    return b
        
@stencil(backend)     
def backup(
    vx: Field[dtype],
    vy: Field[dtype],
    vx_e: Field[dtype],
    vy_e: Field[dtype],
    tracer: Field[dtype],
    tracer_e: Field[dtype]
):
    with computation(PARALLEL), interval(...):
        vx = copya2b(vx_e, vx)
        vy = copya2b(vy_e, vy)
        tracer = copya2b(tracer_e, tracer)
 

    
        
        
