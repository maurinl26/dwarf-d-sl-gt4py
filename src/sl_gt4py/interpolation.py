from gt4py.cartesian.gtscript import Field, stencil
from config import dtype


@stencil("")
def interpolate(
    lx: Field[dtype],
    ly: Field[dtype],
    ix: Field[dtype_int],
    iy: Field[dtype_int],
    tracer: GlobalTable[dtype],  
    cubic_tracer: Field[dtype] 
):
    
    with computation(PARALLEL), interval(...):
        
        first_line = lx * tracer.A[ix, iy] + (1 - lx) * tracer.A[ix + 1, iy]
        sec_line = lx * tracer.A[ix, iy + 1] + (1 - lx) * tracer.A[ix + 1, iy + 1]
        
        cubic_tracer = ly * first_line + (1 - ly) * sec_line
    