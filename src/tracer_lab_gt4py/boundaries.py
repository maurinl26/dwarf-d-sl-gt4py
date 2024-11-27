from gt4py.cartesian.gtscript import Field, stencil, function
from tracer_lab_gt4py.gt4py_config import backend, externals
import numpy as np

#################################################
######### Stencil version #######################
#################################################
@stencil(backend=backend, externals=externals)
def boundaries(
    idx_dep: Field[np.int64],
    idy_dep: Field[np.int64],
):
    
    from __externals__ import (
        bcx_kind, 
        bcy_kind, 
        nx, 
        ny
    )
    
    with computation(PARALLEL), interval(...):
    
        # PRESCRIBED BOUNDARIES
        if __INLINED(bcx_kind == 0):
            idx_dep = max(0, min(idx_dep, nx))
        
        # PERIODIC BOUNDARIES    
        if __INLINED(bcx_kind == 1):
            idx_dep = idx_dep % nx
        
        
        if __INLINED(bcy_kind == 0):
            idy_dep = max(0, min(idy_dep, ny))
        if __INLINED(bcy_kind == 1):
            idy_dep = idy_dep % ny
        
#########################################
####### Function version ################
#########################################    
@function
def bc_1d(
    id_dep: Field[np.float64],
    n: np.int64,
    bc_kind: np.int64
):
    """Computes the departure index given the boundaries
    
    bc_kind == 1 -> periodic
    bc_kind == 0 -> prescribed

    Args:
        id_dep (Field[np.float64]): index
        n (np.int64): size of the domain along the axis
        bc_kind (np.int64): type of boundary
    """
    
    if bc_kind == 0:
        id_dep = max(0, min(id_dep, n))
    if bc_kind == 1:
        id_dep = id_dep % n
        
    return id_dep

        
        
    
    