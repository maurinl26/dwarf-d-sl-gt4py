import numpy as np
import dace

@dace.program
def boundaries(
    indices: np.ndarray,
    n: dace.int32,
    bc_kind: dace.int32
) -> np.ndarray:
    """Apply boundary conditions
    1: periodic
    0: fixed

    Args:
        indices (int): list of indices to shrink
        n (int): size of spatial domain
        bc_kind (int): type of boundaries 

    Returns:
        np.ndarray: processed indices
    """
    
    if bc_kind == 0:
        id_m1 = np.where(indices < n, indices, n - 1)
        id_m1 = np.where(id_m1 >= 0, id_m1, 0)

    else:
        id_m1 = np.where(indices < n, indices, indices % n)
        id_m1 = np.where(id_m1 >= 0, id_m1, id_m1 % n)

    return id_m1