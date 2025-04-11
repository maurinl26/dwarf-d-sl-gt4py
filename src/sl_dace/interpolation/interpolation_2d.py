import numpy as np
import dace
from sl_dace.boundaries import boundaries

@dace.program
def interpolate_lin_2d(
    psi: np.ndarray,
    lx: np.ndarray,
    ly: np.ndarray,
    i_d: np.ndarray,
    j_d: np.ndarray,
    bcx_kind: dace.int32,
    bcy_kind: dace.int32,
    nx: dace.int32,
    ny: dace.int32,
):
    """Perform a 1d linear interpolation

    Args:
        lx (np.ndarray): _description_
        ly (np.ndarray): _description_
        psi (np.ndarray): _description_
        i_d (np.ndarray): _description_
        j_d (np.ndarray): _description_
        bcx_kind (int): _description_
        bcy_kind (int): _description_
        nx (int): _description_
        ny (int): _description_
    """
    # Lagrange lineaire
    # Interp selon x -> field_hat_x
    px = np.array([1 - lx, lx])
    py = np.array([1 - ly, ly])

    # 1. Construire les tableaux d'indices i_d0, i_d1 / j_d0, j_d1
    # Non periodique
    id_0 = boundaries(i_d, nx, bcx_kind)
    id_p1 = boundaries(i_d + 1, nx, bcx_kind)
    
    jd_0 = boundaries(j_d, ny, bcy_kind)
    jd_p1 = boundaries(j_d + 1, ny, bcy_kind)

    # Lookup
    psi_d_i = np.zeros((4, nx, ny))
    for i in range(nx):
        for j in range(ny):
            psi_d_i[0, i, j] = psi[id_0[i, j], jd_0[i, j]]
            psi_d_i[1, i, j] = psi[id_p1[i, j], jd_0[i, j]]

            psi_d_i[2, i, j] = psi[id_0[i, j], jd_p1[i, j]]
            psi_d_i[3, i, j] = psi[id_p1[i, j], jd_p1[i, j]]

    psi_d_j = np.zeros((2, nx, ny))
    psi_d_j[0] = px[0] * psi_d_i[0] + px[1] * psi_d_i[1]
    psi_d_j[1] = px[0] * psi_d_i[2] + px[1] * psi_d_i[3]

    psi_d = py[0] * psi_d_j[0] + py[1] * psi_d_j[1]

    return psi_d

