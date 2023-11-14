import numpy as np

p_1 = lambda l: (1 / 6) * l * (l - 1) * (2 - l)
p0 = lambda l: (1 - l**2) * (1 - l / 2)
p1 = lambda l: (1 / 2) * l * (l + 1) * (2 - l)
p2 = lambda l: (1 / 6) * l * (l**2 - 1)


def interpolate_lin_2d(
    psi: np.ndarray,
    lx: np.ndarray,
    ly: np.ndarray,
    i_d: np.ndarray,
    j_d: np.ndarray,
    bcx_kind: int,
    bcy_kind: int,
    nx: int,
    ny: int
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
    # Polynômes de lagrange
    p0 = lambda l: 1 - l
    p1 = lambda l: l
    
    # Interp selon x -> field_hat_x
    px = np.array([p0(lx), p1(lx)])
    py = np.array([p0(ly), p1(ly)])
    
    # 1. Construire les tableaux d'indices i_d0, i_d1 / j_d0, j_d1
    # Non periodique
    if bcx_kind == 0: 
        id_0 = np.where(i_d < nx, i_d, nx)
        id_0 = np.where(id_0 >= 0, id_0, 0)
        
        id_p1 = np.where(i_d + 1 < nx, i_d + 1, nx)
        id_p1 = np.where(id_p1 >=0, id_p1, 0)
        
    # Periodique
    else:
        id_0 = np.where(i_d < nx, i_d, i_d % nx)
        id_0 = np.where(i_d >= 0, i_d, i_d % nx)
        id_p1 = np.where(i_d + 1 < nx, i_d + 1, 0)
    
    # Non periodique
    if bcy_kind == 0:
        jd_0 = np.where(j_d < ny, j_d, ny)
        jd_0 = np.where(jd_0 >= 0, jd_0, 0)
        jd_p1 = np.where(j_d + 1 < ny, j_d + 1, ny) 
    # Periodique
    else:
        jd_0 = np.where(j_d < ny, j_d, j_d % ny)
        jd_0 = np.where(jd_0 >= 0, jd_0, jd_0 % ny)
        jd_p1 = np.where(j_d + 1 < ny, j_d + 1, 0)
    
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

# Described as arrays
def interpolate_cub_2d():
    NotImplemented


def interpolate_linear_2d(
    lx: np.float64,
    ly: np.float64,
    ii: np.int64,
    jj: np.int64,
    field: np.ndarray,
    bc_kind: int,
):
    """Interpolate sequentially on x axis and on y axis

    Args:
        l (np.float64): _description_
        ii (int): _description_
        field (np.ndarray): _description_

    Returns:
        _type_: _description_
    """
    p0 = lambda l: 1 - l
    p1 = lambda l: l

    px = np.array([p0(lx), p1(lx)])
    py = np.array([p0(ly), p1(ly)])

    if bc_kind == 0:
        padded_field = np.pad(field, (1, 2), "edge")

    # Periodic
    else:
        padded_field = np.pad(field, (1, 2), "wrap")

    return np.dot(np.matmul(padded_field[ii + 1 : ii + 3, jj + 1 : jj + 3], px), py)


def interpolate_cubic_2d(
    lx: np.float64,
    ly: np.float64,
    ii: np.int64,
    jj: np.int64,
    field: np.ndarray,
    bc_kind: int,
) -> np.float64:
    """Interpolate sequentially on x axis and on y axis.

    Args:
        lx (np.float64): _description_
        ly (np.float64): _description_
        ii (int): _description_
        jj (int): _description_
        field (np.ndarray): _description_

    Returns:
        _type_: _description_
    """

    # Padding on interpolation field
    # Fixed
    if bc_kind == 0:
        padded_field = np.pad(field, (1, 2), "edge")

    # Periodic
    else:
        padded_field = np.pad(field, (1, 2), "wrap")

    # Polynomes de lagrange d'ordre 3
    p_1 = lambda l: (1 / 6) * l * (l - 1) * (2 - l)
    p0 = lambda l: (1 - l**2) * (1 - l / 2)
    p1 = lambda l: (1 / 2) * l * (l + 1) * (2 - l)
    p2 = lambda l: (1 / 6) * l * (l**2 - 1)

    px = np.array([p_1(lx), p0(lx), p1(lx), p2(lx)])
    py = np.array([p_1(ly), p0(ly), p1(ly), p2(ly)])

    assert round(sum(px)) == 1
    assert round(sum(py)) == 1

    return np.dot(np.matmul(padded_field[ii : ii + 4, jj : jj + 4].T, px), py)


def cubic_interpolation(
    lx: np.ndarray,
    ly: np.ndarray,
    ii: np.int64,
    jj: np.int64,
    field: np.int64,
    bc_kind: int,
):
    if bc_kind == 0:
        padded_field = np.pad(field, (1, 2), "edge")

    # Periodic
    else:
        padded_field = np.pad(field, (1, 2), "wrap")

    py_1 = (
        p_1(lx) * padded_field[ii, jj]
        + p0(lx) * padded_field[ii + 1, jj]
        + p1(lx) * padded_field[ii + 2, jj]
        + p2(lx) * padded_field[ii + 3, jj]
    )

    py0 = (
        p_1(lx) * padded_field[ii, jj + 1]
        + p0(lx) * padded_field[ii + 1, jj + 1]
        + p1(lx) * padded_field[ii + 2, jj + 1]
        + p2(lx) * padded_field[ii + 3, jj + 1]
    )

    py1 = (
        p_1(lx) * padded_field[ii, jj + 2]
        + p0(lx) * padded_field[ii + 1, jj + 2]
        + p1(lx) * padded_field[ii + 2, jj + 2]
        + p2(lx) * padded_field[ii + 3, jj + 2]
    )

    py2 = (
        p_1(lx) * padded_field[ii, jj + 3]
        + p0(lx) * padded_field[ii + 1, jj + 3]
        + p1(lx) * padded_field[ii + 2, jj + 3]
        + p2(lx) * padded_field[ii + 3, jj + 3]
    )

    p = py_1 * p_1(ly) + py0 * p0(ly) + py1 * p1(ly) + py2 * p2(ly)

    return p