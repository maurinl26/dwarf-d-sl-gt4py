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
        id_0 = np.where(id_0 >= 0, id_0, id_0 % nx)
        
        id_p1 = np.where(i_d + 1 < nx, i_d + 1, 0)
        id_p1 = np.where(id_p1 >= 0, id_p1, id_p1 % nx)
    
    # Non periodique
    if bcy_kind == 0:
        jd_0 = np.where(j_d < ny, j_d, ny)
        jd_0 = np.where(jd_0 >= 0, jd_0, 0)
        
        jd_p1 = np.where(j_d + 1 < ny, j_d + 1, ny) 
        jd_p1 = np.where(jd_p1 >= 0, jd_p1, 0)
        
    # Periodique
    else:
        jd_0 = np.where(j_d < ny, j_d, j_d % ny)
        jd_0 = np.where(jd_0 >= 0, jd_0, jd_0 % ny)
        
        jd_p1 = np.where(j_d + 1 < ny, j_d + 1, 0)
        jd_p1 = np.where(jd_p1 >= 0, jd_p1, jd_p1 % ny)
        
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
def interpolate_cub_2d(
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
    """_summary_
    """
    # Polynomes de Lagrange
    p_1 = lambda l: (1 / 6) * l * (l - 1) * (2 - l)
    p0 = lambda l: (1 - l**2) * (1 - l / 2)
    p1 = lambda l: (1 / 2) * l * (l + 1) * (2 - l)
    p2 = lambda l: (1 / 6) * l * (l**2 - 1)

    px = np.array([p_1(lx), p0(lx), p1(lx), p2(lx)])
    py = np.array([p_1(ly), p0(ly), p1(ly), p2(ly)])
    
    # 1. Tableaux d'indices
    # Non periodique
    if bcx_kind == 0:
        id_m1 = np.where(i_d - 1 < nx, i_d - 1, nx)
        id_m1 = np.where(id_m1 >= 0, id_m1, 0)
        
        id_0 = np.where(i_d < nx, i_d, nx)
        id_0 = np.where(id_0 >=0, id_0, 0)
        
        id_p1 = np.where(i_d + 1 < nx, i_d + 1, nx)
        id_p1 = np.where(id_p1 >= 0, id_p1, 0)
        
        id_p2 = np.where(i_d + 2 < nx, i_d + 2, nx)
        id_p2 = np.where(id_p2 >= 0, id_p2, 0)
        
    else:
        id_m1 = np.where(i_d - 1 < nx, i_d - 1, (i_d - 1) % nx)
        id_m1 = np.where(id_m1 >= 0, id_m1, id_m1 % nx)
        
        id_0 = np.where(i_d < nx, i_d, i_d % nx)
        id_0 = np.where(id_0 >=0, id_0, id_0 % nx)
        
        id_p1 = np.where(i_d + 1 < nx, i_d + 1, (i_d + 1) % nx)
        id_p1 = np.where(id_p1 >= 0, id_p1, id_p1 % nx)
        
        id_p2 = np.where(i_d + 2 < nx, i_d + 2, (i_d + 2) % nx)
        id_p2 = np.where(id_p2 >= 0, id_p2, id_p2 % nx)
        
    if bcy_kind == 0:
        jd_m1 = np.where(j_d - 1 < nx, j_d - 1, nx)
        jd_m1 = np.where(jd_m1 >= 0, jd_m1, 0)
        
        jd_0 = np.where(j_d < nx, j_d, nx)
        jd_0 = np.where(jd_0 >=0, jd_0, 0)
        
        jd_p1 = np.where(j_d + 1 < nx, j_d + 1, nx)
        jd_p1 = np.where(jd_p1 >= 0, jd_p1, 0)
        
        jd_p2 = np.where(j_d + 2 < nx, j_d + 2, nx)
        jd_p2 = np.where(jd_p2 >= 0, jd_p2, 0)
        
    else:
        jd_m1 = np.where(j_d - 1 < nx, j_d - 1, (j_d - 1) % nx)
        jd_m1 = np.where(jd_m1 >= 0, jd_m1, jd_m1 % nx)
        
        jd_0 = np.where(j_d < nx, j_d, j_d % nx)
        jd_0 = np.where(jd_0 >=0, jd_0, jd_0 % nx)
        
        jd_p1 = np.where(j_d + 1 < nx, j_d + 1, (j_d + 1) % nx)
        jd_p1 = np.where(jd_p1 >= 0, jd_p1, jd_p1 % nx)
        
        jd_p2 = np.where(j_d + 2 < nx, j_d + 2, (j_d + 2) % nx)
        jd_p2 = np.where(jd_p2 >= 0, jd_p2, jd_p2 % nx)
        
    # Lookup 
    psi_d_i = np.zeros((16, nx, ny))
    for i in range(nx):
        for j in range(ny):

            psi_d_i[0, i, j] = psi[id_m1[i, j], jd_m1[i, j]]
            psi_d_i[1, i, j] = psi[id_0[i, j], jd_m1[i, j]]
            psi_d_i[2, i, j] = psi[id_p1[i, j], jd_m1[i, j]]
            psi_d_i[3, i, j] = psi[id_p2[i, j], jd_m1[i, j]]
            
            psi_d_i[4, i, j] = psi[id_m1[i, j], jd_0[i, j]]
            psi_d_i[5, i, j] = psi[id_0[i, j], jd_0[i, j]]
            psi_d_i[6, i, j] = psi[id_p1[i, j], jd_0[i, j]]
            psi_d_i[7, i, j] = psi[id_p2[i, j], jd_0[i, j]]
            
            psi_d_i[8, i, j] = psi[id_m1[i, j], jd_p1[i, j]]
            psi_d_i[9, i, j] = psi[id_0[i, j], jd_p1[i, j]]
            psi_d_i[10, i, j] = psi[id_p1[i, j], jd_p1[i, j]]
            psi_d_i[11, i, j] = psi[id_p2[i, j], jd_p1[i, j]]
            
            psi_d_i[12, i, j] = psi[id_m1[i, j], jd_p2[i, j]]
            psi_d_i[13, i, j] = psi[id_0[i, j], jd_p2[i, j]]
            psi_d_i[14, i, j] = psi[id_p1[i, j], jd_p2[i, j]]
            psi_d_i[15, i, j] = psi[id_p2[i, j], jd_p2[i, j]]
            
    psi_d_j = np.zeros((4, nx, ny))
    psi_d_j[0] = px[0] * psi_d_i[0] + px[1] * psi_d_i[1] + px[2] * psi_d_i[2] + px[3] * psi_d_i[3]    
    psi_d_j[1] = px[0] * psi_d_i[4] + px[1] * psi_d_i[5] + px[2] * psi_d_i[6] + px[3] * psi_d_i[7]
    psi_d_j[2] = px[0] * psi_d_i[8] + px[1] * psi_d_i[9] + px[2] * psi_d_i[10] + px[3] * psi_d_i[11]
    psi_d_j[3] = px[0] * psi_d_i[12] + px[1] * psi_d_i[13] + px[2] * psi_d_i[14] + px[3] * psi_d_i[15]
    
    psi_d = py[0] * psi_d_j[0] + py[1] * psi_d_j[1] + py[2] * psi_d_j[2] + py[3] * psi_d_j[3]
    
    return psi_d

