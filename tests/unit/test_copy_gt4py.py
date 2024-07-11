import numpy as np
import gt4py

from sl_gt4py.copy import copy
from sl_gt4py.gt4py_config import dtype, backend, origin
from sl_python.blossey import init_blossey



if __name__ == "__main__":
        
    nx, ny = 50, 50
    nz = 1
    dt, dx, dy = 0.02, 0.02, 0.02
    xcr, ycr = np.meshgrid(np.linspace(0,nx), np.linspace(0,ny))
    
    # Velocities
    vx = gt4py.storage.zeros((nx, ny, nz), dtype, backend=backend, aligned_index=origin)
    vy = gt4py.storage.zeros((nx, ny, nz), dtype, backend=backend, aligned_index=origin)
    vx_tmp = gt4py.storage.zeros((nx, ny, nz), dtype, backend=backend, aligned_index=origin)
    vy_tmp = gt4py.storage.zeros((nx, ny, nz), dtype, backend=backend, aligned_index=origin)
   
    # Init velocities (numpy -> gt4py)
    vx_0, vy_0, vx_p_0, vy_p_0, vx_e_0, vy_e_0 = init_blossey(
         xcr, ycr, 0, dt, dx, dy, nx, ny
    )
    vx[:, :, :] = vx_0[:, :, np.newaxis]
    vy[:, :, :] = vy_0[:, :, np.newaxis]

    copy(
        vx, vx_tmp, vy, vy_tmp
    )
    
    assert(vx.any() == vx_tmp.any())
    assert(vy.any() == vy_tmp.any())
