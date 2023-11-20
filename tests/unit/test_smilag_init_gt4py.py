import numpy as np
import sys
import yaml
import gt4py
from typing import Tuple

sys.path.append("/home/maurinl/sl_gt4py/src")
sys.path.append("/home/maurinl/sl_gt5py/tests")
print(sys.path)

from sl_gt4py.gt4py_config import dtype, dtype_int, backend, origin
from sl_gt4py.smilag_init import sl_init
from config import Config
from sl_python.blossey import init_blossey


if __name__ == "__main__":
        
    nx, ny = 50, 50
    nz = 1
    dt, dx, dy = 0.02, 0.02, 0.02
    xcr, ycr = np.meshgrid(np.linspace(0,nx), np.linspace(0,ny))
    
    # Velocities
    vx = gt4py.storage.zeros((nx, ny, nz), dtype, backend=backend, aligned_index=origin)
    vy = gt4py.storage.zeros((nx, ny, nz), dtype, backend=backend, aligned_index=origin)
    vx_p = gt4py.storage.zeros((nx, ny, nz), dtype, backend=backend, aligned_index=origin)
    vy_p = gt4py.storage.zeros((nx, ny, nz), dtype, backend=backend, aligned_index=origin)
    vx_e = gt4py.storage.zeros((nx, ny, nz), dtype, backend=backend, aligned_index=origin)
    vy_e = gt4py.storage.zeros((nx, ny, nz), dtype, backend=backend, aligned_index=origin)
    
    # Init velocities (numpy -> gt4py)
    vx_0, vy_0, vx_p_0, vy_p_0, vx_e_0, vy_e_0 = init_blossey(
         xcr, ycr, 0, dt, dx, dy, nx, ny
    )
    vx[:, :, :] = vx_0[:, :, np.newaxis]
    vy[:, :, :] = vy_0[:, :, np.newaxis]
    vx_p[:, :, :] = vx_p_0[:, :, np.newaxis]
    vy_p[:, :, :] = vy_p_0[:, :, np.newaxis]
    vx_e[:, :, :] = vx_e_0[:, :, np.newaxis]
    vy_e[:, :, :] = vy_e_0[:, :, np.newaxis]   
    
    # LSETTLS  
    lsettls = True
    sl_init(
        vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, lsettls=lsettls
    )
    
    # LNESC
    lsettls = False
    sl_init(
        vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, lsettls=lsettls
    )

