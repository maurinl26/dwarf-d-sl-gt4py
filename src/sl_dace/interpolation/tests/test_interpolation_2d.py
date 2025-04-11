import pytest
import numpy as np
from sl_dace.interpolation.interpolation_2d import interpolate_lin_2d

def test_interpolate_lin_2d():

    psi = np.ones((50, 50, 10), dtypes=np.float32)
    lx = np.ones((50, 50, 10), dtypes=np.float32)
    ly = np.ones((50, 50, 10), dtypes=np.float32)

    i_d = np.ones((50, 50, 10), dtypes=np.int32)
    j_d = np.ones((50, 50, 10), dtypes=np.int32)

    bcx_kind = 0
    bcy_kind = 0
    nx = 50
    ny = 50

    interpolate_lin_2d(
            psi=psi,
            lx=lx,
            ly=ly,
            i_d=i_d,
            j_d=j_d,
            bcx_kind=bcx_kind,
            bcy_kind=bcy_kind,
            nx=nx,
            ny=ny,
    )