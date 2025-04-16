import pytest
import numpy as np
from sl_dace.interpolation.interpolation_2d import interpolate_lin_2d
import dace

def test_interpolate_lin_2d():

    psi = np.ones((50, 50, 10), dtype=np.float32)
    lx = np.zeros((50, 50, 10), dtype=np.float32)
    ly = np.zeros((50, 50, 10), dtype=np.float32)

    i_d = np.ones((50, 50, 10), dtype=np.int32)
    j_d = np.ones((50, 50, 10), dtype=np.int32)

    psi_dep = np.zeros((50, 50, 10), dtype=np.float32)

    nx, ny, nz = (50, 50, 15)
    h = 5

    c_interpolate_lin_2d = (
        dace.program(interpolate_lin_2d)
        .to_sdfg()
        .compile()
    )

    c_interpolate_lin_2d(
            psi=psi,
            lx=lx,
            ly=ly,
            i_dep=i_d,
            j_dep=j_d,
            psi_dep=psi_dep,
            I=nx,
            J=ny,
            K=nz,
            H=h
    )

    assert psi_dep.mean() == 0