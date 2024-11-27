import sys
import numpy as np

sys.path.append("/home/maurinl/sl_gt4py/src")

from tracer_lab_python.interpolation import interpolate_lin_2d
from utils.config import Config

def r_tilde(x: np.ndarray, y: np.ndarray):
    return 5 * np.sqrt((x - 0.5) ** 2 + (y - 0.5) ** 2)


def tracer_shape(x: np.ndarray, y: np.ndarray, tracer0: np.ndarray):
    r_t = r_tilde(x, y)
    return np.where(r_t < 1, tracer0 + ((1 + np.cos(np.pi * r_t)) / 2) ** 2, tracer0)


if __name__ == "__main__":
    
    model_starttime = 0
    model_endtime = 1
    nstep = 50
    nitmp = 4
    dt = (model_endtime - model_starttime) / nstep
    xmin, xmax = 0, 1
    ymin, ymax = 0, 1
    nx, ny = 50, 50
    U, V = 1, 1
    lsettls = True
    bcx_kind, bcy_kind = 1, 1

    config = Config(dt, xmin, xmax, nx, ymin, ymax, ny, bcx_kind, bcy_kind)
    
    psi = tracer_shape(config.xcr, config.ycr, 1)

    field = interpolate_lin_2d(
        psi,
        np.zeros((nx, ny)),
        np.zeros((nx, ny)),
        config.I,
        config.J,
        1,
        1,
        nx,
        ny
    )
    
    print(f"Interp Lin : min {np.min(field)}, max {np.max(field)}, mean {np.mean(field)}")
    print(f"Errors : L1 {np.max(np.abs(field - psi)):.02f}, L2 {(1 / (nx * ny)) * np.sum((field - psi) ** 2):.02f}")
    
    field = interpolate_lin_2d(
        psi,
        0.5 * np.ones((nx, ny)),
        0.5 * np.ones((nx, ny)),
        config.I,
        config.J,
        1,
        1,
        nx,
        ny
    )
    
    print(f"Interp Lin : min {np.min(field)}, max {np.max(field)}, mean {np.mean(field)}")
    print(f"Errors : L1 {np.max(np.abs(field - psi)):.02f}, L2 {(1 / (nx * ny)) * np.sum((field - psi) ** 2):.02f}")
    
    field = interpolate_lin_2d(
        psi,
        np.zeros((nx, ny)),
        np.zeros((nx, ny)),
        config.I + 1,
        config.J,
        1,
        1,
        nx,
        ny
    )
    
    print(f"Interp Lin : min {np.min(field)}, max {np.max(field)}, mean {np.mean(field)}")
    print(f"Errors : L1 {np.max(np.abs(field - psi)):.02f}, L2 {(1 / (nx * ny)) * np.sum((field - psi) ** 2):.02f}")