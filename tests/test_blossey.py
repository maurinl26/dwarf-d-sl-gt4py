import logging
import time
from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("/home/maurinl/FVM_GT4Py_slim/src")
sys.path.append("/home/maurinl/sl_gt4py/src")

from sl_python.interpolation import interpolate_cubic_2d
from sl_python.sl_2D import sl_init, sl_xy
from config import Config
from sl_python.plot import plot_2D_scalar_field

logging.basicConfig(stream=sys.stdout, level=logging.INFO)

# Pointwise functions
def r(x, y):
    return np.sqrt((x - 0.5) ** 2 + (y - 0.5) ** 2)


def r_tilde(x, y):
    return 5 * np.sqrt((x - 0.3) ** 2 + (y - 0.5) ** 2)


def theta(x, y):
    return np.arctan((y - 0.5) / (x - 0.5))


def u(r, t, T):
    return (4 * np.pi * r / T) * (
        1 + np.cos(2 * np.pi * t / T) * (1 - (4 * r) ** 6) / (1 + (4 * r) ** 6)
    )


def cfl_1d(u: np.float64, dx: np.float64, dt: np.float64):
    return u * dt / dx


def tracer_shape(x: np.float64, y: np.float64, tracer0: np.float64):
    r_t = r_tilde(x, y)
    return np.where(r_t < 1, tracer0 + ((1 + np.cos(np.pi * r_t)) / 2) ** 2, tracer0)


# Field functions
def polar_coordinates(xcr: np.ndarray, ycr: np.ndarray):
    R = r(xcr, ycr)
    Theta = theta(xcr, ycr)
    return R, Theta

# Taken from Christian Kühnlein
def blossey_stf(
    xcr: np.ndarray,
    ycr: np.ndarray,
    mt: np.float64
)-> np.ndarray:
    """Compute blossey stream function as described by 
    Durran & Blossey (Selective monotonicity preservation in scalar advection)

    Args:
        t (np.float64): _description_
        T (np.float64): _description_
        r (np.ndarray): _description_

    Returns:
        _type_: _description_
    """
    rsq = (xcr - 0.5) ** 2.0 + (ycr - 0.5) ** 2.0

    term1 = 0.5 * rsq
    term2 = np.log(1.0 - 16.0 * rsq + 256.0 * rsq**2.0) / 96.0
    term3 = -np.log(1.0 + 16.0 * rsq) / 48.0
    term4 = -np.sqrt(3.0) * np.arctan((-1.0 + 32.0 * rsq) / np.sqrt(3.0)) / 48.0

    stf = 4.0 * np.pi * (term1 + np.cos(2.0 * np.pi * mt) * (term1 + term2 + term3 + term4))

    return stf

# Taken from Christian
def stream_function_xy(
    stf: np.ndarray,
    dx: np.float,
    dy: np.float
):
    
    dstf_dx, dstf_dy = np.gradient(stf)
    u = -dstf_dy
    v = dstf_dx
    
    return u, v

def blossey_velocity(
    xcr: np.ndarray, ycr: np.ndarray, mt, dx, dy
) -> Tuple[np.ndarray]:
    rcr, thetacr = polar_coordinates(xcr, ycr)

    stf = blossey_stf(xcr, ycr, mt)
    u, v = stream_function_xy(str, dx, dy)

    return u, v


def init_blossey(xcr: np.ndarray, ycr: np.ndarray, mt, dx, dy, nx: int, ny: int):
    
    # Vitesses
    vx, vy = blossey_velocity(xcr, ycr, mt, dx, dy)
    vx_p, vy_p = vx.copy(), vy.copy()
    vx_e, vy_e = np.empty((nx, ny)), np.empty((nx, ny))

    return vx, vy, vx_p, vy_p, vx_e, vy_e


def blossey_tracer(xcr, ycr):
    # Tracer
    tracer = tracer_shape(xcr, ycr, 0)
    tracer_e = tracer.copy()
    return tracer, tracer_e


def backup(tracer, tracer_e):
    tracer = tracer_e.copy()
    return tracer

# Driver
def sl_driver(
    config: Config, 
    vx: np.ndarray, 
    vy: np.ndarray, 
    vx_e: np.ndarray, 
    vy_e: np.ndarray, 
    vx_p: np.ndarray, 
    vy_p: np.ndarray, 
    tracer: np.ndarray, 
    tracer_e: np.ndarray, 
    lsettls: bool, 
    nitmp: int
):
    tracer_ref = tracer.copy()

    # Advection
    vx_e, vy_e, vx, vy = sl_init(
        vx_e=vx_e, vy_e=vy_e, vx=vx, vy=vy, vx_p=vx_p, vy_p=vy_p, lsettls=lsettls
    )

    for jstep in range(nitmp):
    
        t = jstep * config.dt
        logging.info(f"Step : {jstep}")
        logging.info(f"Step : {jstep}")

        tracer = backup(tracer=tracer, tracer_e=tracer_e)

        vx, vy = blossey_velocity(vx, vy, mt, dx, dy)
        
        plot_blossey(config.xcr, config.ycr, vx, vy, tracer)
        

        # Estimations
        vx_e, vy_e, tracer_e = sl_xy(
            config=config,
            vx=vx,
            vy=vy,
            vx_e=vx_e,
            vy_e=vy_e,
            tracer=tracer,
            tracer_e=tracer_e,
            interpolation_function=interpolate_cubic_2d,
            nsiter=4
        )


        # Diagnostics and outputs
        courant_xmax = np.max(cfl_1d(vx, config.dx, config.dt))
        courant_ymax = np.max(cfl_1d(vy, config.dy, config.dt))
    
        np.save(f"tracer_{jstep}", tracer)
        plot_blossey(config.xcr, config.ycr, vx, vy, tracer)

        logging.info(f"Maximum courant number : {max(courant_xmax, courant_ymax)}")

    e_inf = np.max(np.abs(tracer - tracer_ref))
    e_2 = np.sqrt((1 / (config.nx * config.ny)) * np.sum((tracer - tracer_ref) ** 2))

    logging.info(f"Error E_inf : {e_inf}")
    logging.info(f"Error RMSE : {e_2}")
    
    
    plot_blossey(config.xcr, config.ycr, vx, vy, tracer)
    
    
def plot_blossey(
    xcr: np.ndarray, ycr: np.ndarray, vx: np.ndarray, vy: np.ndarray, tracer: np.ndarray
    ):
    
    # Shape
    fig, ax = plt.subplots()
    
    # Vent
    ax.quiver(xcr, ycr, vx, vy, color="C0", angles="xy", scale_units="xy", scale=1, width=0.001)
    ax.set(xlim=(0, 1), ylim=(0, 1))
    
    levels = [0.05 + i * 0.1 for i in range(0, 10)]
    ax.contour(xcr, ycr, tracer, colors="black", vmin=0.05, vmax=0.95, levels=levels)
    
    plt.show()

if __name__ == "__main__":
    
    # Shift in config file
    T = 1
    t = 0

    model_starttime = 0
    model_endtime = 1
    nitmp = 200
    dt = (model_endtime - model_starttime) / nitmp
    xmin, xmax = 0, 1
    ymin, ymax = 0, 1
    nx = 50
    ny = 50
    lsettls = False
    bcx_kind, bcy_kind = 1, 1

    config = Config(1, xmin, xmax, nx, ymin, ymax, ny, bcx_kind, bcy_kind)

    vx, vy, vx_p, vy_p, vx_e, vy_e = init_blossey(config.xcr, config.ycr, t, T)
    tracer, tracer_e = blossey_tracer(config.xcr, config.ycr)
    plot_blossey(config.xcr, config.ycr, tracer, 50)


    # Advection encapsulation
    start_time = time.time()
    sl_driver(config, vx, vy, vx_e, vy_e, vx_p, vy_p, tracer, tracer_e, lsettls, nitmp)
    duration = time.time() - start_time
    logging.info(f"Duration : {duration} s")