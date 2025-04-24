import numpy as np
from typing import Tuple
from gt4py.cartesian.gtscript import PARALLEL, interval, computation, Field, sqrt, cos, pi, arctan, stencil


def radius(x: Field[IJK, float], y: Field[IJK, float], radius: Field[IJK, float]):

    with computation(PARALLEL), interval(...):
        radius = sqrt((x - 0.5) ** 2 + (y - 0.5) ** 2)

def radius_tilde(x: Field[IJK, float], y: Field[IJK, float], radius_tilde: Field[IJK, float]):
    with computation(PARALLEL), interval(...):
        radius_tilde = 5 * sqrt((x - 0.3) ** 2 + (y - 0.5) ** 2)


def theta(x: Field[IJK, float], y: Field[IJK, float], theta: Field[IJK, float]):
    with computation(PARALLEL), interval(...):
        theta = arctan((y - 0.5) / (x - 0.5))

def u_velocity(radius: Field[IJK, float], t: Field[IJK, float], u: Field[IJK, float], T: float):
    with computation(PARALLEL), interval(...):
        u = (4 * pi * radius / T) * (
                1 + cos(2 * pi * t / T) * (1 - (4 * radius) ** 6) / (1 + (4 * radius) ** 6)
        )

def tracer_shape(radius_tilde: Field[IJK, float], tracer: Field[IJK, float], tracer0: float):
    with computation(PARALLEL), interval(...):
        tracer = tracer0
        tracer += (
            (1 + cos(pi * radius_tilde) / 2) ** 2
            if radius_tilde < 1.0
            else 0
        )

