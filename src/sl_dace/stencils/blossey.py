import dace
from sl_dace.utils.dims import I, J, K
from sl_dace.utils.typingx import dtype_float, dtype_int

@dace.program
def radius(x: dtype_float[I, J, K], y: dtype_float[I, J, K], radius: dtype_float[I, J, K]):

    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        radius[i,j,k] = sqrt((x[i, j, k] - 0.5) ** 2 + (y[i, j, k] - 0.5) ** 2)

@dace.program
def radius_tilde(x: dtype_float[I, J, K], y: dtype_float[I, J, K], radius_tilde: dtype_float[I, J, K]):
    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        radius_tilde[i, j, k] = 5 * sqrt((x[i, j, k] - 0.3) ** 2 + (y[i, j, k] - 0.5) ** 2)

@dace.program
def theta(x: dtype_float[I, J, K], y: dtype_float[I, J, K], theta: dtype_float[I, J, K]):
    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        theta[i, j, k]  = arctan((y[i, j, k]  - 0.5) / (x[i, j, k]  - 0.5))

@dace.program
def u_velocity(radius: dtype_float[I, J, K], t: dtype_float[I, J, K], u: dtype_float[I, J, K], T: float):
    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        u[i, j, k]  = (4 * pi * radius[i, j, k]  / T) * (
                1 + cos(2 * pi * t[i, j, k]  / T) * (1 - (4 * radius[i, j, k] ) ** 6) / (1 + (4 * radius[i, j, k] ) ** 6)
        )

@dace.program
def tracer_shape(radius_tilde: dtype_float[I, J, K], tracer: dtype_float[I, J, K], tracer0: float):
    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        tracer[i, j, k]  = tracer0[i, j, k]
        tracer += (
            (1 + cos(pi * radius_tilde[i, j, k] ) / 2) ** 2
            if radius_tilde[i, j, k]  < 1.0
            else 0
        )

