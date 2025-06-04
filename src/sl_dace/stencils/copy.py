import dace
from sl_dace.utils.typingx import dtype_int, dtype_float
from sl_dace.utils.dims import I, J, K

@stencil_collection("copy")
def copy(
    vx_tmp: dtype_float[I, J, K],
    vy_tmp: dtype_float[I, J, K],
    vx: dtype_float[I, J, K],
    vy: dtype_float[I, J, K]
):

    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        vx_tmp = vx
        vy_tmp = vy

def backup(
    tracer: dtype_float[I, J, K],
    tracer_e: dtype_float[I, J, K]
):
    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        tracer[i, j, k] = tracer_e[i, j, k]