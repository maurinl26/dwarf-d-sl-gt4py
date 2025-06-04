import dace
from sl_dace.utils.dims import I, J, K

# todo: prescribe boundaries
# 3. PPM coefficients
def ppm_coefficients_x(
        psi: dtype_float[I, J, K],
        psih: dtype_float[I, J, K + 1],
        a0: dtype_float[I, J, K],
        a1: dtype_float[I, J, K],
        a2: dtype_float[I, J, K],
):
    """
    Compute parabolic reconstruction on cell
    Limiter in case of local extrema.

    :param psi: tracer at mass point
    :param psih: tracer at face point
    :param a0: parabolic coefficient (cst)
    :param a1: parabolic coefficient (x)
    :param a2: parabolic coefficient (x**2)
    """
    for i, j, k in dace.map[0:I, 0:J, 0:K]:
        a0[i, j, k] = psih[i-1, j, k]
        a1[i, j, k] = -4 * psih[i, j, k] - 2 * psih[i+1, j, k] + 6 * psi[i, j, k]
        a2[i, j, k] = 3 * psih[i, j, k] + 3 * psih[i+1, j, k] - 6 * psi[i, j, k]

        # Constant reconstruction if local extrema
        if -a1 / (2 * a2) > 0:
            a0[i, j, k] = psi[i, j, k]
            a1[i, j, k] = 0
            a2[i, j, k] = 0
