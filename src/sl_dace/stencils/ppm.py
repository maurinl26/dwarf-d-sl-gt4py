from gt4py.cartesian.gtscript import Field, IJK, computation, PARALLEL, interval

# todo: prescribe boundaries
# 3. PPM coefficients
def ppm_coefficients_x(
        psi: Field[IJK, float],
        psih: Field[IJK, float],
        a0: Field[IJK, float],
        a1: Field[IJK, float],
        a2: Field[IJK, float]
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
    with computation(PARALLEL), interval(...):
        a0[0, 0, 0] = psih[-1, 0, 0]
        a1[0, 0, 0] = -4 * psih[0, 0, 0] - 2 * psih[1, 0, 0] + 6 * psi[0, 0, 0]
        a2[0, 0, 0] = 3 * psih[0, 0, 0] + 3 * psih[1, 0, 0] - 6 * psi[0, 0, 0]

        # Constant reconstruction if local extrema
        if -a1 / (2 * a2) > 0:
            a0[0, 0, 0] = psi
            a1[0, 0, 0] = 0
            a2[0, 0, 0] = 0
