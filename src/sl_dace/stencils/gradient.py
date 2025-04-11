from gt4py.cartesian.gtscript import computation, interval, PARALLEL

def c2_x(f: Field[float], df_dx: Field[float], dx: float):

    with computation(PARALLEL), interval(...):
        df_dx = (f[-1, 0, 0] - f[-1, 0, 0]) / (2 * dx)


def c2_y(f: Field[float], df_dy: Field[float], dy: float):

    with computation(PARALLEL), interval(...):
        df_dy = (f[0, -1, 0] - f[0, 1, 0]) / (2 * dy)
