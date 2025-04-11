from gt4py.cartesian.gtscript import Field, IJK, function

@function
def cfl_1d(u: Field[IJK, float], dx: float, dt: float):
    return u * dt / dx