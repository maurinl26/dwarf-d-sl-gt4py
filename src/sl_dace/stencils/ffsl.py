from ifs_physics_common.framework.stencil import stencil_collection, function_collection
from gt4py.cartesian.gtscript import computation, PARALLEL, interval, Field, IJK, floor, exp, function


# 1. remap velocity
@stencil_collection("velocity_on_faces_x")
def velocity_on_faces_x(
        vx: Field[IJK, float],
        vxh: Field[IJK, float],
):
    """
    compute face velocity given mass point velocity
    with linear interpolation

    faces points are indexed from 0 to nx + 1
    mass points are indexed from 0 to nx

    mass point i is between face point i and face point i + 1

    :param vx: velocity on x (mass point)
    :param vy: velocity on y (mass point)
    :param vxh: velocity on x (face point)
    :param vyh: velocity on y (face point)
    """
    with computation(PARALLEL), interval(...):
        vxh[0, 0, 0] = 0.5 * (vx[-1, 0, 0] + vx[0, 0, 0])

@stencil_collection("velocity_on_faces_y")
def velocity_on_faces_y(
        vy: Field[IJK, float],
        vyh: Field[IJK, float]
):
    """
        compute face velocity given mass point velocity
        with linear interpolation

        faces points are indexed from 0 to nx + 1
        mass points are indexed from 0 to nx

        mass point i is between face point i and face point i + 1

        :param vy: velocity on y (mass point)
        :param vyh: velocity on y (face point)
        """

    with computation(PARALLEL), interval(...):
        vyh[0, 0, 0] = 0.5 * (vy[0, -1, 0] + vy[0, 1, 0])


# 1.5 convert velocity to integer and fractional cfl
@stencil_collection("split_cfl_x")
def split_cfl_x(
        vxh: Field[IJK, float],
        cxh_int: Field[IJK, int],
        cxh_frac: Field[IJK, float],
        dx: float,
        dt: float
):
    """
    Compute the fractional and the integer parts of the cfl on the faces
    of the domain.

    :param vxh: velocity on faces
    :param cxh_int: integer part of the cfl
    :param cxh_frac: fractional part of the cfl
    :param dx: grid spacing on x
    :param dt: timestep
    """
    with computation(PARALLEL), interval(...):
        cxh =  - vxh * dt / dx
        cxh_int = floor(cxh)
        chx_frac = cxh - floor(cxh)


@stencil_collection("split_cfl_y")
def split_cfl_y(
        vyh: Field[IJK, float],
        cyh_int: Field[IJK, int],
        cyh_frac: Field[IJK, float],
        dy: float,
        dt: float
):
    """

    :param vyh:
    :param cyh_int:
    :param cyh_frac:
    :param dy:
    :param dt:
    """
    with computation(PARALLEL), interval(...):
        cyh =  - vyh * dt / dy
        cyh_int = floor(cyh)
        chy_frac = cyh - floor(cyh)


# 2. Interpolate h (advected tracer)
@stencil_collection("fourth_order_facet_interpolation_x")
def fourth_order_facet_interpolation_x(
        psihx: Field[IJK, float],
        psi: Field[IJK, float]
):
    """
    4th order interpolation to compute facet values given mass point values

    :param psih: tracer field at facet points
    :param psi:  tracer field at mass points
    """
    with computation(PARALLEL), interval(...):
        psihx[0, 0, 0] = (psi[-2, 0, 0] + 7 * psi[-1, 0, 0] + 7 * psi[0, 0, 0] + psi[1, 0, 0]) / 12


@stencil_collection("fourth_order_facet_interpolation_y")
def fourth_order_facet_interpolation_y(
        psihy: Field[IJK, float],
        psi: Field[IJK, float]
):
    with computation(PARALLEL), interval(...):
        psihy[0, 0, 0] = (psi[0, -2, 0] + 7 * psi[0, -1, 0] + 7 * psi[0, 0, 0] + psi[0, 1, 0]) / 12


# 4. PPM limiter
@stencil_collection("monotonic_limiter_x")
def monotonic_limiter_x(
    psihx: Field[IJK, float],
    psi: Field[IJK, float]
):
    """
    Monotonic limiter for face values.
    Face value lies between neighbouring mass point values

    :param psih: tracer field at facet point
    :param psi: tracer field at mass point
    """
    with computation(PARALLEL), interval(...):
        psihx[0, 0, 0] = min(
            max(psi[-1, 0, 0], psi[0, 0, 0]),
            max(psihx[0, 0, 0], min(psi[-1, 0, 0], psi[0, 0, 0]))
        )


@function_collection("sigmoid")
@function
def sigmoid(
    x: float
):
    """ Sigmoid activation function """
    return exp(x) / (1 + exp(x))


@stencil_collection("soft_monotonic_limiter_x")
def soft_monotonic_limiter_x(
    psihx: Field[IJK, float],
    psi: Field[IJK, float]
):
    """ Differentiable limiter with min max """
    with computation(PARALLEL), interval(...):
        normalized_value = sigmoid(
            (psihx[0, 0, 0] - psi[0, 0, 0]) /
            (psi[-1, 0, 0] - psi[0, 0, 0])
        )
        psihx[0, 0, 0] = (
            psi[0, 0, 0] + (psi[-1, 0, 0] - psi[0, 0, 0]) * normalized_value
        )

@stencil_collection("monotonic_limiter_y")
def monotonic_limiter_y(
    psihy: Field[IJK, float],
    psi: Field[IJK, float]
):
    with computation(PARALLEL), interval(...):
        psihy[0, 0, 0] = min(
            max(psi[0, -1, 0], psi[0, 0, 0]),
            max(psihy[0, 0, 0], min(psi[0, -1, 0], psi[0, 0, 0]))
        )


@stencil_collection("soft_monotonic_limiter_y")
def soft_monotonic_limiter_y(
    psihx: Field[IJK, float],
    psi: Field[IJK, float]
):
    """ Differentiable limiter with min max """
    with computation(PARALLEL), interval(...):
        normalized_value = sigmoid(
            (psihx[0, 0, 0] - psi[0, 0, 0]) /
            (psi[0, -1, 0] - psi[0, 0, 0])
        )
        psihx[0, 0, 0] = (
            psi[0, 0, 0]
            + (psi[0, -1, 0] - psi[0, 0, 0]) * normalized_value
        )


# note : integer part of the flux is computed in numpy
# sum of integer and fractional part also in numpy

# 4.5 sum of fractional and integer flux
@stencil_collection("integer_and_fractional_flux_sum")
def integer_and_fractional_flux_sum(
        fhx: Field[IJK, float],
        fhx_int: Field[IJK, float],
        fhx_frac: Field[IJK, float]
):
    """
    Sum integer and fractional part of the flux.

    :param fhx: total flux
    :param fhx_int: integer part of the flux
    :param fhx_frac: fractional part of the flux
    """
    with computation(PARALLEL), interval(...):
        fhx = fhx_frac + fhx_int


# 5. density update with updated flux
@stencil_collection("inner_density_update_x")
def inner_density_update_x(
        rho: Field[IJK, float],
        rho_ix: Field[IJK, float],
        fhx: Field[IJK, float],
        ds_yz: float,
        dv: float,
        dt: float,
):
    """
    Compute the inner density update on x axis

    :param rho: density to update
    :param fhx: flux at cell face
    :param ds_yz: section of cell face
    :param dv: volume of cell
    :param dt: time step
    """
    with computation(PARALLEL), interval(...):
        rho_ix = rho - dt * (fhx[0, 0, 0] * ds_yz - fhx[0, 0, 0] * ds_yz) / dv


@stencil_collection("inner_density_update_y")
def inner_density_update_y(
        rho: Field[IJK, float],
        rho_iy: Field[IJK, float],
        fhy: Field[IJK, float],
        ds_xz: float,
        dv: float,
        dt: float
):
    """
    Compute the inner density update on y axis

    :param rho: density
    :param fhy: flux on face on y
    :param ds_xz: cell's normal section to y axis
    :param dv: volume of cell
    :param dt: time step
    """
    with computation(PARALLEL), interval(...):
        rho_iy = rho - dt * (fhy[0, 0, 0] * ds_xz - fhy[0, 0, 0] * ds_xz) / dv


# todo : implement SWIFT outer steps for splitting
def swift_outer_density_update(
        rho_ay: Field[IJK, float],
        sigma_x
):
    ...


