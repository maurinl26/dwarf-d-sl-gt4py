# JAX implementation of the semi-lagrangian advection scheme

from sl_jax.sl_2D import (
    dep_search_1d,
    sl_init,
    backup,
    lagrangian_search,
    sl_xy,
)

from sl_jax.interpolation.interpolation_2d import (
    interpolate_lin_2d,
    interpolate_cub_2d,
    max_interpolator_2d,
    min_interpolator_2d,
)

from sl_jax.boundaries import boundaries

from sl_jax.diagnostics import (
    diagnostic_lipschitz,
    diagnostic_overshoot,
    diagnostic_undershoot,
)

from sl_jax.filter import (
    overshoot_filter,
    undershoot_filter,
)

from sl_jax.periodic_filters import (
    periodic_overshoot_filter,
    periodic_undershoot_filter,
)

__all__ = [
    # Main functions
    "sl_xy",
    "lagrangian_search",
    "dep_search_1d",
    "sl_init",
    "backup",
    # Interpolation
    "interpolate_lin_2d",
    "interpolate_cub_2d",
    "max_interpolator_2d",
    "min_interpolator_2d",
    # Boundaries
    "boundaries",
    # Diagnostics
    "diagnostic_lipschitz",
    "diagnostic_overshoot",
    "diagnostic_undershoot",
    # Filters
    "overshoot_filter",
    "undershoot_filter",
    "periodic_overshoot_filter",
    "periodic_undershoot_filter",
]
