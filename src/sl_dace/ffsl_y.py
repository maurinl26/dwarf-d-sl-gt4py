from collections.abc import dict_values
from ifs_physics_common.framework.storage import managed_temporary_storage
from gt4py.cartesian.gtscript import stencil
from functools import partial, cached_property
import numpy as np
from config import Config
from sl_dace.stencils.ffsl import (
    velocity_on_faces_y,
    split_cfl_y,
    fourth_order_facet_interpolation_y,
    monotonic_limiter_y,
    fractional_flux_y,
    inner_density_update_y,
    fractional_and_integer_flux_sum
)
from sl_dace.ppm import ppm_coefficients_y
from sl_dace.interpolation.flux_integral import (
    integer_flux_integral_y,
    fractional_flux_integral_y
)
from typing import Union



class FluxFormSemiLagY:

    def __init__(self, config: Config):
        ...

    def __call__(self):
        ...
