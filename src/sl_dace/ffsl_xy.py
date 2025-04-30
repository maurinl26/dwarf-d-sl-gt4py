from collections.abc import dict_values
from ifs_physics_common.framework.storage import managed_temporary_storage
from gt4py.cartesian.gtscript import stencil
from functools import partial, cached_property
import numpy as np
from config import Config
from sl_dace.stencils.ffsl import (
    velocity_on_faces_x,
    velocity_on_faces_y,
    split_cfl_x,
    split_cfl_y,
    fourth_order_facet_interpolation_x,
    fourth_order_facet_interpolation_y,
    monotonic_limiter_x,
    monotonic_limiter_y,
    fractional_flux_x,
    inner_density_update_x,
    fractional_and_integer_flux_sum
)
from sl_dace.ppm import ppm_coefficients_x
from sl_dace.interpolation.flux_integral import (
    integer_flux_integral_x,
    fractional_flux_integral_x
)
from typing import Union


class FluxFormSemiLagXY:

    def __init__(self,
                 config: Config,
                 splitting_mode: Union["COSMIC", "SWIFT"] = "SWIFT"
                 ):
        self.splitting_mode = splitting_mode


    def __call__(self):

        match self.splitting_mode:
            case "SWIFT":
                ...
                #FluxFormSemiLagX(sigma)
                #FluxFormSemiLagY(sigma)

                # FluxFormSemiLagX(rho)
                # FluxFormSemiLagY(rho)

                # FluxFormSemiLagY(rhox / sigmax)
                # FluxFormSemiLagX(rhoy / sigmay)
                # todo: partial operators for pure advection


            case "COSMIC":
                ...





