from collections.abc import dict_values
from config import Config
from sl_dace.ffsl_x import flux_from_semi_lag_x
import dace


def flux_form_semi_lag_xy():
    ...

    # FluxFormSemiLagX(sigma)
    # FluxFormSemiLagY(sigma)

    # FluxFormSemiLagX(rho)
    # FluxFormSemiLagY(rho)

    # FluxFormSemiLagY(rhox / sigmax)
    # FluxFormSemiLagX(rhoy / sigmay)
    # todo: partial operators for pure advection

