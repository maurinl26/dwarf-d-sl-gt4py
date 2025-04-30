from collections.abc import dict_values
from ifs_physics_common.framework.storage import managed_temporary_storage
from config import Config


class FluxFormSemiLagXY:

    def __init__(self,
                 config: Config,
                 ):



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

