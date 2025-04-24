from typing import Tuple, Union
from gt4py.cartesian.gtscript import stencil

from sl_dace.stencils.sl_init import settls_init, nesc_init

class SLInit:

    def __init__(self, option = Union["lsettls", "lnesc"]):

        match option:
            case "lsettls":
                self.sl_init = stencil(
                    backend="dace:cpu",
                    definition=nesc_init,
                    name="nesc"
                )
            case "lnesc":
                self.sl_init = stencil(
                    backend="dace:cpu",
                    definition=settls_init,
                    name="settls"
                 )

    # todo: implement call
    def __call__(self):
        ...

