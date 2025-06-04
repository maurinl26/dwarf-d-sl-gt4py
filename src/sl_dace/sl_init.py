from typing import Union

import dace
from sl_dace.utils.sdfg import build_sdfg
from sl_dace.stencils.sl_init import settls_init, nesc_init

class SLInit:

    def __init__(self, option = Union["lsettls", "lnesc"]):

        match option:
            case "lsettls":
                self.d_sl_init = build_sdfg(nesc_init)
            case "lnesc":
                self.sl_init = build_sdfg(settls_init)

    @dace.method
    def __call__(self,
                 state: dict):


        ...

