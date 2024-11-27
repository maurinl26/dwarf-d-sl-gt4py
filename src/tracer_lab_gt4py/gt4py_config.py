
import numpy as np


# gridtools4py settings
backend = "numpy"  # options: "numpy", "gt:cpu_ifirst", "gt:cpu_kfirst", "gt:gpu", "dace:cpu", "dace:gpu"
backend_opts = {'verbose': True} if backend.startswith('gt') else {}
dtype_float = np.float64
dtype_int = np.int32
origin = (0, 0, 0)
rebuild = True

#####################
#### Danger zone ####
#####################
nx, ny = 50, 50
bcx_kind, bcy_kind = 1, 1
origin = (0, 0, 0)

externals = {
    "nx": nx,
    "ny": ny,
    "bcx_kind":bcx_kind,
    "bcy_kind":bcy_kind
}