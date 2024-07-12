import numpy as np

dtype = np.float64
dtype_int = np.int32

backend = "gt:cpu_ifirst"
backend_opts = {"verbose": True} if backend.startswith("gt") else {}
origin = (0, 0, 0)
rebuild = False

# Functions
externals = {
    
}