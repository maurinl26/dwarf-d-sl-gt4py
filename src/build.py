import numpy as np

dtype = np.float64

backend = "numpy"
backend_opts = {"verbose": True} if backend.startswith("gt") else {}
origin = (0, 0, 0)
rebuild = False

# Functions
externals = {
    
}