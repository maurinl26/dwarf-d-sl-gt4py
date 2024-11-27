from gt4py.storage import from_array, ones, zeros
import numpy as np
from tracer_lab_gt4py.interpolation import elarmes_1d_linear 

if __name__ == "__main__":
    
    backend = "gt:cpu_ifirst"
    
    nx, ny, nz = 50, 50, 1
    
    xmin, xmax = 0, 200
    ymin, ymax = 0, 200
    
    shape = (nx, ny, nz)
    origin = (0, 0, 0)
    
    
    
    xc = np.arange(xmin, xmax, nx)
    yc = np.arange(ymin, ymax, ny)
    
    # Grid coordinates
    xcr, ycr = np.meshgrid(xc, yc, indexing="ij")
    
    # Grid indexing
    i_axis = np.arange(0, nx)
    j_axis = np.arange(0, ny)
    
    idx_x, idx_y = np.meshgrid(i_axis, j_axis, indexing="ij")
    
    #############################################
    ######## Numpy fields #######################
    #############################################
    
    tracer0 = np.zeros((nx, ny)) # just one level
    tracer9 = np.zeros((nx, ny, nz))
    weight_x_dep = 0.5 * np.ones((nx, ny, nz))
    weight_y_dep = 0.5 * np.ones((nx, ny, nz))
    idx_x_dep = idx_x
    idx_y_dep = idx_y
    
    
    #############################################
    ########## GT4Py storage ####################
    #############################################
    tracer0_in =  from_array(tracer0, np.float64, backend=backend, aligned_index=(0, 0))
    tracer9_out = from_array(tracer9, np.float64, backend=backend, aligned_index=origin)
    weight_x_dep_in = from_array(weight_x_dep, np.float64, backend=backend, aligned_index=origin)
    weight_y_dep_in = from_array(weight_y_dep, np.float64, backend=backend, aligned_index=origin)
    idx_x_dep_in = from_array(idx_x_dep, int, backend=backend, aligned_index=(0,0))
    idx_y_dep_in = from_array(idx_y_dep, int, backend=backend, aligned_index=(0,0))
    
    
    #############################################
    ######### Stencil call  #####################
    #############################################
    elarmes_1d_linear(
        tracer_e=tracer9_out,
        tracer=tracer0_in, 
        weight_x=weight_x_dep_in,
        weight_y=weight_y_dep_in,
        idx_dep=idx_x_dep_in,
        idy_dep=idx_y_dep_in,   
        origin=origin,
    )
    
    
    
    
    