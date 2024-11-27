from tracer_lab_gt4py.slag_2D_xy import slag_search
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
    
    idx_x = idx_x[:,:, np.newaxis]
    idx_y = idx_y[:,:, np.newaxis]
        
    #############################################
    ######## Numpy fields #######################
    #############################################
    
    
    weight_x_dep = 0.5 * np.ones((nx, ny, nz))
    weight_y_dep = 0.5 * np.ones((nx, ny, nz))
    idx_x_dep = idx_x.copy()
    idx_y_dep = idx_y.copy()
    vx = np.zeros((nx, ny, nz))
    vx_e = np.zeros((nx, ny, nz))
    vy = np.zeros((nx, ny, nz))
    vy_e = np.zeros((nx, ny, nz))
    
    
    #############################################
    ########## GT4Py storage ####################
    #############################################

    weight_x_dep_out = from_array(weight_x_dep, np.float64, backend=backend, aligned_index=origin)
    weight_y_dep_out = from_array(weight_y_dep, np.float64, backend=backend, aligned_index=origin)
    idx_x_dep_out = from_array(idx_x_dep, np.int32, backend=backend, aligned_index=origin)
    idx_y_dep_out = from_array(idx_y_dep, np.int32, backend=backend, aligned_index=origin)
    
    idx_x_arr_in = from_array(idx_x_dep, np.int32, backend=backend, aligned_index=origin)
    idx_y_arr_in = from_array(idx_y_dep, np.int32, backend=backend, aligned_index=origin)
    
    vx_e_in = from_array(vx_e, np.float64, backend=backend, aligned_index=origin)
    vy_e_in = from_array(vy_e, np.float64, backend=backend, aligned_index=origin)
    vx_tmp_in = from_array(vx, np.float64, backend=backend, aligned_index=origin)
    vy_tmp_in = from_array(vy, np.float64, backend=backend, aligned_index=origin)
    
    
    #############################################
    ######### Stencil call  #####################
    #############################################
    slag_search(
        weight_x_dep=weight_x_dep_out, 
        weight_y_dep=weight_y_dep_out, 
        idx_dep=idx_x_dep_out,
        idy_dep=idx_y_dep_out,
        vx_e=vx_e_in,
        vx_tmp=vx_tmp_in,
        vy_e=vy_e_in,
        vy_tmp=vy_tmp_in,
        idx_arr=idx_x_arr_in,
        idy_arr=idx_y_arr_in,
        dx=1.0,
        dy=1.0,
        dth=0.5
    )