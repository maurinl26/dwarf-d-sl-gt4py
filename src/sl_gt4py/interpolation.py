from typing import Tuple
from gt4py.cartesian import gtscript
import numpy as np

# TODO : implementer l'interpolation lineaire à gauche (-1, 0, 1, 2)
@gtscript.function
def lagrange_cubic_x(
    l3: gtscript.Field[np.float64],
    psi: gtscript.Field[np.float64],
    x: gtscript.Field[np.float64],
    l: np.float64
):
    p1 = psi[0, 0, 0] * ( 
            l / (x[0, 0, 0] - x[-1, 0, 0]) - 1 ) * (
            l / (x[0, 0, 0] - x[-2, 0, 0]) - 1 ) * (
            l / (x[0, 0, 0] - x[1, 0, 0])  - 1 )
            
    p2 = psi[1, 0, 0] * ( 
            l / (x[1, 0, 0] - x[-1, 0, 0])  -1 ) * ( 
            l / (x[1, 0, 0] - x[-2, 0, 0]) - 1 ) * (
            l / (x[1, 0, 0] - x[0, 0, 0])  - 1 )
    p3 = psi[-1, 0, 0] * ( 
            l / (x[-1, 0, 0] - x[0, 0, 0])  -1 ) * (
            l / (x[-1, 0, 0] - x[-2, 0, 0]) - 1 ) * (
            l / (x[-1, 0, 0] - x[1, 0, 0])  - 1 )
    p4 = psi[-2, 0, 0] * ( 
            l / (x[-2, 0, 0] - x[-1, 0, 0])  -1 ) * (
            l / (x[-2, 0, 0] - x[0, 0, 0]) - 1 ) * (
            l / (x[-2, 0, 0] - x[1, 0, 0])  - 1 )
            
    l3 = p1 + p2 + p3 + p4
    
@gtscript.function
def lagrange_cubic_y(
    l3: gtscript.Field[np.float64],
    psi: gtscript.Field[np.float64],
    x: gtscript.Field[np.float64],
    l: np.float64
):
    p1 = psi[0, 0, 0] * ( 
            l / (x[0, 0, 0] - x[0, -1, 0])  -1 ) * (
            l / (x[0, 0, 0] - x[0, -2, 0]) - 1 ) * (
            l / (x[0, 0, 0] - x[0, 1, 0])  - 1 )
            
    p2 = psi[0, 1, 0] * ( 
            l / (x[0, 1, 0] - x[0, -1, 0])  -1 ) * ( 
            l / (x[0, 1, 0] - x[0, -2, 0]) - 1 ) * (
            l / (x[0, 1, 0] - x[0, 0, 0])  - 1 )
    p3 = psi[0, -1, 0] * ( 
            l / (x[0, -1, 0] - x[0, 0, 0])  -1 ) * (
            l / (x[0, -1, 0] - x[0,-2, 0]) - 1 ) * (
            l / (x[0, -1, 0] - x[0, 1, 0])  - 1 )
    p4 = psi[0, -2, 0] * ( 
            l / (x[0, -2, 0] - x[0, -1, 0])  -1 ) * (
            l / (x[0, -2, 0] - x[0, 0, 0]) - 1 ) * (
            l / (x[0, -2, 0] - x[0, 1, 0])  - 1 )
            
    l3 = p1 + p2 + p3 + p4
    
@gtscript.function
def lagrange_cubic_z(
    l3: gtscript.Field[np.float64],
    psi: gtscript.Field[np.float64],
    x: gtscript.Field[np.float64],
    l: np.float64
):
    p1 = psi[0, 0, 0] * ( 
            l / (x[0, 0, 0] - x[0, 0, -1])  -1 ) * (
            l / (x[0, 0, 0] - x[0, 0, -2]) - 1 ) * (
            l / (x[0, 0, 0] - x[0, 0, 1])  - 1 )
            
    p2 = psi[0, 0, 1] * ( 
            l / (x[0, 0, 1] - x[0, 0, -1])  -1 ) * ( 
            l / (x[0, 0, 1] - x[0, 0, -2]) - 1 ) * (
            l / (x[0, 0, 1] - x[0, 0, 0])  - 1 )
    p3 = psi[0, 0, -1] * ( 
            l / (x[0, 0, -1] - x[0, 0, 0])  -1 ) * (
            l / (x[0, 0, -1] - x[0, 0, -2]) - 1 ) * (
            l / (x[0, 0, -1] - x[0, 0, 1])  - 1 )
    p4 = psi[0, 0, -2] * ( 
            l / (x[0, 0, -2] - x[0, 0, -1])  -1 ) * (
            l / (x[0, 0, -2] - x[0, 0, 0]) - 1 ) * (
            l / (x[0, 0, -2] - x[0, 0, 1])  - 1 )
            
    l3 = p1 + p2 + p3 + p4
    
@gtscript.function
def lagrange_cubic_z_stretched(
        l3: gtscript.Field[np.float64],
        psi: gtscript.Field[np.float64],
        x: gtscript.Field[np.float64],
        l: np.float64
):
        
        return NotImplemented
    
@gtscript.function
def lagrange_linear(
        lx: np.float64,
        ly: np.float64,
        lz: np.float64,
        psi: gtscript.Field[np.float64],
        X = Tuple[gtscript.Field]
):
        psi_x = lagrange_cubic_x(psi_tmp, psi, X[0], lx)
        psi_y = lagrange_cubic_y(psi_y, psi_x, X[1], ly)
        psi_z = lagrange_cubic_z(psi_z, psi, X[2], lz)

    
#  TODO :  implementer l'interpolation lineaire a gauche (0, 1)
@gtscript.function
def lagrange_linear_x(
    l1: gtscript.Field[np.float64],
    psi: gtscript.Field[np.float64],
    x: gtscript.Field[np.float64],
    l: np.float64
):
    l1 = psi[0, 0, 0] * ( l / (x[-1, 0, 0] - x[0, 0, 0]) - 1)

