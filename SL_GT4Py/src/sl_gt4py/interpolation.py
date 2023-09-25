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
            l / (x[0, 0, 0] - x[-1, 0, 0])  -1 ) * (
            l / (x[0, 0, 0] - x[-2, 0, 0]) - 1 ) * (
            l / (x[0, 0, 0] - x[1, 0, 0])  - 1 )
            
    p2 = psi[1, 0, 0] * ( 
            l / (x[1, 0, 0] - x[-1, 0, 0])  -1 ) * ( 
            l / (x[1, 0, 0] - x[-2, 0, 0]) - 1 ) * (
            l / (x[1], 0, 0] - x[0, 0, 0])  - 1 )
    p3 = psi[-1, 0, 0] * ( 
            l / (x[-1, 0, 0] - x[0, 0, 0])  -1 ) * (
            l / (x[-1, 0, 0] - x[-2, 0, 0]) - 1 ) * (
            l / (x[-1, 0, 0] - x[1, 0, 0])  - 1 )
    p4 = psi[-2, 0, 0] * ( 
            l / (x[-2, 0, 0] - x[-1, 0, 0])  -1 ) * (
            l / (x[-2, 0, 0] - x[0, 0, 0]) - 1 ) * (
            l / (x[-2, 0, 0] - x[1, 0, 0])  - 1 )
            
    l3 = p1 + p2 + p3 + p4
    
#  TODO :  implementer l'interpolation lineaire a gauche (0, 1)
@gtscript.function
def lagrange_linear_x(
    l1: gtscript.Field[np.float64],
    psi: gtscript.Field[np.float64],
    x: gtscript.Field[np.float64],
    l: np.float64
):
    l1 = psi[0, 0, 0] * ( l / (x[-1, 0, 0] - x[0, 0, 0]) - 1)

