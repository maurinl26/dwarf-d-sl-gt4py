import numpy as np
import matplotlib.pyplot as plt
from sl_python.interpolation import interpolate_cubic_2d

def diagnostic_interpolation(field: np.ndarray, i: int, j: int, j_step: int, dx: float):
    """Trace les polynomes de lagrange associés

    Args:
        tracer (np.ndarray): _description_
        i (int): _description_
        j (int): _description_
    """
    
    lx = np.linspace(-1, 2, 30)
    ly = np.linspace(-1, 2, 30)
    
    # Polynomes de lagrange d'ordre 3
    p_1 = lambda l: (1/6) * l * (l - 1) * (2 - l) 
    p0 = lambda l: (1 - l**2) * (1 - l / 2)
    p1 = lambda l: (1/2) * l * (l + 1) * (2 - l)
    p2 = lambda l: (1/6) * l * (l**2 - 1) 
    
    g = lambda l: 20 * np.exp(-(l + (i - 6 - j_step * 1.5))**2 / 4)

    p0_x = [p_1(l)*field[i-1, j] + p0(l)*field[i, j] + p1(l)*field[i+1, j] + p2(l)*field[i+2, j] for l in lx]
    g_x = [g(l) for l in lx]
    
    plt.plot(lx, p0_x)
    plt.plot(lx, g_x)
    plt.show()