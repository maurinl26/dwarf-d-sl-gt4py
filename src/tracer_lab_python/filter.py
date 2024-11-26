import logging
import numpy as np

from config import Config
from tracer_lab_python.interpolation import max_interpolator_2d, min_interpolator_2d
from tracer_lab_python.periodic_filters import periodic_overshoot_filter, periodic_undershoot_filter

logging.getLogger(__name__)

def filter_driver(
    tracer: np.ndarray,
    tracer_e: np.ndarray,
    dep_idx_x: np.ndarray,
    dep_idx_y: np.ndarray,
    config: Config,
) -> np.ndarray:
    """Performs undershoot and overshoot filters.

    Args:
        tracer (np.ndarray): tracer field 
        dep_idx_x (np.ndarray): index of departure point
        dep_idx_y (np.ndarray): index of depature point
        bcx_kind (int): boundary kind
        bcy_kind (int): boundary kind
        nx (int): number of steps
        ny (int): number of steps
    """
    
    bcx_kind, bcy_kind = config.bcx_kind, config.bcy_kind
    nx, ny = config.nx, config.ny
    
    tracer_sup = max_interpolator_2d(
        tracer, dep_idx_x, dep_idx_y, bcx_kind, bcy_kind, nx, ny
        )

    tracer_inf = min_interpolator_2d(
        tracer, dep_idx_x, dep_idx_y, bcx_kind, bcy_kind, nx, ny
        )

    if bcx_kind == 1 and bcy_kind == 1:
        tracer_e = periodic_overshoot_filter(
                tracer_e, tracer_sup, nx, ny
            )
        tracer_e = periodic_undershoot_filter(
                tracer_e, tracer_inf, nx, ny
            )

    if bcx_kind == 0 and bcy_kind == 0:
        tracer_e = overshoot_filter(tracer_e, tracer_sup, nx, ny)
        tracer_e = undershoot_filter(tracer_e, tracer_inf, nx, ny)
        
    return tracer_e


def overshoot_filter(
    tracer_e: np.ndarray, tracer_sup: np.ndarray, nx: int, ny: int
):
    """Filter overshoots on interpolated field
    based on previous field.

    Args:
        tracer_e (np.ndarray): Grid point tracer at t + dt
        tracer (np.ndarray): Grid point tracer at t
        i_d (int): departure point index
        j_d (int): departure point index

        nx (int): _description_
        ny (int): _description_
    """

    # Sup du champ autour du point interpolé
    # Noyaux d'interpolation max
    
    

    # Surplus
    overshoot = np.maximum(tracer_e - tracer_sup, 0)

    ###### Bottom left corner
    i, j = 0, 0
    slot_4 = max(tracer_sup[i, j] - tracer_e[i + 1, j], 0)
    slot_5 = max(tracer_sup[i, j] - tracer_e[i + 1, j + 1], 0)
    slot_6 = max(tracer_sup[i, j] - tracer_e[i, j + 1], 0)

    total_disp = slot_4 + slot_5 + slot_6

    if total_disp > overshoot[i, j]:
        tracer_e[i, j] -= overshoot[i, j]
        tracer_e[i + 1, j] += (slot_4 / total_disp) * overshoot[i, j]
        tracer_e[i + 1, j + 1] += (slot_5 / total_disp) * overshoot[i, j]
        tracer_e[i, j + 1] += (slot_6 / total_disp) * overshoot[i, j]

    ####### Bottom right corner
    # Distribution en cas d'overshoot
    i, j = nx - 1, 0

    slot_6 = max(tracer_sup[i, j] - tracer_e[i, j + 1], 0)
    slot_7 = max(tracer_sup[i, j] - tracer_e[i - 1, j + 1], 0)
    slot_8 = max(tracer_sup[i, j] - tracer_e[i - 1, j], 0)

    total_disp = slot_6 + slot_7 + slot_8

    if total_disp > overshoot[i, j]:
        tracer_e[i, j] -= overshoot[i, j]
        tracer_e[i, j + 1] += (slot_6 / total_disp) * overshoot[i, j]
        tracer_e[i - 1, j + 1] += (slot_7 / total_disp) * overshoot[i, j]
        tracer_e[i - 1, j] += (slot_8 / total_disp) * overshoot[i, j]

    ######## Upper left corner
    i, j = 0, ny - 1
    slot_2 = max(tracer_sup[i, j] - tracer_e[i, j - 1], 0)
    slot_3 = max(tracer_sup[i, j] - tracer_e[i + 1, j - 1], 0)
    slot_4 = max(tracer_sup[i, j] - tracer_e[i + 1, j], 0)

    total_disp = slot_2 + slot_3 + slot_4

    if total_disp > overshoot[i, j]:
        tracer_e[i, j] -= overshoot[i, j]
        tracer_e[i, j - 1] += (slot_2 / total_disp) * overshoot[i, j]
        tracer_e[i + 1, j - 1] += (slot_3 / total_disp) * overshoot[i, j]
        tracer_e[i + 1, j] += (slot_4 / total_disp) * overshoot[i, j]

    ######## Upper right corner
    i, j = nx - 1, ny - 1
    # Distribution en cas d'overshoot
    slot_1 = max(tracer_sup[i, j] - tracer_e[i - 1, j - 1], 0)
    slot_2 = max(tracer_sup[i, j] - tracer_e[i, j - 1], 0)
    slot_8 = max(tracer_sup[i, j] - tracer_e[i - 1, j], 0)

    total_disp = slot_1 + slot_2 + slot_8

    if total_disp > overshoot[i, j]:
        tracer_e[i, j] -= overshoot[i, j]
        tracer_e[i - 1, j - 1] += (slot_1 / total_disp) * overshoot[i, j]
        tracer_e[i, j - 1] += (slot_2 / total_disp) * overshoot[i, j]
        tracer_e[i - 1, j] += (slot_8 / total_disp) * overshoot[i, j]

    ####### left border
    i = 0
    for j in range(1, ny - 1):
        slot_2 = max(tracer_sup[i, j] - tracer_e[i, j - 1], 0)
        slot_3 = max(tracer_sup[i, j] - tracer_e[i + 1, j - 1], 0)
        slot_4 = max(tracer_sup[i, j] - tracer_e[i + 1, j], 0)
        slot_5 = max(tracer_sup[i, j] - tracer_e[i + 1, j + 1], 0)
        slot_6 = max(tracer_sup[i, j] - tracer_e[i, j + 1], 0)

        total_disp = slot_2 + slot_3 + slot_4 + slot_5 + slot_6

        if total_disp > overshoot[i, j]:
            tracer_e[i, j] -= overshoot[i, j]
            tracer_e[i, j - 1] += (slot_2 / total_disp) * overshoot[i, j]
            tracer_e[i + 1, j - 1] += (slot_3 / total_disp) * overshoot[i, j]
            tracer_e[i + 1, j] += (slot_4 / total_disp) * overshoot[i, j]
            tracer_e[i + 1, j + 1] += (slot_5 / total_disp) * overshoot[i, j]
            tracer_e[i, j + 1] += (slot_6 / total_disp) * overshoot[i, j]

    ######## Right border
    i = nx - 1
    for j in range(1, ny - 1):
        slot_1 = max(tracer_sup[i, j] - tracer_e[i - 1, j - 1], 0)
        slot_2 = max(tracer_sup[i, j] - tracer_e[i, j - 1], 0)
        slot_6 = max(tracer_sup[i, j] - tracer_e[i, j + 1], 0)
        slot_7 = max(tracer_sup[i, j] - tracer_e[i - 1, j + 1], 0)
        slot_8 = max(tracer_sup[i, j] - tracer_e[i - 1, j], 0)

        total_disp = slot_1 + slot_2 + slot_6 + slot_7 + slot_8

        if total_disp > overshoot[i, j]:
            tracer_e[i, j] -= overshoot[i, j]
            tracer_e[i - 1, j - 1] += (slot_1 / total_disp) * overshoot[i, j]
            tracer_e[i, j - 1] += (slot_2 / total_disp) * overshoot[i, j]
            tracer_e[i, j + 1] += (slot_6 / total_disp) * overshoot[i, j]
            tracer_e[i - 1, j + 1] += (slot_7 / total_disp) * overshoot[i, j]
            tracer_e[i - 1, j] += (slot_8 / total_disp) * overshoot[i, j]

    ####### Bottom border
    j = 0
    for i in range(1, nx - 1):
        slot_4 = max(tracer_sup[i, j] - tracer_e[i + 1, j], 0)
        slot_5 = max(tracer_sup[i, j] - tracer_e[i + 1, j + 1], 0)
        slot_6 = max(tracer_sup[i, j] - tracer_e[i, j + 1], 0)
        slot_7 = max(tracer_sup[i, j] - tracer_e[i - 1, j + 1], 0)
        slot_8 = max(tracer_sup[i, j] - tracer_e[i - 1, j], 0)

        total_disp = slot_4 + slot_5 + slot_6 + slot_7 + slot_8

        if total_disp > overshoot[i, j]:
            tracer_e[i, j] -= overshoot[i, j]

            tracer_e[i + 1, j] += (slot_4 / total_disp) * overshoot[i, j]
            tracer_e[i + 1, j + 1] += (slot_5 / total_disp) * overshoot[i, j]
            tracer_e[i, j + 1] += (slot_6 / total_disp) * overshoot[i, j]
            tracer_e[i - 1, j + 1] += (slot_7 / total_disp) * overshoot[i, j]
            tracer_e[i - 1, j] += (slot_8 / total_disp) * overshoot[i, j]

    ####### Upper border
    j = ny - 1
    for i in range(1, nx - 1):
        slot_1 = max(tracer_sup[i, j] - tracer_e[i - 1, j - 1], 0)
        slot_2 = max(tracer_sup[i, j] - tracer_e[i, j - 1], 0)
        slot_3 = max(tracer_sup[i, j] - tracer_e[i + 1, j - 1], 0)
        slot_4 = max(tracer_sup[i, j] - tracer_e[i + 1, j], 0)
        slot_8 = max(tracer_sup[i, j] - tracer_e[i - 1, j], 0)

        total_disp = slot_1 + slot_2 + slot_3 + slot_4 + slot_8

        if total_disp > overshoot[i, j]:
            tracer_e[i, j] -= overshoot[i, j]

            tracer_e[i - 1, j - 1] += (slot_1 / total_disp) * overshoot[i, j]
            tracer_e[i, j - 1] += (slot_2 / total_disp) * overshoot[i, j]
            tracer_e[i + 1, j - 1] += (slot_3 / total_disp) * overshoot[i, j]
            tracer_e[i + 1, j] += (slot_4 / total_disp) * overshoot[i, j]
            tracer_e[i - 1, j] += (slot_8 / total_disp) * overshoot[i, j]

    # Inner domain
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            # Distribution en cas d'overshoot
            slot_1 = max(tracer_sup[i, j] - tracer_e[i - 1, j - 1], 0)
            slot_2 = max(tracer_sup[i, j] - tracer_e[i, j - 1], 0)
            slot_3 = max(tracer_sup[i, j] - tracer_e[i + 1, j - 1], 0)
            slot_4 = max(tracer_sup[i, j] - tracer_e[i + 1, j], 0)
            slot_5 = max(tracer_sup[i, j] - tracer_e[i + 1, j + 1], 0)
            slot_6 = max(tracer_sup[i, j] - tracer_e[i, j + 1], 0)
            slot_7 = max(tracer_sup[i, j] - tracer_e[i - 1, j + 1], 0)
            slot_8 = max(tracer_sup[i, j] - tracer_e[i - 1, j], 0)

            total_disp = (
                slot_1 + slot_2 + slot_3 + slot_4 + slot_5 + slot_6 + slot_7 + slot_8
            )

            # inner domain
            if total_disp > overshoot[i, j]:
                tracer_e[i, j] -= overshoot[i, j]
                tracer_e[i - 1, j - 1] += (slot_1 / total_disp) * overshoot[i, j]
                tracer_e[i, j - 1] += (slot_2 / total_disp) * overshoot[i, j]
                tracer_e[i + 1, j - 1] += (slot_3 / total_disp) * overshoot[i, j]
                tracer_e[i + 1, j] += (slot_4 / total_disp) * overshoot[i, j]
                tracer_e[i + 1, j + 1] += (slot_5 / total_disp) * overshoot[i, j]
                tracer_e[i, j + 1] += (slot_6 / total_disp) * overshoot[i, j]
                tracer_e[i - 1, j + 1] += (slot_7 / total_disp) * overshoot[i, j]
                tracer_e[i - 1, j] += (slot_8 / total_disp) * overshoot[i, j]
                
    return tracer_e
                
def undershoot_filter(
    tracer_e: np.ndarray, tracer_inf: np.ndarray, nx: int, ny: int
):
    """Filter undershoots on interpolated field
    based on previous field.

    Args:
        tracer_e (np.ndarray): Grid point tracer at t + dt
        tracer (np.ndarray): Grid point tracer at t
        i_d (int): departure point index
        j_d (int): departure point index

        nx (int): _description_
        ny (int): _description_
    """
    
    # 
    undershoot = np.maximum(tracer_inf - tracer_e, np.finfo(np.float64).tiny)

    ###### Bottom left corner
    i, j = 0, 0
    slot_4 = max(tracer_e[i + 1, j] - tracer_e[i, j], 0)
    slot_5 = max(tracer_e[i + 1, j + 1] - tracer_e[i, j], 0)
    slot_6 = max(tracer_e[i, j + 1] - tracer_e[i, j], 0)

    total_disp = slot_4 + slot_5 + slot_6

    if total_disp > undershoot[i, j]:
        tracer_e[i, j]  += undershoot[i, j]
        tracer_e[i + 1, j] -= (slot_4 / total_disp) * undershoot[i, j]
        tracer_e[i + 1, j + 1] -= (slot_5 / total_disp) * undershoot[i, j]
        tracer_e[i, j + 1] -= (slot_6 / total_disp) * undershoot[i, j]

    ####### Bottom right corner
    # Distribution en cas d'undershoot
    i, j = nx - 1, 0

    slot_6 = max(tracer_e[i, j + 1] - tracer_inf[i, j], 0)
    slot_7 = max(tracer_e[i - 1, j + 1] - tracer_inf[i, j], 0)
    slot_8 = max(tracer_e[i - 1, j] - tracer_inf[i, j], 0)

    total_disp = slot_6 + slot_7 + slot_8

    if total_disp > undershoot[i, j]:
        tracer_e[i, j]  += undershoot[i, j]
        tracer_e[i, j + 1] -= (slot_6 / total_disp) * undershoot[i, j]
        tracer_e[i - 1, j + 1] -= (slot_7 / total_disp) * undershoot[i, j]
        tracer_e[i - 1, j] -= (slot_8 / total_disp) * undershoot[i, j]

    ######## Upper left corner
    i, j = 0, ny - 1
    slot_2 = max(tracer_e[i, j - 1] - tracer_inf[i, j], 0)
    slot_3 = max(tracer_e[i + 1, j - 1] - tracer_inf[i, j], 0)
    slot_4 = max(tracer_e[i + 1, j] - tracer_inf[i, j], 0)

    total_disp = slot_2 + slot_3 + slot_4

    if total_disp > undershoot[i, j]:
        tracer_e[i, j]  += undershoot[i, j]
        tracer_e[i, j - 1] -= (slot_2 / total_disp) * undershoot[i, j]
        tracer_e[i + 1, j - 1] -= (slot_3 / total_disp) * undershoot[i, j]
        tracer_e[i + 1, j] -= (slot_4 / total_disp) * undershoot[i, j]

    ######## Upper right corner
    i, j = nx - 1, ny - 1
    # Distribution en cas d'undershoot
    slot_1 = max(tracer_e[i - 1, j - 1] - tracer_inf[i, j], 0)
    slot_2 = max(tracer_e[i, j - 1] - tracer_inf[i, j], 0)
    slot_8 = max(tracer_e[i - 1, j] - tracer_inf[i, j], 0)

    total_disp = slot_1 + slot_2 + slot_8

    if total_disp > undershoot[i, j]:
        tracer_e[i, j]  += undershoot[i, j]
        tracer_e[i - 1, j - 1] -= (slot_1 / total_disp) * undershoot[i, j]
        tracer_e[i, j - 1] -= (slot_2 / total_disp) * undershoot[i, j]
        tracer_e[i - 1, j] -= (slot_8 / total_disp) * undershoot[i, j]

    ####### left border
    i = 0
    for j in range(1, ny - 1):
        slot_2 = max(tracer_e[i, j - 1] - tracer_inf[i, j], 0)
        slot_3 = max(tracer_e[i + 1, j - 1] - tracer_inf[i, j], 0)
        slot_4 = max(tracer_e[i + 1, j] - tracer_inf[i, j], 0)
        slot_5 = max(tracer_e[i + 1, j + 1] - tracer_inf[i, j], 0)
        slot_6 = max(tracer_e[i, j + 1] - tracer_inf[i, j], 0)

        total_disp = slot_2 + slot_3 + slot_4 + slot_5 + slot_6

        if total_disp > undershoot[i, j]:
            tracer_e[i, j]  += undershoot[i, j]
            tracer_e[i, j - 1] -= (slot_2 / total_disp) * undershoot[i, j]
            tracer_e[i + 1, j - 1] -= (slot_3 / total_disp) * undershoot[i, j]
            tracer_e[i + 1, j] -= (slot_4 / total_disp) * undershoot[i, j]
            tracer_e[i + 1, j + 1] -= (slot_5 / total_disp) * undershoot[i, j]
            tracer_e[i, j + 1] -= (slot_6 / total_disp) * undershoot[i, j]

    ######## Right border
    i = nx - 1
    for j in range(1, ny - 1):
        slot_1 = max(tracer_e[i - 1, j - 1] - tracer_inf[i, j], 0)
        slot_2 = max(tracer_e[i, j - 1] - tracer_inf[i, j], 0)
        slot_6 = max(tracer_e[i, j + 1] - tracer_inf[i, j], 0)
        slot_7 = max(tracer_e[i - 1, j + 1] - tracer_inf[i, j], 0)
        slot_8 = max(tracer_e[i - 1, j] - tracer_inf[i, j], 0)

        total_disp = slot_1 + slot_2 + slot_6 + slot_7 + slot_8

        if total_disp > undershoot[i, j]:
            tracer_e[i, j]  += undershoot[i, j]
            tracer_e[i - 1, j - 1] -= (slot_1 / total_disp) * undershoot[i, j]
            tracer_e[i, j - 1] -= (slot_2 / total_disp) * undershoot[i, j]
            tracer_e[i, j + 1] -= (slot_6 / total_disp) * undershoot[i, j]
            tracer_e[i - 1, j + 1] -= (slot_7 / total_disp) * undershoot[i, j]
            tracer_e[i - 1, j] -= (slot_8 / total_disp) * undershoot[i, j]

    ####### Bottom border
    j = 0
    for i in range(1, nx - 1):
        slot_4 = max(tracer_e[i + 1, j] - tracer_inf[i, j], 0)
        slot_5 = max(tracer_e[i + 1, j + 1] - tracer_inf[i, j], 0)
        slot_6 = max(tracer_e[i, j + 1] - tracer_inf[i, j], 0)
        slot_7 = max(tracer_e[i - 1, j + 1] - tracer_inf[i, j], 0)
        slot_8 = max(tracer_e[i - 1, j] - tracer_inf[i, j], 0)

        total_disp = slot_4 + slot_5 + slot_6 + slot_7 + slot_8

        if total_disp > undershoot[i, j]:
            tracer_e[i, j]  += undershoot[i, j]
            tracer_e[i + 1, j] -= (slot_4 / total_disp) * undershoot[i, j]
            tracer_e[i + 1, j + 1] -= (slot_5 / total_disp) * undershoot[i, j]
            tracer_e[i, j + 1] -= (slot_6 / total_disp) * undershoot[i, j]
            tracer_e[i - 1, j + 1] -= (slot_7 / total_disp) * undershoot[i, j]
            tracer_e[i - 1, j] -= (slot_8 / total_disp) * undershoot[i, j]

    ####### Upper border
    j = ny - 1
    for i in range(1, nx - 1):
        slot_1 = max(tracer_e[i - 1, j - 1] - tracer_inf[i, j], 0)
        slot_2 = max(tracer_e[i, j - 1] - tracer_inf[i, j], 0)
        slot_3 = max(tracer_e[i + 1, j - 1] - tracer_inf[i, j], 0)
        slot_4 = max(tracer_e[i + 1, j] - tracer_inf[i, j], 0)
        slot_8 = max(tracer_e[i - 1, j] - tracer_inf[i, j], 0)

        total_disp = slot_1 + slot_2 + slot_3 + slot_4 + slot_8

        if total_disp > undershoot[i, j]:
            tracer_e[i, j]  += undershoot[i, j]
            tracer_e[i - 1, j - 1] -= (slot_1 / total_disp) * undershoot[i, j]
            tracer_e[i, j - 1] -= (slot_2 / total_disp) * undershoot[i, j]
            tracer_e[i + 1, j - 1] -= (slot_3 / total_disp) * undershoot[i, j]
            tracer_e[i + 1, j] -= (slot_4 / total_disp) * undershoot[i, j]
            tracer_e[i - 1, j] -= (slot_8 / total_disp) * undershoot[i, j]

    # Inner domain
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            # Distribution en cas d'undershoot
            slot_1 = max(tracer_e[i - 1, j - 1] - tracer_inf[i, j], 0)
            slot_2 = max(tracer_e[i, j - 1] - tracer_inf[i, j], 0)
            slot_3 = max(tracer_e[i + 1, j - 1] - tracer_inf[i, j], 0)
            slot_4 = max(tracer_e[i + 1, j] - tracer_inf[i, j], 0)
            slot_5 = max(tracer_e[i + 1, j + 1] - tracer_inf[i, j], 0)
            slot_6 = max(tracer_e[i, j + 1] - tracer_inf[i, j], 0)
            slot_7 = max(tracer_e[i - 1, j + 1] - tracer_inf[i, j], 0)
            slot_8 = max(tracer_e[i - 1, j] - tracer_inf[i, j], 0)

            total_disp = (
                slot_1 + slot_2 + slot_3 + slot_4 + slot_5 + slot_6 + slot_7 + slot_8
            )

            # inner domain
            if total_disp > undershoot[i, j]:
                tracer_e[i, j]  += undershoot[i, j]
                tracer_e[i - 1, j - 1] -= (slot_1 / total_disp) * undershoot[i, j]
                tracer_e[i, j - 1] -= (slot_2 / total_disp) * undershoot[i, j]
                tracer_e[i + 1, j - 1] -= (slot_3 / total_disp) * undershoot[i, j]
                tracer_e[i + 1, j] -= (slot_4 / total_disp) * undershoot[i, j]
                tracer_e[i + 1, j + 1] -= (slot_5 / total_disp) * undershoot[i, j]
                tracer_e[i, j + 1] -= (slot_6 / total_disp) * undershoot[i, j]
                tracer_e[i - 1, j + 1] -= (slot_7 / total_disp) * undershoot[i, j]
                tracer_e[i - 1, j] -= (slot_8 / total_disp) * undershoot[i, j]
                
    return tracer_e
