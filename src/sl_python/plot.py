import matplotlib.pyplot as plt
import numpy as np

def plot_2D_scalar_field(xcr: np.ndarray, ycr: np.ndarray, scalar_field: np.ndarray, levels: int):
    """Plot 2D scalar field

    Args:
        xc (np.ndarray): x coordinates
        yc (np.ndarray): y coordinares
        scalar_field (np.ndarray): scalar field over xc, yc
    """
    # Shape

    fig, ax = plt.subplots(nrows=1)

    # Vent vertical
    cnt = ax.contourf(
        xcr,
        ycr,
        scalar_field,
        vmin=-15,
        vmax=15,
        cmap="bwr",
        extend="neither",
        levels=levels,
    )
    fig.colorbar(cnt, ax=ax, extend="neither")

    return fig


    