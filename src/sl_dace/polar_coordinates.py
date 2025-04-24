import numpy as np
from typing import Tuple
from gt4py.cartesian.gtscript import stencil

from sl_dace.stencils.blossey import radius, theta

class PolarCoordinates:

    def __init__(self, grid: Tuple[int]):
        self.grid = grid

        # X and Y axis
        self.xcr = np.arange(0, self.grid[0])
        self.ycr = np.arange(0, self.grid[1])

        # Horizontal coordinates on a grid
        self.horizontal_coordinates = np.meshgrid((self.xcr, self.ycr))

        # Stencils
        self.radius  = stencil(
            backend="dace:cpu",
            definition=radius,
            name="radius"
        )
        self.theta = stencil(
            backend="dace:cpu",
            definition=theta,
            name="theta"
        )

    def __call__(self):
        self.radius(self.xcr, self.ycr, radius, domain=self.grid, origin=(0, 0, 0))
        self.theta(self.xcr, self.ycr, theta, domain=self.grid, origin=(0, 0, 0))