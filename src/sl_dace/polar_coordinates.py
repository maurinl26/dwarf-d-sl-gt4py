import numpy as np
from config import Config

import dace
from sl_dace.utils.sdfg import build_sdfg
from sl_dace.stencils.blossey import radius, theta

class PolarCoordinates:

    def __init__(self, config: Config):
        self.config = config
        self.grid = self.config.domain

        # X and Y axis
        self.xcr = np.arange(0, self.grid[0])
        self.ycr = np.arange(0, self.grid[1])

        # Horizontal coordinates on a grid
        self.horizontal_coordinates = np.meshgrid((self.xcr, self.ycr))

        # Stencils
        self.d_radius  = build_sdfg(radius)
        self.d_theta = build_sdfg(theta)


    def __call__(self,
                 radius: np.ndarray,
                 theta: np.ndarray
                 ):
        self.d_radius(self.xcr, self.ycr, radius, domain=self.grid, origin=(0, 0, 0))
        self.d_theta(self.xcr, self.ycr, theta, domain=self.grid, origin=(0, 0, 0))