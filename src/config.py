from dataclasses import dataclass
from build import dtype

import numpy as np

@dataclass
class Config:
    
    dt: dtype
    
    xmin: dtype
    xmax: dtype
    nx: int
        
    ymin: dtype
    ymax: dtype
    ny: int
    
    bcx_kind: int
    bcy_kind: int
    
    def __post_init__(self):
        
        self.dth = self.dt / 2
        self.indices()
        self.coordinates()
                
    def indices(self):
        
        i_indices = np.arange(self.nx)
        j_indices = np.arange(self.ny)
        
        I, J = np.meshgrid(i_indices, j_indices)
        self.I, self.J = I.T, J.T
        
    def coordinates(self):
    
        self.xc = np.linspace(self.xmin, self.xmax, self.nx)
        self.yc = np.linspace(self.ymin, self.ymax, self.ny)
        
        self.dx = self.xc[1] - self.xc[0]
        self.dy = self.yc[1] - self.yc[0]
        
        xcr, ycr = np.meshgrid(self.xc, self.yc)
        self.xcr, self.ycr = xcr.T, ycr.T 
        
        
        
        