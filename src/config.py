from dataclasses import dataclass, field
import logging

import numpy as np

@dataclass
class Config:
    
    dt: float
    dth: float = field(init=False)
    
    xmin: float
    xmax: float
    nx: int
        
    ymin: float
    ymax: float
    ny: int
    

    bcx_kind: int
    bcy_kind: int
    
    lsettls: bool
    
    model_starttime: float 
    model_endtime: float

    nz: int = field(default=1)

    filter: bool = field(default=False)

    nstep: float = field(init=False)
    nitmp: float = field(default=4)

    gt4py_backend : str = field(default="numpy")

    def __post_init__(self):
        
        self.dth = self.dt / 2
        self.nstep = np.ceil((self.model_endtime - self.model_starttime) / self.dt)

        logging.info(f"Time step dt : {self.dt:.06f} s")
        logging.info(f"N steps : {self.nstep:.06f} s")

        self.indices()
        self.coordinates()

    @property
    def domain(self):
        return self.nx, self.ny, self.nz

    def indices(self):
        
        self.i_indices = np.arange(self.nx)
        self.j_indices = np.arange(self.ny)
        
        I, J = np.meshgrid(self.i_indices, self.j_indices)
        self.I, self.J = I.T, J.T

    def coordinates(self):
    
        self.xc = np.linspace(self.xmin, self.xmax, self.nx)
        self.yc = np.linspace(self.ymin, self.ymax, self.ny)
        
        self.dx = self.xc[1] - self.xc[0]
        self.dy = self.yc[1] - self.yc[0]
        
        xcr, ycr = np.meshgrid(self.xc, self.yc)
        self.xcr, self.ycr = xcr.T, ycr.T

