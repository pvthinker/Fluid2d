from param import Param
from numpy import exp,pi,sin,cos, random

class Forcing(Param):
    """ define the forcing """

    def __init__(self,param,grid):
        self.list_param=['sizevar']
        param.copy(self,self.list_param)

        self.list_param=['nh','msk','fill_halo']
        grid.copy(self,self.list_param)

        self.intensity=1e-1
        self.gamma = 1e-1
        self.t=0

        self.forc = random.normal(size=self.sizevar)*self.msk
        self.fill_halo(self.forc)

    def add_forcing(self,x,t,dxdt):
        """ add the forcing term on x[0]=the vorticity """
        
        dt = t-self.t
        self.t=t

        self.forc = (1-self.gamma*dt)*self.forc + dt*self.gamma * random.normal(size=self.sizevar)*self.msk
        self.fill_halo(self.forc)

        dxdt[0] += self.intensity * self.forc
        
