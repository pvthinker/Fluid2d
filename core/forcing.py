from param import Param
from numpy import exp,pi,sin,cos

class Forcing(Param):
    """ define the forcing """

    def __init__(self,param,grid):

        yr = (grid.yr-param.Ly*.5)/param.Ly

        sigma = .1

        #        self.forc = 2e-4*(yr/sigma)*exp( -yr**2 / (2*sigma**2) )*grid.msk
        #self.forc = -1e-4*cos(yr*pi)*grid.msk

        Q = 4e-1
        nh = param.nh
        nxl=grid.nxl
        self.forc = yr*0.        
        self.forc[-nh-1,:] = -Q
        self.forc[nh,:] = +Q
        self.forc *= grid.msk

    def add_forcing(self,x,t,dxdt):
        """ add the forcing term on x[0]=the vorticity """
        dxdt[4] += self.forc #* dt
        
