from param import Param
from numpy import exp,pi,sin,cos

class Forcing(Param):
    """ define the forcing """

    def __init__(self,param,grid):

        yr = (grid.yr-param.Ly*.5)/param.Ly

        sigma = .1 # jet width relative to Ly

        tau0 = 1e-4

        # jet-like forcing (localized in y)
        self.forc = tau0*(yr/sigma)*exp( -yr**2 / (2*sigma**2) )*grid.msk
        
        # basin-scale forcing: double gyre configuration
        #
        # uncomment the forcing below
        #
        #self.forc = tau0*sin(yr*pi)*grid.msk


        total = grid.domain_integration(self.forc)
        
        self.forc -= (total/grid.area)*grid.msk



    def add_forcing(self,x,t,dxdt):
        """ add the forcing term on x[0]=the vorticity """
        dxdt[0] += self.forc 
        

