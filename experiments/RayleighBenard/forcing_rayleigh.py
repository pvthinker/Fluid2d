from param import Param
from numpy import exp,pi,sin,cos

class Forcing(Param):
    """ define the forcing """

    def __init__(self,param,grid):

        yr = (grid.yr-param.Ly*.5)/param.Ly

        Q = 1e-2
        nh = param.nh
        nxl=grid.nxl
        self.forc = yr*0.

        type_forcing = 'coolroof'#'horizontalconvection'

        if type_forcing=='coolroof':
            if grid.j0 == grid.npy-1:
                self.forc[-nh-1,:] = -Q
            if grid.j0 == 0:
                self.forc[nh,:] = +Q

        if type_forcing=='spreadcooling':
            # bottom heating
            if grid.j0 == 0:
                self.forc[nh,:] = +Q
            # uniform cooling
            self.forc -= Q/param.ny

        if type_forcing == 'horizontalconvection':
            if grid.j0 == grid.npy-1:
                xr = grid.xr
                self.forc[xr<param.Lx/2] = Q
                self.forc[xr>param.Lx/2] = -Q
                self.forc[:-nh-1,:]=0.
            

        
        self.forc *= grid.msk

        # transform the surface flux into a volume flux
        self.forc *= (1./grid.dx) 

    def add_forcing(self,x,t,dxdt):
        """ add the forcing term on x[0]=the vorticity """
        dxdt[4] += self.forc 
        
