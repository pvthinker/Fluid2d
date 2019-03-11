from param import Param
import numpy as np


class Forcing(Param):
    """ define the forcing """

    def __init__(self, param, grid):

        yr = (grid.yr-param.Ly*.5)/param.Ly

        sigma = .1  # jet width relative to Ly

        tau0 = 1e-4

        # jet-like forcing (localized in y)
        # uncomment the forcing below
        # self.forc = tau0*(yr/sigma)*np.exp(-yr**2 / (2*sigma**2))*grid.msk

        # basin-scale forcing: double gyre configuration
        self.forc = tau0*np.sin(yr*np.pi)*grid.msk

        total = grid.domain_integration(self.forc)

        self.forc -= (total/grid.area)*grid.msk

    def add_forcing(self, x, t, dxdt):
        """ add the forcing term on x[0]=the vorticity """
        dxdt[0] += self.forc
