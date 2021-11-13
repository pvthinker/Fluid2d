from param import Param
import numpy as np


class Forcing(Param):
    """ define the forcing """

    def __init__(self, param, grid):

        self.list_param = ['deltab', 'dt']
        param.copy(self, self.list_param)

        self.list_param = ['j0', 'npy', 'nh']
        grid.copy(self, self.list_param)

        xr = (grid.xr-param.Lx*.5)/param.Lx
        yr = (grid.yr-param.Ly*.5)/param.Ly

        nh = param.nh
        self.forc = yr*0.
        N2 = 1.
        angle = 25.  # in degrees
        # angle close to 0° => phase propagates horizontally,
        #                      beams are vertical
        # angle close to 90° => phase propagates vertically,
        #                      beams are horizontal
        self.period = 2*np.pi/(np.sqrt(N2)*np.cos(angle*np.pi/180))
        amplitude = .1
        spacing = 0.02
        width = .02
        y0 = 0.

        def gaussian(x0, y0, delta):
            d2 = (xr-x0)**2 + (yr-y0)**2
            return np.exp(-d2/(2*delta**2))

        self.forc[:, :] = (gaussian(-spacing/2, y0, width)
                           - gaussian(+spacing/2, y0, width))
        self.forc *= grid.msk*amplitude/width

    def add_forcing(self, x, t, dxdt,coef=None):
        """ add the forcing term on x[0]=the vorticity """
        dxdt[0] += self.forc * np.sin(2*np.pi*t/self.period)

        # attempt to have a forcing white noise in time
        # coef = np.random.uniform()*2-1
        # dxdt[0] += self.forc * coef
