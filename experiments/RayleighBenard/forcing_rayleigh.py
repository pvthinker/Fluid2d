from param import Param


class Forcing:
    """ define the forcing """

    def __init__(self, param, grid):

        self.list_param = ['deltab', 'dt']
        param.copy(self, self.list_param)

        self.list_param = ['j0', 'npy', 'nh']
        grid.copy(self, self.list_param)

        yr = (grid.yr-param.Ly*.5)/param.Ly

        Q = 1e-2
        nh = param.nh
        self.forc = yr*0.

        type_forcing = 'coolroof'
        #type_forcing = 'horizontalconvection'

        # type_forcing = 'prescribed_buoyancy'

        if type_forcing == 'coolroof':
            if grid.j0 == grid.npy-1:
                self.forc[-nh-1, :] = -Q
            if grid.j0 == 0:
                self.forc[nh, :] = +Q

        elif type_forcing == 'spreadcooling':
            # bottom heating
            if grid.j0 == 0:
                self.forc[nh, :] = +Q
            # uniform cooling
            self.forc -= Q/param.ny

        elif type_forcing == 'horizontalconvection':
            if grid.j0 == grid.npy-1:
                xr = grid.xr
                self.forc[xr < param.Lx/2] = Q
                self.forc[xr > param.Lx/2] = -Q
                self.forc[:-nh-1, :] = 0.

        self.type_forcing = type_forcing

        self.forc *= grid.msk

        # transform the surface flux into a volume flux
        self.forc *= (1./grid.dx)

    def add_forcing(self, x, t, dxdt, coef=1.):
        """ add the forcing term on x[0]=the vorticity """
        if self.type_forcing == 'prescribed_buoyancy':
            #
            # for this experiment you need to use the fixed time step
            # whose value is controlled by max(visco,diffus)
            #
            if self.j0 == self.npy-1:
                dxdt[4][-self.nh-1, :] -= (x[4][-self.nh-1, :]
                                           + self.deltab*.5)/self.dt
            if self.j0 == 0:
                dxdt[4][self.nh, :] -= (x[4][self.nh, :]
                                        - self.deltab*.5)/self.dt

        else:
            dxdt[4] += self.forc
        dxdt[4] *= coef
