from fluid2d import Fluid2d
from param import Param
from grid import Grid
import numpy as np


class Forcing:
    """ define the forcing """

    def __init__(self, param, grid):
        yr = (grid.yr - param.Ly/2) / param.Ly

        # jet width relative to param.Ly
        sigma = 0.1

        # forcing intensity
        tau0 = 1e-4

        # type of forcing: comment/uncomment choice A vs. B

        # choice A/ jet-like forcing (localized in y)
        # self.forc = tau0 * (yr/sigma)*np.exp(-yr**2/(2*sigma**2)) * grid.msk

        # choice B/ basin-scale forcing: double gyre configuration
        self.forc = tau0 * np.sin(yr * np.pi) * grid.msk

        total = grid.domain_integration(self.forc)
        self.forc -= (total / grid.area) * grid.msk

    def add_forcing(self, x, t, dxdt):
        """ add the forcing term on x[0]=the vorticity """
        dxdt[0] += self.forc


param = Param('default.xml')
param.modelname = 'quasigeostrophic'
param.expname = 'dbl_gyre_00'

# domain and resolution
param.nx = 64*2
param.ny = 64
param.npy = 1
param.Lx = 2.
param.Ly = param.Lx/2
param.geometry = 'closed'

# time
param.tend = 2000.
param.cfl = 1.5
param.adaptable_dt = True
param.dt = 1.
param.dtmax = 100.

# discretization
param.order = 5


# output
param.var_to_save = ['pv', 'psi', 'pvanom']
param.list_diag = 'all'
param.freq_his = 50
param.freq_diag = 10

# plot
param.plot_var = 'pv'
param.freq_plot = 10
a = 0.5
param.cax = [-a, a]
param.plot_interactive = True
param.colorscheme = 'imposed'
param.generate_mp4 = True

# physics
param.beta = 1.
param.Rd = 1.  # 2*grid.dx
param.forcing = True
param.noslip = False
param.diffusion = False

grid = Grid(param)
param.Kdiff = 0.5e-4*grid.dx

# add an island
# grid.msk[28:32,34:38]=0
# grid.finalize_msk()


f2d = Fluid2d(param, grid)
model = f2d.model

# set the forcing
model.forc = Forcing(param, grid)

xr, yr = grid.xr, grid.yr
vor = model.var.get('pv')


# set an initial tracer field

vor[:] = 1e-2*np.random.normal(size=np.shape(vor))*grid.msk
y = vor[:]*1.
model.ope.fill_halo(y)
vor[:] = y

model.add_backgroundpv()

model.set_psi_from_pv()

f2d.loop()
