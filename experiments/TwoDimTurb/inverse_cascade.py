from fluid2d import Fluid2d
from param import Param
from grid import Grid
import numpy as np


class Forcing:
    """ define the forcing """

    def __init__(self, param, grid):
        self.list_param = ['sizevar']
        param.copy(self, self.list_param)

        self.list_param = ['nh', 'msk', 'fill_halo']
        grid.copy(self, self.list_param)

        self.intensity = 1e-1
        self.gamma = 1e-1
        self.t = 0

        self.forc = np.random.normal(size=self.sizevar) * self.msk
        self.fill_halo(self.forc)

    def add_forcing(self, x, t, dxdt):
        """ add the forcing term on x[0]=the vorticity """

        dt = t - self.t
        self.t = t

        # Define the forcing as the sum of the forcing in the last step and
        # a new random forcing
        self.forc = (
            (1 - dt * self.gamma) * self.forc
            + dt * self.gamma * np.random.normal(size=self.sizevar) * self.msk
        )
        self.fill_halo(self.forc)

        dxdt[0] += self.intensity * self.forc


param = Param('default.xml')
param.modelname = 'euler'
param.expname = 'twodimturb_00'

# domain and resolution
param.nx = 64*2
param.ny = param.nx
param.Ly = param.Lx
param.npx = 1
param.npy = 1
param.geometry = 'perio'

# time
param.tend = 100.
param.cfl = 1.5
param.adaptable_dt = True
param.dt = 1.
param.dtmax = .5

# discretization
param.order = 5

# output
param.var_to_save = ['vorticity', 'psi', 'tracer']
param.list_diag = 'all'
param.freq_plot = 5
param.freq_his = 10
param.freq_diag = 1

# plot
param.plot_interactive = True
param.plot_var = 'vorticity'
param.plot_psi = True
param.cax = np.array([-1, 1])*.5
param.colorscheme = 'imposed'
param.generate_mp4 = True

# physics
# you may activate the forcing defined above
# it is a white noise forcing with some time correlation
# such forcing is often used in turbulence studies
param.forcing = True
param.decay = False  # set it to False if forcing == True
param.noslip = False
param.diffusion = False

# add a passive tracer
param.additional_tracer = ['tracer']

grid = Grid(param)
param.Kdiff = 5e-4*grid.dx

f2d = Fluid2d(param, grid)
model = f2d.model

# set the forcing
model.forc = Forcing(param, grid)

xr, yr = grid.xr, grid.yr
vor = model.var.get('vorticity')
trac = model.var.get('tracer')

# set an initial small scale random vorticity field (white noise)
np.random.seed(1)
noise = np.random.normal(size=np.shape(yr))*grid.msk
# noise=model.ope.gmg.generate_random_large_scale_field()

grid.fill_halo(noise)
noise -= grid.domain_integration(noise)*grid.msk/grid.area
vor[:] = noise

# initialize the passive tracer with square tiles (better idea?)
trac[:] = np.round(xr*6) % 2 + np.round(yr*6) % 2


model.set_psi_from_vorticity()
f2d.loop()
