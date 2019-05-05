from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'quasigeostrophic'
param.expname = 'dbl_gyre_00'

# domain and resolution
param.nx = 64*2
param.ny = 64
param.npy = 1
param.Lx = 2.
param.Ly = param.Lx/2
param.geometry = 'square'

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
param.forcing_module = 'forcing_dbl_gyre'
param.noslip = False
param.diffusion = False

grid = Grid(param)
param.Kdiff = 0.5e-4*grid.dx

# add an island
# grid.msk[28:32,34:38]=0
# grid.finalize_msk()


f2d = Fluid2d(param, grid)
model = f2d.model

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
