"""
Lock exchange experiment
"""
from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'lockexch'

# domain and resolution

param.ny = 64
param.nx = param.ny*8
param.npx = 1
param.npy = 1
param.Lx = 8.
param.Ly = 1.
param.geometry = 'closed'

# time
param.tend = 8.
param.cfl = 1.
param.adaptable_dt = True
param.dt = .1
param.dtmax = .1

# discretization
param.order = 5


# output
param.plot_var = 'buoyancy'
param.var_to_save = ['vorticity', 'buoyancy', 'psi']
param.list_diag = 'all'
param.freq_his = .2
param.freq_diag = .1

# plot
param.plot_interactive = True
param.freq_plot = 10
param.colorscheme = 'imposed'
param.cax = [-1, 1]
param.generate_mp4 = False

# physics
param.gravity = 1.
param.forcing = False
param.diffusion = False
param.noslip = False
# by default hydroepsilon == 1 -> non-hydrostatic case
# set it to value << 1 to switch toward a more hydrostatic case
#param.hydroepsilon = 0.01

grid = Grid(param)
param.Kdiff = 2e-3*grid.dx**2

xr, yr = grid.xr, grid.yr
grid.finalize_msk()

f2d = Fluid2d(param, grid)
model = f2d.model


buoy = model.var.get('buoyancy')


x0 = param.Lx*0.5
width = 3*grid.dx
buoy[:] = -np.tanh((xr-x0)/width) * grid.msk


f2d.loop()
