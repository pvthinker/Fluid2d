"""
Numerical experiment of the gravity current done in the lab
at FDSE, with Paul Billant
"""
from fluid2d import Fluid2d
from param import Param
from grid import Grid
import numpy as np

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'gravcurr'

# domain and resolution

scaling = 1.

param.ny = 64
param.nx = 8*param.ny
param.npx = 1
param.npy = 1
param.Lx = 8.
param.Ly = 1.
param.geometry = 'closed'

# time
param.tend = 15./scaling
param.cfl = .5
param.adaptable_dt = True
param.dt = .1
param.dtmax = .1

# discretization
param.order = 5


# output
param.plot_var = 'buoyancy'
param.var_to_save = ['vorticity', 'buoyancy', 'psi']
param.list_diag = 'all'
param.freq_his = .1/scaling
param.freq_diag = .1/scaling

# plot
param.plot_interactive = False
param.freq_plot = 10
param.colorscheme = 'imposed'
param.cax = [0, 3]
param.generate_mp4 = True

# physics
param.gravity = 1.
param.forcing = False
param.diffusion = False
param.noslip = False

grid = Grid(param)
param.Kdiff = 2e-3*grid.dx**2

xr, yr = grid.xr, grid.yr
# add a mask
#hb = np.exp(-(xr/param.Lx-0.5)**2 * 50)*0.1
#grid.msk[(yr < hb)] = 0
grid.finalize_msk()

f2d = Fluid2d(param, grid)
model = f2d.model


rho = model.var.get('buoyancy')


# #
# sigma = 3*grid.dx
# buoy[:] = -0.5*(np.tanh( (xr-param.Lx*0.9)/sigma)+1)
# buoy *= grid.msk

dx = grid.dx
rho[:] = 0.5-0.5*np.tanh((xr-0.85*param.Lx)/(dx/2))
ycut = np.tanh((yr-param.Lx*0.1)/(dx/2))
rho[ycut>0] = 1
#rho *= 0.25

rho *= scaling**2
#model.set_psi_from_vorticity()


f2d.loop()
