"""
Gravity current along a topographic slope
"""
from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'gravcurr_00'

# domain and resolution

param.ny = 128
param.nx = param.ny*2
param.npx = 1
param.npy = 1
param.Lx = 2.
param.Ly = 1.
param.geometry = 'closed'

# time
param.tend = 5.
param.cfl = 1.2
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
param.cax = [0, 3]
param.generate_mp4 = True

# physics
param.gravity = 1.
param.forcing = False
param.diffusion = True
param.noslip = False

grid = Grid(param)
param.Kdiff = 2e-3*grid.dx**2

xr, yr = grid.xr, grid.yr
# add a mask
hb = np.exp(-(xr/param.Lx-0.5)**2 * 50)*0.7
grid.msk[(yr < hb)] = 0
grid.finalize_msk()

f2d = Fluid2d(param, grid)
model = f2d.model


buoy = model.var.get('buoyancy')


def sigmoid(x, delta):
    return 1/(1+np.exp(-(x-0.5)/delta))


def stratif():
    sigma = 3*grid.dx  # width of the interface
    b = sigmoid(xr/param.Lx, sigma/param.Lx)
    return (1-b) * grid.msk


# to have a look, exchange the type of flow
buoy[:] = (1-stratif() - 0.5)*grid.msk

buoy[:] = (yr * (1 + 1*(1+np.tanh((xr-param.Lx/2)/0.1)))) * grid.msk

# add noise to trigger the instability
noise = np.random.normal(size=np.shape(yr), scale=1.)*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area
grid.fill_halo(noise)

buoy += 1e-8*noise

model.set_psi_from_vorticity()


f2d.loop()
