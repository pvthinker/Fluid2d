from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'taylorinst_00'

# domain and resolution
param.nx = 64*2

param.ny = param.nx
param.npx = 1
param.npy = 1
param.Lx = 1.
param.Ly = param.Lx
param.geometry = 'xchannel'

# time
param.tend = 10.
param.cfl = 1.5
param.adaptable_dt = True
param.dt = 0.004
param.dtmax = 1e-2

# discretization
param.order = 5

# output
param.plot_var = 'buoyancy'
param.var_to_save = ['vorticity', 'buoyancy', 'psi']
param.list_diag = 'all'
param.freq_his = 1.
param.freq_diag = 1.

# plot
param.plot_interactive = True
param.freq_plot = 10
param.colorscheme = 'imposed'
param.cax = [-.6, .6]
param.generate_mp4 = True

# physics
param.gravity = 1.
param.forcing = False
param.diffusion = False
param.noslip = False

grid = Grid(param)

xr, yr = grid.xr, grid.yr
Lx, Ly = param.Lx, param.Ly

# to set a wall with a small aperture, uncomment
# grid.msk[((xr < 0.4*Lx) | (xr > 0.6*Lx)) &
#          (yr > 0.5-grid.dx) & (yr < 0.5+grid.dx)] = 0
# grid.msknoslip = grid.msk.copy()
# grid.finalize_msk()

param.Kdiff = 2e-3*grid.dx


f2d = Fluid2d(param, grid)
model = f2d.model

xr, yr = grid.xr, grid.yr
buoy = model.var.get('buoyancy')


def sigmoid(x, delta):
    return 1/(1+np.exp(-(x-0.5)/delta))


def stratif():
    sigma = grid.dx/2  # width of the interface
    b = sigmoid(yr/param.Ly, sigma/param.Lx)
    return b * grid.msk


buoy[:] = (1-stratif() - 0.5)*grid.msk
# add noise to trigger the instability
noise = np.random.normal(size=np.shape(yr), scale=1.)*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area
grid.fill_halo(noise)

buoy += 1e-3*noise

model.set_psi_from_vorticity()


f2d.loop()
