from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'khi_0'

# domain and resolution
ratio = 2
param.ny = 64*2
param.nx = param.ny*ratio
param.Ly = 1.
param.Lx = 1*ratio
param.npx = 1
param.npy = 1
param.geometry = 'xchannel'

# time
param.tend = 5.
param.cfl = 1.5
param.adaptable_dt = True
param.dt = 0.01
param.dtmax = 0.02

# discretization
param.order = 5
param.timestepping = 'RK3_SSP'

# output
param.var_to_save = ['vorticity', 'psi', 'u', 'v', 'buoyancy', 'banom']
param.list_diag = 'all'
param.freq_plot = 5
param.freq_his = .25
param.freq_diag = .1

# plot
param.plot_interactive = True
param.plot_var = 'buoyancy'
param.cax = [-8., 8.]
param.colorscheme = 'imposed'
param.generate_mp4 = True
param.cmap = 'inferno'

# physics
param.forcing = False
param.noslip = False
param.diffusion = False
param.forcing = False

param.gravity = 1.

nh = param.nh

grid = Grid(param)


f2d = Fluid2d(param, grid)
model = f2d.model

xr, yr = grid.xr, grid.yr
vor = model.var.get('vorticity')
buoy = model.var.get('buoyancy')

# control parameters of the experiment
N = 4.  # Brunt Vaisala frequency squared
S = 20.  # vertical shear

# linear stratification
buoy[:, :] = grid.yr0*grid.msk*N**2
model.bref[:] = buoy

print('Ri = %4.2f' % (N**2/S**2))

# we take the mean jet profile as a tanh
# U = tanh( grid.yr0 * sqrt(S2) )
# but we control it via the vorticity distribution

vor[:, :] = -S / np.cosh(grid.yr0 * S)**2.*grid.msk


# add noise to trigger the instability
np.random.seed(42)  # set the seed to ensure reproducibility
noise = np.random.normal(size=np.shape(yr))*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area
grid.fill_halo(noise)


vor[:, :] += 1e0*noise

vor *= grid.msk


model.set_psi_from_vorticity()


f2d.loop()
