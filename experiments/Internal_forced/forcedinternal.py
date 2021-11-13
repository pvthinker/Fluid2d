from param import Param
from grid import Grid
from fluid2d import Fluid2d
from numpy import exp, sqrt, pi, cos, sin, where, random, shape, tanh, cumsum, cosh
import numpy as np
from restart import Restart
from island import Island

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'ForcedInternal_0'

# domain and resolution
ratio = 1
param.ny = 64*2
param.nx = param.ny*ratio
param.Ly = 1.
param.Lx = param.Ly*ratio
param.npx = 1
param.npy = 1
param.geometry = 'closed'

param.hydroepsilon = 1.

# time
param.tend = 50.
param.cfl = 0.1
param.adaptable_dt = False
param.dt = 1e-2
param.dtmax = .1

# discretization
param.order = 5
param.timestepping = 'RK3_SSP'  # 'LFAM3'

# output
param.var_to_save = ['vorticity', 'psi', 'u', 'v', 'buoyancy', 'tracer']
param.list_diag = ['ke', 'pe', 'energy', 'brms', 'vorticity', 'enstrophy']
param.freq_plot = 20
param.freq_his = .2
param.freq_diag = .1

# plot
param.plot_interactive = True
param.plot_var = 'banom'
param.cax = np.asarray([-0.1, 0.1])*1e-1
param.colorscheme = 'imposed'

# physics
param.noslip = False
param.diffusion = False
param.forcing = True
param.forcing_module = 'forcing'
param.additional_tracer = ['tracer']
param.generate_mp4 = True

param.gravity = 1.

nh = param.nh

grid = Grid(param)
xr, yr = grid.xr, grid.yr

# add a topography (mask)
# hb = exp(-(xr/param.Lx-0.2)**2 * 50)*0.35
# grid.msk[(yr < hb)] = 0
# grid.finalize_msk()


f2d = Fluid2d(param, grid)
model = f2d.model

xr, yr = grid.xr, grid.yr

buoy = model.var.get('buoyancy')

# control parameters of the experiment
N2 = 1.  # Brunt Vaisala frequency squared


# stratification

# continuous
buoy[:, :] = grid.yr * grid.msk * N2

# exponentially decreasing stratification (promoting breaking)
y = grid.yr/param.Ly
# buoy[:, :] = np.exp(y*4)/np.exp(4)*param.Ly * grid.msk * N2

nlayers = 30
# layers
#buoy[:, :] = np.floor(grid.yr*nlayers/param.Ly)*(param.Ly/nlayers)*N2*grid.msk

model.bref[:, :] = buoy


model.set_psi_from_vorticity()

# 2/ set an initial tracer field


def vortex(x0, y0, sigma):
    x = sqrt((xr-param.Lx*x0)**2+(yr-param.Ly*y0)**2)
    y = 2-2./(1+exp(-x/sigma))
    return y


sigma = 0.01*param.Lx

state = model.var.get('tracer')

tracer_config = 1

if tracer_config == 0:  # single localized patch
    state[:, :] = 2*vortex(0.3, 0.4, sigma)

if tracer_config == 1:  # tiles
    state[:, :] = np.round(xr*6) % 2 + np.round(yr*6) % 2

if tracer_config == 2:  # a closed line
    x = xr/param.Lx
    y = yr/param.Ly
    msk = np.zeros_like(yr)
    msk[(x > 0.4) & (x < 0.7) & (y > 0.4) & (y < 0.7)] = 1
    z = np.roll(msk, -1, axis=1)+np.roll(msk, -1, axis=0) + \
        np.roll(msk, +1, axis=1)+np.roll(msk, +1, axis=0)-4*msk
    state[:, :] = 0.
    state[z > 0] = 4.


state[:, :] *= grid.msk


f2d.loop()
