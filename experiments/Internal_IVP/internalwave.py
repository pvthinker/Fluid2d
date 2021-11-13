from param import Param
from grid import Grid
from fluid2d import Fluid2d
from numpy import exp, sqrt, pi, cos, sin, where, random, shape, tanh, cumsum, cosh
import numpy as np
from restart import Restart
from island import Island

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'Internal_0'

# domain and resolution
ratio = 2
param.ny = 64*2
param.nx = param.ny*ratio
param.Ly = 1.
param.Lx = param.Ly*ratio
param.npx = 1
param.npy = 1
param.geometry = 'closed'

param.hydroepsilon = 1

# time
param.tend = 5.
param.cfl = 0.2
param.adaptable_dt = False
param.dt = 1e-2
param.dtmax = 1.

# discretization
param.order = 5
param.timestepping = 'RK3_SSP'  # 'LFAM3'

# output
param.var_to_save = ['vorticity', 'psi', 'u', 'v', 'buoyancy', 'banom']
param.list_diag = ['ke', 'pe', 'energy', 'brms', 'vorticity', 'enstrophy']
param.freq_plot = 5
param.freq_his = .2
param.freq_diag = .1

# plot
param.plot_interactive = True
param.plot_var = 'banom'
param.cax = np.asarray([-0.1, 0.1])
param.colorscheme = 'imposed'  # 'fixed'
param.generate_mp4 = True

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
N2 = 10.  # Brunt Vaisala frequency squared
nlayers = 10
background_stratif = 'continuous'
delta = param.Lx*0.03
delta = grid.dx*4
amplitude = .2

if background_stratif == 'continuous':
    buoy[:, :] = grid.yr * grid.msk * N2

elif background_stratif == 'layered':

    y = grid.yr*nlayers/param.Ly
    buoy[:, :] = np.floor(y) * (param.Ly/nlayers)*N2*grid.msk
else:
    print('background stratification is undefined')
    exit(0)
model.bref[:, :] = buoy

y = -grid.yr0
perturb = amplitude*exp(-(grid.yr0**2+grid.xr0**2)/(2*delta**2))
buoy += perturb

# we perturb the fluid by introducing a localized dipole of vorticity
x = grid.xr0
#vor[:, :] = perturb*sin(x*2*pi*2)/delta * grid.msk


model.set_psi_from_vorticity()


f2d.loop()
