from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

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
param.forcing = True
param.decay = False  # set it to False if forcing == True
param.noslip = False
param.diffusion = False
#
# you may activate the forcing and use the forcing below
# it's a white noise forcing with some time correlation
# such forcing is often used in turbulence studies
param.forcing_module = 'forcing_euler'

# add a passive tracer
param.additional_tracer = ['tracer']

grid = Grid(param)
param.Kdiff = 5e-4*grid.dx

f2d = Fluid2d(param, grid)
model = f2d.model


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
