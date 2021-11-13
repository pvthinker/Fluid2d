# footprint of an advection scheme (space+time)

from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np
import matplotlib.pyplot as plt


param = Param('default')
param.modelname = 'advection'
param.expname = 'footprint_0'
param.nx = 32*2
param.ny = 32*2
param.npy = 1
param.Lx = 1.
param.Ly = 1.
param.geometry = 'disc'
param.tend = 40.
param.cfl = 1.0
param.order = 5
param.timestepping = 'RK3_SSP'

param.plot_var = 'tracer'
param.var_to_save = ['tracer']
param.list_diag = ['mean', 'rms']
param.freq_his = 10
param.freq_diag = 1
param.freq_plot = 1
param.colorscheme = 'minmax'
param.cax = [-1.2, 1.2]
param.plot_interactive = True
param.adaptable_dt = True
param.dt = .1
param.dtmax = 1e4
param.diffusion = True
param.Kdiff = 0.

grid = Grid(param)


# grid.msk[-60:-40,40:60]=0
# grid.msk[:64,64]=0

f2d = Fluid2d(param, grid)
model = f2d.model

xr, yr = grid.xr, grid.yr
# 1/ set the stream function
sigma = 2*param.Lx
state = model.var.get('tracer')
state[:] = 1.  # body rotation
state[:] = state[:]*grid.msk
model.set_psi_from_tracer()


nh = f2d.nh

# 2/ set one grid cell to 1
state[:] = 0.
state[nh+param.ny//3, nh+param.nx//3] = 1.
state[:] *= grid.msk


# f2d.loop()


t = 0.
var = f2d.model.var

# advect it for one time step only

f2d.model.diagnostics(var, 0.)
f2d.set_dt(0)
for k in range(1):
    f2d.model.step(t, f2d.dt)

print('everything was ok')

# look at the result

z2d = f2d.model.var.get(param.plot_var)[nh:-nh, nh:-nh]
msk = grid.msk[nh:-nh, nh:-nh]
z2d = np.log10(1e-6+abs(z2d))

f2d.plotting.set_cax(z2d)
cax = f2d.plotting.cax
z2d[msk == 0] = np.nan
fig = plt.figure()
im = plt.imshow(z2d, vmin=cax[0], vmax=cax[1], origin='lower')
plt.title('combined space-time stencil (in log-scale)')
plt.colorbar(im)
plt.show()
