from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'quasigeostrophic'
param.expname = 'turb_qg_0'

# domain and resolution
param.nx = 64*2
param.ny = 64
param.npy = 1
param.Lx = 2.
param.Ly = param.Lx/2
param.geometry = 'xchannel'

# time
param.tend = 5000.
param.cfl = 0.6
param.adaptable_dt = True
param.dt = 1.
param.dtmax = 10.

# discretization
param.order = 5

# output
param.var_to_save = ['pv', 'psi', 'u', 'v']
param.list_diag = ['ke', 'pv', 'pv2']
param.freq_plot = 5
param.freq_his = 10
param.freq_diag = 1

# plot
param.plot_interactive = True
param.plot_psi = True
param.cmap = 'inferno'
param.plot_var = 'pv'
param.cax = np.array([-1, 1])*.5
param.colorscheme = 'imposed'

# physics
param.beta = 1.
param.Rd = 50.  # 2*grid.dx
param.forcing = False
param.noslip = False
param.diffusion = False
#
# you may activate the forcing and use the forcing below
# it's a white noise forcing with some time correlation
# such forcing is often used in turbulence studies
# param.forcing_module='forcing_euler'

grid = Grid(param)
param.Kdiff = 5e-4*grid.dx

f2d = Fluid2d(param, grid)
model = f2d.model


xr, yr = grid.xr, grid.yr
pv = model.var.get('pv')

amplitude = 1.

# set an initial small scale random vorticity field (white noise)
noise = np.random.normal(size=np.shape(yr))*grid.msk
# noise=model.ope.gmg.generate_random_large_scale_field()

grid.fill_halo(noise)
noise -= grid.domain_integration(noise)*grid.msk/grid.area
pv[:] = noise*amplitude

sigma = 2*grid.dx
Bu = (param.Rd / sigma)**2
U = min(1., Bu)*np.abs(amplitude)*sigma
Lbeta = np.sqrt(U/param.beta)
lamb = np.abs(amplitude)/(param.beta*sigma)
Rh = (Lbeta/sigma)**2

print('='*40)
print('Burger Bu = %f' % Bu)
print('Rhines Rh = %f' % Rh)
print('Lambda La = %f' % lamb)
print('U = %f' % U)
print('Lbeta = %f' % Lbeta)
print('='*40)

model.add_backgroundpv()

model.set_psi_from_pv()
f2d.loop()
