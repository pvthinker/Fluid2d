from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'RB_00'

# domain and resolution
param.nx = 64*2
param.ny = param.nx/4
param.npx = 1
param.Lx = 4.
param.Ly = 1.
param.geometry = 'xchannel'

# time
param.tend = 40.
param.cfl = 1.
param.adaptable_dt = True
param.dt = .1
param.dtmax = .1
param.exacthistime = False

# discretization
param.order = 5
param.aparab = 0.02

# output
param.plot_var = 'buoyancy'
param.var_to_save = ['vorticity', 'buoyancy', 'v', 'psi']
param.list_diag = 'all'
param.freq_his = 1
param.freq_diag = 1

# plot
param.plot_interactive = True
param.freq_plot = 10
param.colorscheme = 'imposed'
param.cax = [-.5, .5]
param.generate_mp4 = True

# physics
param.gravity = 1.
param.forcing = True
param.forcing_module = 'forcing_rayleigh'
param.diffusion = True
param.noslip = True

grid = Grid(param)
# Prandtl number is Kvorticity / Kbuoyancy

prandtl = 1.


deltab = 600  # this is directly the Rayleigh number is L=visco=diffus=1

param.deltab = deltab  # make it known to param

L = param.Ly
visco = .002*grid.dy
diffus = visco / prandtl

# Rayleigh number is
Ra = deltab * L ** 3 / (visco * diffus)
print('Rayleigh number is %4.1f' % Ra)

param.Kdiff = {}
param.Kdiff['vorticity'] = visco
param.Kdiff['buoyancy'] = diffus  # param.Kdiff['vorticity'] / prandtl

# time step is imposed by the diffusivity
# param.dt = 0.25*grid.dx**2 / max(visco, diffus)
# param.dtmax = param.dt
# print('dt = %f' % param.dt)

f2d = Fluid2d(param, grid)
model = f2d.model

xr, yr = grid.xr, grid.yr

buoy = model.var.get('buoyancy')
# add noise to trigger the instability
noise = np.random.normal(size=np.shape(yr))*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area
grid.fill_halo(noise)

buoy += 1e-1*noise

model.set_psi_from_vorticity()

f2d.loop()
