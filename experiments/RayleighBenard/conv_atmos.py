"""
Convective boundary layer

ingredients:

   - linear reference buoyancy profile (stable)
   - uniform heating from the bottom
   - no cooling in the column

   - viscosity and diffusion (taken equal)
"""
from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'RB_128'

# domain and resolution
param.nx = 128*2
param.ny = 128
param.npx = 1
param.Lx = 2.
param.Ly = 1.
param.geometry = 'xchannel'

# time
param.tend = 25.
param.cfl = 0.9 # <- 1.0 is too large

assert param.cfl <= 0.9

param.adaptable_dt = True
param.dt = 1e2
param.dtmax = 0.2
param.exacthistime = True

# discretization
param.order = 5

# output
param.plot_var = 'buoyancy'
param.var_to_save = ['vorticity', 'buoyancy', 'v', 'psi']
param.list_diag = 'all'
param.freq_his = 0.2
param.freq_diag = 0.2

# plot
param.plot_interactive = True
param.freq_plot = 10
param.colorscheme = 'imposed'
param.cax = [0, 2]
param.generate_mp4 = False

# physics
param.gravity = 1.
param.forcing = True
param.forcing_module = 'embedded'
param.diffusion = True
param.noslip = True

grid = Grid(param)
# Prandtl number is Kvorticity / Kbuoyancy

prandtl = 1.

visco = .05*grid.dy
diffus = visco / prandtl


param.Kdiff = {}
param.Kdiff['vorticity'] = visco
param.Kdiff['buoyancy'] = diffus

param.heatflux = 1e-2

class Forcing:
    """ define the forcing """

    def __init__(self, param, grid):

        self.list_param = ['deltab', 'dt']
        param.copy(self, self.list_param)

        self.list_param = ['j0', 'npy', 'nh']
        grid.copy(self, self.list_param)


        Q = param.heatflux
        nh = param.nh
        self.nh = nh
        self.forc = np.zeros_like(grid.yr)

        if grid.j0 == 0:
            self.forc[nh, :] = +Q

        self.forc *= grid.msk

        # transform the surface flux into a volume flux
        self.forc *= (1./grid.dx)
        self.dz = grid.dy
        self.K = param.Kdiff['buoyancy']
        
    def add_forcing(self, x, t, dxdt,coef=1):
        """ add the forcing term on x[0]=the vorticity """
        dxdt[4] += self.forc*coef
        nh=self.nh
        dxdt[4][-nh-1] += self.K*(x[4][-nh-1]-x[4][-nh-2])/self.dz**2*coef


f2d = Fluid2d(param, grid)
model = f2d.model

model.forc = Forcing(param, grid)

xr, zr = grid.xr, grid.yr

buoy = model.var.get('buoyancy')
# linear reference stratification
buoy[:] = zr*2

# add noise to trigger the instability
np.random.seed(42)
noise = np.random.normal(size=np.shape(xr))*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area
grid.fill_halo(noise)

buoy += 1e-3*noise

model.set_psi_from_vorticity()

f2d.loop()
#print(buoy.shape, nh)
