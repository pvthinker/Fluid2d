from param import Param
from grid import Grid
from fluid2d import Fluid2d
from numpy import exp, sqrt, pi, cos, where, random, shape, tanh, cumsum, cosh
from restart import Restart
from island import Island

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'LW_0'

# domain and resolution
ratio = 2
param.ny = 64*2
param.nx = param.ny*ratio
param.Ly = 1.
param.Lx = param.Ly*ratio
param.npx = 1
param.npy = 1
param.geometry = 'xchannel'

# time
param.tend = 10.
param.cfl = 0.9
param.adaptable_dt = True
param.dt = 1e-2
param.dtmax = 1.

# discretization
param.order = 5
param.timestepping = 'RK3_SSP'

# output
param.var_to_save = ['vorticity', 'psi', 'u', 'v', 'buoyancy', 'banom']
param.list_diag = ['ke', 'pe', 'vorticity', 'enstrophy']
param.freq_plot = 20
param.freq_his = .5
param.freq_diag = .1

# plot
param.plot_interactive = True
param.plot_var = 'v'
param.cax = [-.25, .25]
param.colorscheme = 'imposed'
param.generate_mp4 = True
param.plot_psi = True

# physics
param.forcing = False
param.noslip = False
param.diffusion = False
param.forcing = False
param.isisland = True

param.gravity = 1.

nh = param.nh


grid = Grid(param)

# control parameters of the experiment
N2 = 10.  # Brunt Vaisala frequency squared
hseamount = 0.1
Lseamount = 0.1
U = .5

# add a seamount
z = grid.yr-hseamount*exp(-(grid.xr0+0.5)**2/(Lseamount**2/2))
grid.msk[z <= 0] = 0


if grid.j0 == param.npy-1:
    msk = grid.msk.copy()*0
    msk[-nh:-1, :] = 1
    idx = where(msk == 1)
    grid.island.add(idx, -U)


grid.finalize_msk()


f2d = Fluid2d(param, grid)
model = f2d.model

xr, yr = grid.xr, grid.yr

buoy = model.var.get('buoyancy')


# linear stratification
buoy[:, :] = grid.yr0*grid.msk*N2
model.bref[:, :] = buoy


model.set_psi_from_vorticity()


f2d.loop()
