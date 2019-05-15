from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

""" Reentrant channel with prescribed transport through the channel. The transport is controlled via psi0, the streamfunction at the North Wall -- psi=0 on the South Wall. An island with psi=psi0/2 is also set in the middle of the channel. Even though there is no explicit forcing, this setup creates an island wake. Parameter controls include:

- Rd/L, relative size of the deformation radius
- beta Rd^2 L / psi0, relative importance of the jet speed vs Rossby wave speed
""" 

param = Param('default.xml')
param.modelname = 'quasigeostrophic'
param.expname = 'channel'

# domain and resolution
param.nx = 64*2
param.ny = 64
param.npy = 1
param.Lx = 2.
param.Ly = param.Lx/2
param.geometry = 'xperio'

# time
param.tend = 200.
param.cfl = 1.2
param.adaptable_dt = True
param.dt = 1.
param.dtmax = 100.

# discretization
param.order = 5


# output
param.var_to_save = ['pv', 'u', 'v', 'psi', 'pvanom']
param.list_diag = 'all'
param.freq_his = 10
param.freq_diag = 5.

# plot
param.plot_var = 'pv'
param.freq_plot = 10
a = 0.5
param.cax = [-a, a]
param.plot_interactive = True
param.colorscheme = 'imposed'
param.generate_mp4 = False

# physics
param.beta = 1.
param.Rd = 10*param.Lx/param.nx
param.forcing = False
param.forcing_module = 'forcing_dbl_gyre'
param.noslip = False
param.diffusion = False
param.isisland=True

psi0 = -5e-4 # this sets psi on the Northern wall (psi=0 on the Southern wall)


grid = Grid(param)

nh = grid.nh

def disc(param, grid, x0, y0, sigma):
    """ function to set up a circular island
    note that this function returns the indices, not the mask """
    xr, yr = grid.xr, grid.yr
    r = np.sqrt( (xr-param.Lx*x0)**2+(yr-param.Ly*y0)**2 )
    return np.where(r<=sigma)

sigma = 0.08
idx = disc(param, grid, 0.125*param.Lx, 0.5, sigma)
grid.msk[idx] = 0
grid.msknoslip[idx] = 0
grid.island.add(idx, psi0*.5)


if grid.j0 == param.npy-1:
    msk = grid.msk.copy()*0
    msk[-nh:-1, :] = 1
    idx = np.where(msk == 1)
    grid.island.add(idx, psi0)

# tell the code to deactivate the no-slip along the outer walls
if grid.j0 == 0:
    grid.msknoslip[:nh, :] = 1
if grid.j0 == param.npy-1:
    grid.msknoslip[-nh:, :] = 1

param.Kdiff = 0.5e-4*grid.dx

f2d = Fluid2d(param, grid)
model = f2d.model

xr, yr = grid.xr, grid.yr
vor = model.var.get('pv')



# set an initial tracer field

vor[:] = 1e-4*np.random.normal(size=np.shape(vor))*grid.msk
y = vor[:]*1.
model.ope.fill_halo(y)
vor[:] = y

model.add_backgroundpv()

model.set_psi_from_pv()

f2d.loop()
