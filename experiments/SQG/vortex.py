from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param("toto")
param.modelname = 'sqg'
param.expname = 'sqg_vortex'

# domain and resolution
param.nx = 64*2
param.ny = param.nx
param.npy = 1
param.Lx = 1.
param.Ly = 1.
param.geometry = 'perio'

# time
param.tend = 40
param.cfl = 0.8
param.adaptable_dt = True
param.dt = 1.
param.dtmax = 100.

# discretization
param.order = 5

# output
param.ageostrophic = False
param.var_to_save = ['pv', 'psi', 'u', 'v', 'vorticity']
param.list_diag = ['pv', 'pv2']
param.freq_his = 1
param.freq_diag = 1

# plot
param.plot_var = 'pv'
param.freq_plot = 10
a = 0.2
param.cax = [-a, a]
param.plot_interactive = False
param.colorscheme = 'symmetric'
param.generate_mp4 = False

# physics
param.beta = 0.
param.Rd = 10.  # 2*grid.dx
param.forcing = False
param.noslip = False
param.diffusion = False

grid = Grid(param)
param.Kdiff = 0.1e-3*grid.dx



f2d = Fluid2d(param, grid)
model = f2d.model

xr, yr = grid.xr, grid.yr
vor = model.var.get('pv')


def vortex(param, grid,
           x0, y0, sigma,
           vortex_type,
           ratio=1):

    xr, yr = grid.xr, grid.yr
    # ratio controls the ellipticity, ratio=1 is a disc
    x = np.sqrt((xr-param.Lx*x0)**2+(yr-param.Ly*y0)**2*ratio**2)

    y = x.copy()*0.

    if vortex_type in ('gaussian', 'cosine', 'step'):
        if vortex_type == 'gaussian':
            y = np.exp(-x**2/(sigma**2))

        if vortex_type == 'cosine':
            y = np.cos(x/sigma*np.pi/2)
            y[x > sigma] = 0.

        if vortex_type == 'step':
            y[x <= sigma] = 1.
    else:
        print('this kind of vortex (%s) is not defined' % vortex_type)

    return y


# 2/ set an initial tracer field
vtype = 'gaussian'
# vortex width
sigma = 0.1*param.Lx
amplitude = -.5

vortex_config = 'dipole'

if vortex_config == 'single':
    vor[:] = amplitude*vortex(param, grid, 0.5, 0.5, sigma, vtype, ratio=1)

elif vortex_config == 'dipole':
    d = 3*sigma/param.Lx
    vor[:] += +amplitude*vortex(param, grid, 0.5, 0.5-d/2, sigma, vtype)
    vor[:] += +amplitude*vortex(param, grid, 0.5, 0.5+d/2, sigma, vtype)


model.set_psi_from_pv()

f2d.loop()
