from param import Param
from grid import Grid
from fluid2d import Fluid2d
from numpy import exp, sqrt, pi, cos, random, shape

param = Param('default.xml')
param.modelname = 'quasigeostrophic'
param.expname = 'vortex_beta_0'

# domain and resolution
param.nx = 64*2
param.ny = 64
param.npy = 1
param.Lx = 2.
param.Ly = param.Lx/2
param.geometry = 'xchannel'

# time
param.tend = 10e3
param.cfl = 1.
param.adaptable_dt = True
param.dt = 1.
param.dtmax = 100.

# discretization
param.order = 5


# output
param.var_to_save = ['pv', 'psi']
param.list_diag = ['ke', 'pv', 'pv2']
param.freq_his = 50
param.freq_diag = 10

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
param.Rd = 0.01  # 2*grid.dx
param.forcing = False
param.noslip = False
param.diffusion = False

grid = Grid(param)
param.Kdiff = 0.1e-3*grid.dx

# add an island
# grid.msk[28:32,34:38]=0
# grid.finalize()


f2d = Fluid2d(param, grid)
model = f2d.model

xr, yr = grid.xr, grid.yr
vor = model.var.get('pv')


def vortex(param, grid, x0, y0, sigma, vortex_type, ratio=1):
    xr, yr = grid.xr, grid.yr
    # ratio controls the ellipticity, ratio=1 is a disc
    x = sqrt((xr-param.Lx*x0)**2+(yr-param.Ly*y0)**2*ratio**2)

    y = x.copy()*0.

    if vortex_type in ('gaussian', 'cosine', 'step'):
        if vortex_type == 'gaussian':
            y = exp(-x**2/(sigma**2))

        if vortex_type == 'cosine':
            y = cos(x/sigma*pi/2)
            y[x > sigma] = 0.

        if vortex_type == 'step':
            y[x <= sigma] = 1.
    else:
        print('this kind of vortex (%s) is not defined' % vortex_type)

    return y

# 2/ set an initial tracer field


vtype = 'gaussian'
# vortex width
sigma = 0.08*param.Lx

vortex_config = 'single'

if vortex_config == 'single':
    vtype = 'gaussian'
    sigma = 0.05*param.Lx
    amplitude = -0.2
    vor[:] = amplitude*vortex(param, grid, 0.5, 0.5, sigma, vtype, ratio=1)


model.add_backgroundpv()

model.set_psi_from_pv()

f2d.loop()
