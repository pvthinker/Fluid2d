from fluid2d import Fluid2d
from param import Param
from grid import Grid
import numpy as np
import wavepackets as wp

# If the code immediately stops with

# Traceback (most recent call last):
#   File "vortex.py", line 1, in <module>
#     from param import Param
# ImportError: No module named param

# it means that you forgot to do
# source activate.sh in your terminal


param = Param('default.xml')
param.modelname = 'euler'
param.expname = 'merging_00'

# domain and resolution
param.nx = 64*2
param.ny = param.nx
param.Ly = param.Lx
param.npx = 1
param.npy = 1
param.geometry = 'closed'

# time
param.tend = 20.
param.cfl = 1.2
param.adaptable_dt = True
param.dt = 0.01
param.dtmax = 100

# discretization
param.order = 5
param.timestepping = 'RK3_SSP'
param.exacthistime = True

# output
param.var_to_save = ['vorticity', 'psi', 'tracer']
param.list_diag = 'all'
param.freq_plot = 10
param.freq_his = .2
param.freq_diag = 0.02

# plot
param.freq_plot = 10
param.plot_interactive = True
param.plot_psi = True
param.plot_var = 'vorticity'
param.cax = np.array([-2, 2.])*2
param.colorscheme = 'imposed'
param.generate_mp4 = True

# physics
param.noslip = False
param.diffusion = False
param.additional_tracer = ['tracer']

grid = Grid(param)

param.Kdiff = 5e-2*grid.dx

xr, yr = grid.xr, grid.yr


f2d = Fluid2d(param, grid)
model = f2d.model


vor = model.var.get('vorticity')


def vortex(param, grid, x0, y0, sigma,
           vortex_type, ratio=1):
    """Setup a compact distribution of vorticity

    at location x0, y0 vortex, width is sigma, vortex_type controls
    the radial vorticity profile, ratio controls the x/y aspect ratio
    (for ellipses)

    """
    xr, yr = grid.xr, grid.yr
    # ratio controls the ellipticity, ratio=1 is a disc
    x = np.sqrt((xr-param.Lx*x0)**2+(yr-param.Ly*y0)**2*ratio**2)

    y = x.copy()*0.

    if vortex_type in ('gaussian', 'cosine', 'step', 'square', 'triangle'):
        if vortex_type == 'gaussian':
            y = np.exp(-x**2/(sigma**2))

        elif vortex_type == 'cosine':
            y = np.cos(x/sigma*np.pi/2)
            y[x > sigma] = 0.

        elif vortex_type == 'step':
            y[x <= sigma] = 1.

        elif vortex_type == 'square':
            y[msk] = wp.square(xr-param.Lx*x0, yr-param.Ly*y0, sigma)

        elif vortex_type == 'triangle':
            y = wp.triangle(xr-param.Lx*x0, yr-param.Ly*y0, sigma)
    else:
        raise ValueError('this kind of vortex (%s) is not defined' % vortex_type)

    return y


# 2/ set an initial tracer field
vtype = 'gaussian'
# vortex width
sigma = 0.05*param.Lx

vortex_config = 'corotating'

if vortex_config == 'single':
    vor[:] = vortex(param, grid, 0.5, 0.5, sigma,
                    vtype, ratio=1)


elif vortex_config == 'dipole':
    x0 = 0.5
    vor[:] = -vortex(param, grid, x0, 0.42, sigma,
                     vtype, ratio=1)
    vor[:] += vortex(param, grid, x0, 0.58, sigma,
                     vtype, ratio=1)

elif vortex_config == 'corotating':
    dist = 4*sigma
    vor[:] = vortex(param, grid, 0.5, 0.5+dist/2, sigma, vtype)
    vor[:] += vortex(param, grid, 0.5, 0.5-dist/2, sigma, vtype)



vor[:] = vor*grid.msk

if False:
    np.random.seed(1)  # this guarantees the results reproducibility
    noise = np.random.normal(size=np.shape(yr))*grid.msk
    noise -= grid.domain_integration(noise)*grid.msk/grid.area
    grid.fill_halo(noise)

    noise_amplitude = 1e-3

    vor += noise*noise_amplitude

model.set_psi_from_vorticity()

state = model.var.get('tracer')
state[:] = np.round(xr*6) % 2 + np.round(yr*6) % 2
state *= grid.msk

# % normalization of the vorticity so that enstrophy == 1.
model.diagnostics(model.var, 0)
enstrophy = model.diags['enstrophy']
# print('enstrophy = %g' % enstrophy)
vor[:] = vor[:] / np.sqrt(enstrophy)
model.set_psi_from_vorticity()


f2d.loop()
