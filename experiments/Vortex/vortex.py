from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np
import ana_profiles as ap

# If the code immediately stops with

# Traceback (most recent call last):
#   File "vortex.py", line 1, in <module>
#     from param import Param
# ImportError: No module named param

# it means that you forgot to do
# source activate.sh in your terminal


param = Param('default.xml')
param.modelname = 'euler'
param.expname = 'vortex_08'

# domain and resolution
param.nx = 64*2
param.ny = param.nx
param.Ly = param.Lx
param.npx = 1
param.npy = 1
param.geometry = 'closed'

# time
param.tend = 10
param.cfl = 1.
param.adaptable_dt = True
param.dt = 0.01
param.dtmax = 100

# discretization
param.order = 3
param.timestepping = 'RK3_SSP'
param.exacthistime = True
param.diag_fluxes = True

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
param.cax = np.array([-2, 2.])*5
param.colorscheme = 'imposed'
param.generate_mp4 = True

# physics
param.noslip = False
param.diffusion = False
param.additional_tracer = ['tracer']

grid = Grid(param)

param.Kdiff = 5e-2*grid.dx

xr, yr = grid.xr, grid.yr

# it's time to modify the mask and add obstacles  if you wish, 0 is land

msk_config = 'none'  # the other possibility is 'T-wall' or 'bay'

if msk_config == 'bay':
    x0, y0, radius = 0.5, 0.35, 0.2
    y1 = 0.5
    msk2 = ap.vortex(xr, yr, param.Lx, param.Ly,
                     x0, y0, radius, 'step')

    grid.msk[yr < y1] = 0
    grid.msk += np.asarray(msk2, dtype=int)
    grid.msk[grid.msk < 0] = 0
    grid.msk[grid.msk > 1] = 1
    grid.msk[0:1, :] = 0
    grid.finalize_msk()

elif msk_config == 'T-wall':
    i0, j0 = param.nx//2, param.ny//2
    di = int(0.25*param.Lx/grid.dx)
    grid.msk[:j0, i0] = 0
    grid.msk[j0, i0-di:i0+di] = 0
    grid.finalize_msk()

else:
    # do nothing
    pass

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

    if vortex_type in ('gaussian', 'cosine', 'step'):
        if vortex_type == 'gaussian':
            y = np.exp(-x**2/(sigma**2))

        elif vortex_type == 'cosine':
            y = np.cos(x/sigma*np.pi/2)
            y[x > sigma] = 0.

        elif vortex_type == 'step':
            y[x <= sigma] = 1.
    else:
        print('this kind of vortex (%s) is not defined' % vortex_type)

    return y


# 2/ set an initial tracer field
vtype = 'gaussian'
# vortex width
sigma = 0.0*param.Lx

vortex_config = 'dipole2'

if vortex_config == 'single':
    vtype = 'gaussian'
    sigma = 0.03*param.Lx
    vor[:] = vortex(param, grid, 0.4, 0.54, sigma,
                    vtype, ratio=1)

elif vortex_config == 'dipolebay':
    vtype = 'gaussian'
    sigma = 0.03*param.Lx
    y2 = 0.53
    vor[:] = vortex(param, grid, 0.15, y2, sigma,
                    vtype, ratio=1)
    vor[:] -= vortex(param, grid, -0.15, y2, sigma,
                     vtype, ratio=1)

elif vortex_config == 'dipole2':
    vtype = 'gaussian'
    sigma = 0.05*param.Lx
    x0 = 0.7
    vor[:] = -vortex(param, grid, x0, 0.42, sigma,
                     vtype, ratio=1)
    vor[:] += vortex(param, grid, x0, 0.58, sigma,
                     vtype, ratio=1)

elif vortex_config == 'rankine':
    vtype = 'step'
    ring_sigma = 0.2*param.Lx
    ring_amp = 1.
    vor[:] = ring_amp * vortex(param, grid, 0.5, 0.5, ring_sigma,
                               vtype, ratio=1)
    # sigma ring, core = 0.2, 0.135 yields a tripole (with step distribution)
    # sigma ring, core = 0.2, 0.12 yields a dipole (with step distribution)
    core_sigma = 0.173*param.Lx
    core_amp = ring_amp*(ring_sigma**2-core_sigma**2.)/core_sigma**2.
    vor[:] -= (core_amp+ring_amp)*vortex(param, grid, 0.5, 0.5, core_sigma,
                                         vtype, ratio=1)

elif vortex_config == 'dipole':
    vtype = 'gaussian'
    sigma = 0.04*param.Lx
    vor[:] = vortex(param, grid, 0.3, 0.52, sigma, vtype)
    vor[:] -= vortex(param, grid, 0.3, 0.48, sigma, vtype)

elif vortex_config == 'chasing':
    sigma = 0.03*param.Lx
    vtype = 'step'
    vor[:] = vortex(param, grid, 0.3, 0.6, sigma, vtype)
    vor[:] -= vortex(param, grid, 0.3, 0.4, sigma, vtype)
    vor[:] += vortex(param, grid, 0.1, 0.55, sigma, vtype)
    vor[:] -= vortex(param, grid, 0.1, 0.45, sigma, vtype)

elif vortex_config == 'corotating':
    sigma = 0.06*param.Lx
    dist = 0.25*param.Lx
    vtype = 'gaussian'
    vor[:] = vortex(param, grid, 0.5, 0.5+dist/2, sigma, vtype)
    vor[:] += vortex(param, grid, 0.5, 0.5-dist/2, sigma, vtype)

elif vortex_config == 'collection':
    vtype = 'cosine'
    x0 = [0.3, 0.4, 0.6, 0.8]
    y0 = [0.5, 0.5, 0.5, 0.5]
    amplitude = [1, -2, -1, 2]
    width = np.array([1, 0.5, 1, 0.5])*0.04*param.Lx

    for x, y, a, s in zip(x0, y0, amplitude, width):
        vor[:] += a*vortex(param, grid, x, y, s, vtype)

elif vortex_config == 'unequal':
    # Melander, Zabusky, McWilliams 1987
    # Asymmetric vortex merger in two dimensions: Which vortex is 'victorious'?
    s1 = 0.04*param.Lx
    a1 = 1.
    s2 = 0.1*param.Lx
    a2 = 0.2
    vtype = 'cosine'
    vor[:] = a1*vortex(param, grid, 0.5, 0.6, s1, vtype)
    vor[:] += a2*vortex(param, grid, 0.5, 0.4, s2, vtype)


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
