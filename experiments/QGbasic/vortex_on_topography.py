from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'quasigeostrophic'
param.expname = 'vortex_topo_0'

# domain and resolution
ratio = 2
param.nx = 64*4
param.ny = param.nx//ratio
param.npx = 1
param.npy = 1
param.Lx = 2.
param.Ly = param.Lx/ratio
param.geometry = 'perio'

# time
param.tend = 50.
param.cfl = 0.8
param.adaptable_dt = True
param.dt = 1.
param.dtmax = 100.

# discretization
param.order = 5
param.flux_splitting_method = 'minmax'

# output
param.ageostrophic = True
param.var_to_save = ['pv', 'psi', 'u', 'v', 'pvanom', 'vorticity',
                     'btorque', 'ua', 'va']
param.list_diag = ['ke', 'pv', 'pv2']
param.freq_his = 1
param.freq_diag = .5

# plot
param.plot_var = 'pv'
param.freq_plot = 10
a = 4.
param.cax = [-a, a]
param.plot_interactive = True
param.plot_psi = True
param.plot_pvback = True
param.colorscheme = 'symmetric'
param.imshow_interpolation = 'bilinear'
param.generate_mp4 = True

# physics
param.beta = 1.
param.Rd = 0.1  # 2*grid.dx
param.bottom_torque = True
param.forcing = False
param.forcing_module = 'forcing_dipoles'
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
vtype = 'cosine'
# vortex width
sigma = 0.06*param.Lx
amplitude = -1.

vortex_config = 'dipole'
topo_config = 'seamount'

remove_pvtopo = True

"""Set the topography"""

if topo_config == 'seamount':
    # define a seamount
    vtype = 'gaussian'
    sigma = 0.08*param.Lx
    htopo = -vortex(param, grid, 0.2, 0.5, sigma, vtype)
    sigma = 0.03*param.Lx
    htopo += 3*vortex(param, grid, 0.22, 0.46, sigma, vtype)

elif topo_config == 'margin':
    htopo = 0.3*(1+np.tanh((yr/param.Ly-0.85)/.1))

elif topo_config == 'canyon':
    steep = 0.02
    stepx = 1+0.2*(np.tanh((xr/param.Lx-0.2)/steep)
                   - np.tanh((xr/param.Lx-0.3)/steep))
    htopo = 0.3*(1+np.tanh((yr/param.Ly+0.5*(stepx-2))/.2))

elif topo_config == 'ridge':
    htopo = 0.3*(np.exp(-(yr/param.Ly-0.6)**2/(2*.05**2)))


"""Set the PV distribution"""

if vortex_config == 'single':
    x0, y0 = 0.5, 0.6
    vor[:] = amplitude*vortex(param, grid, x0, y0, sigma, vtype, ratio=1)

elif vortex_config == 'dipole':
    d = 4*sigma/param.Lx
    x0, y0 = 0.4, 0.5
    vor[:] += -amplitude*vortex(param, grid, x0, y0-d/2, sigma, vtype)
    vor[:] += +amplitude*vortex(param, grid, x0, y0+d/2, sigma, vtype)

elif vortex_config == 'dipole2':
    d = 4*sigma/param.Lx
    vor[:] += -amplitude*vortex(param, grid, 0.5-d/2, 0.2, sigma, vtype)
    vor[:] += +amplitude*vortex(param, grid, 0.5+d/2, 0.2, sigma, vtype)

elif vortex_config == 'jet':
    sigma = 0.05
    yy = (yr-0.5*param.Ly)/sigma
    vor[:] = amplitude*yy*np.exp(-yy**2/2)

elif vortex_config == 'followtopo':
    vor[:] = htopo*amplitude

elif vortex_config == 'random':
    # set an initial small scale random vorticity field (white noise)
    noise = np.random.normal(size=np.shape(yr))*grid.msk
    # noise=model.ope.gmg.generate_random_large_scale_field()

    grid.fill_halo(noise)
    noise -= grid.domain_integration(noise)*grid.msk/grid.area
    vor[:] = noise*amplitude*grid.msk


if remove_pvtopo:
    vor[:] -= htopo


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


# hack the background PV = by-pass the default choice
# (meridional gradient) and set it to htopo
pvback = f2d.model.pvback
pvback[:, :] += htopo
vor += pvback

f2d.plotting.pvback = pvback

# default way
# model.add_backgroundpv()

model.set_psi_from_pv()

f2d.loop()
