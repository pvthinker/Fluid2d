import numpy as np
from param import Param
from grid import Grid
from fluid2d import Fluid2d
from scipy.special import erf
from restart import Restart
from island import Island

param = Param('default.xml')
param.modelname = 'thermalwind'
param.expname = 'SI_00'


# domain and resolution
ratio = 2
param.ny = 64*2
param.nx = param.ny*ratio
param.Ly = 1.
param.Lx = param.Ly*ratio
param.npx = 1
param.npy = 1
param.geometry = 'closed'

# time
param.tend = 150.
param.cfl = 1.0
param.adaptable_dt = True
param.dt = .1
param.dtmax = 2.

# discretization
param.order = 3
param.timestepping = 'RK3_SSP'

# output
param.var_to_save = ['vorticity', 'psi',
                     'V', 'buoyancy', 'qE', 'tracer']  # qE is Ertel PV
param.list_diag = 'all'
param.freq_plot = 4
param.freq_his = 5.
param.freq_diag = .1
param.diag_fluxes = True

# plot
param.generate_mp4 = True
param.plot_interactive = True
param.plot_var = 'qE'
param.cax = np.asarray([-0.1, .1])*2
param.colorscheme = 'imposed'
# comment param.plotting_module to restore the default plotting
param.plotting_module = 'plotting_SI'

# physics
param.forcing = False
param.noslip = False
param.diffusion = False
param.forcing = False
param.generate_mp4 = True

param.gravity = 1.
param.f0 = 0.1
param.additional_tracer = ['tracer']


nh = param.nh

grid = Grid(param)

param.Kdiff = 5e-1*grid.dx**2

f2d = Fluid2d(param, grid)
model = f2d.model

xr, yr = grid.xr, grid.yr
vor = model.var.get('vorticity')
buoy = model.var.get('buoyancy')
V = model.var.get('V')
qE = model.var.get('qE')
trac = model.var.get('tracer')

sigma = 0.2
delta = .2
V0 = .2  # jet speed

N2 = .15  # vertical   buoyancy gradient
delta = .4
sigma = .4

# horizontal buoyancy gradient (with the gaussian jet below)
M2 = param.f0*V0/delta

# condition for SI, q>0 (Stone 1966, Taylor-Ferrari 2009)
q = M2**2/param.f0**2 - N2

# geostrophic Richardson number (Whitt & Thomas 2013)
Rig = param.f0**2 * N2 / M2**2

# Richardson angle (Thomas et al 2012)
phiRiB = np.arctan(- 1/Rig)*180/np.pi
print('Richardson angle : %g' % phiRiB)

x = (grid.xr0/sigma)

# derivative of erf(x) is exp(-x^2) * 2 / sqrt(pi)

# thermal wind balance reads  \partial_x b - f \partial_z V = 0 (if V is along y)

# background stratification
z = grid.yr0/param.Ly

bback = N2*z  # <- uniform stratification
#bback[z<0] = 10*N2*z[z<0]

config = 'symmetric'

if config == 'centrifugal_surface':
    # surface intensified jet
    z = (grid.yr-param.Ly)/delta
    V[:, :] = V0*np.exp(-x**2) * np.exp(z)
    buoy[:, :] = V0*(sigma/delta)*(param.f0) * erf(x) * \
        (np.sqrt(pi)/2) * np.exp(z) + bback

elif config == 'centrifugal_interior':
    # it's centrifugal instability because the vertical component of the vorticity is negative (Thomas et al 2012)

    # interior jet (with max speed at mid-depth)
    z = (grid.yr-param.Ly*.5)/delta
    V[:, :] = V0*np.exp(-x**2) * np.exp(-z**2)
    buoy[:, :] = V0*(sigma/delta)*(param.f0) * erf(x) * \
        (np.sqrt(np.pi)/2) * (-2*z)*np.exp(-z**2) + bback

elif config == 'noname':
    z = (grid.yr-param.Ly*.5)/delta
#    V[:, :] = V0*np.tanh(x) * exp(-z**2)
#    buoy[:, :] = V0 * np.log(np.cosh(x)) * (-2*z)*exp(-z**2) + bback
    V[:, :] = V0/np.cosh(x)**2 * np.exp(-z**2)
    buoy[:, :] = (V0*param.f0) * np.tanh(x) * (-2*z)*np.exp(-z**2) + bback

elif config == 'symmetric':
    z = grid.yr0/delta
    x = (grid.xr0)/sigma
    alpha = 1.
    V[:, :] = V0*(1+np.tanh(z))*(1-alpha*2*x*np.exp(-x**2))
    buoy[:, :] = (V0*param.f0) / np.cosh(z)**2 * \
        (x+alpha*np.exp(-x**2)) + bback

elif config == 'noname2':
    def phi(x_, z_):
        return np.exp(-.5*(x_+z_)**2)
    z = grid.yr0
    x = grid.xr0
    eps = 0.1*grid.dx

    V0 = .5
    V[:, :] = V0*(phi(x+eps, z)-phi(x-eps, z))/(2*eps) + 0.5*x
    buoy[:, :] = (V0*param.f0)*(phi(x, z+eps)-phi(x, z-eps))/(2*eps)+bback

else:
    z = grid.yr0/delta
    x = grid.xr0/sigma
    alpha = 5.
    phi = x + alpha*z
    V[:, :] = V0*phi
    buoy[:, :] = (V0*param.f0) * (alpha+0.5*x**2) + bback

V *= grid.msk
buoy *= grid.msk


vor[:, :] = 0.

# add noise to trigger the instability
np.random.seed(1)
noise = np.random.normal(size=np.shape(yr))*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area
grid.fill_halo(noise)

vor[:, :] += 1e-1*noise*grid.msk

vor *= grid.msk


model.set_psi_from_vorticity()

model.compute_pv()
# set the tracer equals to initial PV!
trac[:, :] = qE

f2d.loop()
