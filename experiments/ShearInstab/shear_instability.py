from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'euler'
param.expname = 'shear_instab_00'

# domain and resolution
param.nx = 64*2*2
param.ny = 64*2
param.Lx = 2.
param.Ly = param.Lx/2
param.geometry = 'xchannel'

# time
param.tend = 50.
param.cfl = 1.5
param.adaptable_dt = True
param.dt = 1.
param.dtmax = 10.

# discretization
param.order = 5

# output
param.var_to_save = ['vorticity', 'psi', 'u', 'v', 'tracer']
param.list_diag = 'all'
param.freq_plot = 10
param.freq_his = .5
param.freq_diag = .1

# plot
param.plot_interactive = True
param.plot_var = 'vorticity'
param.plot_psi = False
param.cax = [-1.1, 1.1]
param.colorscheme = 'imposed'
param.generate_mp4 = True

# physics
param.forcing = False
param.noslip = False
param.diffusion = False
param.additional_tracer = ['tracer']

grid = Grid(param)
param.Kdiff = 2e-4*grid.dx

# it's time to modify the mask and add obstacles  if you wish, 0 is land
# grid.msk[:55,35]=0
# grid.msk[:,:4]=0
# grid.finalize_msk()

f2d = Fluid2d(param, grid)
model = f2d.model

xr, yr = grid.xr, grid.yr
vor = model.var.get('vorticity')


def vortex(x0, y0, sigma):
    x = np.sqrt((xr-param.Lx*x0)**2+(yr-param.Ly*y0)**2)
    y = 4./(1+np.exp(x/(sigma/2)))
    # y = cos(x/sigma*pi/2)
    # y[x>sigma]=0.
    return y


# 2/ set an initial shear
sigma = 0.07
yy = (yr/param.Ly-0.5)

# comment/uncomment choice A vs. B vs. C

# choice A/ corresponding to a gaussian shaped jet
# vor[:] = (yy/sigma)*exp( - (yy/sigma)**2/2 )

# choice B/ or a cosine shaped jet
# vor[:] = np.sin(yy*np.pi/sigma)
# vor[abs(yy/sigma) > 1] = 0.

# choice C/ or a piecewise constant jet
vor[:] = np.sign(yy)
vor[abs(yy/sigma) > 1] = 0.

# add noise to trigger the instability
np.random.seed(42)
noise = np.random.normal(size=np.shape(yr))*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area
grid.fill_halo(noise)

vor += noise*1e-2

model.set_psi_from_vorticity()

state = model.var.get('tracer')
state[:] = np.round(xr*6) % 2 + np.round(yr*6) % 2

f2d.loop()
