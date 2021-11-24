# test all about advection schemes
from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np


param = Param('default')
param.modelname = 'advection'
param.expname = 'advection_0'

# domain and resolution
param.nx = 64*2
param.ny = 64*2
param.npx = 1
param.npy = 1
param.Lx = 1.
param.Ly = 1.
param.geometry = 'disc'  # 'square', 'perio', 'disc', 'xchannel'

# time
param.tend = 10.
param.cfl = 1.5
param.adaptable_dt = True
param.dt = 1e4
param.dtmax=1e4
param.asselin_cst = 0.1 # to go with the Leap-Frog time stepping

# *** discretization ***
param.order = 5  # order of the spatial discretization, between 1 and 5
# Even orders are centered fluxes
# Odd orders  are upwinded fluxes

param.timestepping = 'RK3_SSP'  # LFAM3, Heun, EF, LF, RK3_SSP, AB2, AB3.
# See core/timestepping.py to see the
# list, you may add your own in there
#
# table of CFL (from Lemarie et al 2015)
#
# caution: their order=4 discretization is not the same as in Fluid2d
#
# |-------------+-------+-------+-------|
# |    order    |   2   |   3   |   4   |
# |-------------+-------+-------+-------|
# | LFRA nu=0.1 | 0.904 | 0.472 | 0.522 |
# | LFAM3       | 1.587 | 0.871 | 0.916 |
# | AB2 eps=0.1 | 0.503 | 0.554 | 0.29  |
# | AB3         | 0.724 | 0.397 | 0.418 |
# | RK3         |  1.73 | 1.626 |     1 |
# |-------------+-------+-------+-------|
#

#
# ********* PARAMETERS CONTROLLING THE FLOW AND THE TRACER ********
#
flow_config = 1 # controls the flow
# 0: translation (it needs param.geometry = 'perio')
# 1: body rotation
# 2: shear
# 3: vortex

tracer_config = 0  # controls the tracer
# 0: isolated patch, width controlled with patch_width
# 1: tiles
# 2: dx width square contour
patch_width = 0.02*param.Lx


# output
param.var_to_save = ['tracer', 'psi']
param.list_diag = ['mean', 'rms']
param.freq_his = .1
param.freq_diag = .1

# plot
param.plot_var = 'tracer'  # 'psi' or 'tracer'
param.freq_plot = 30
param.colorscheme = 'imposed'  # 'imposed' or 'minmax'=self adjust
param.cax = [-0.2, 1.2]
param.plot_interactive = True
param.generate_mp4 = True
param.plotting_module = 'plotting_adv'
param.plot_psi = True

# physics
param.diffusion = False


grid = Grid(param)
param.Kdiff = 1e-3*grid.dx

# put obstacles: uncomment below
# grid.msk[-60:-40,40:60]=0
# grid.msk[:64,64]=0

f2d = Fluid2d(param, grid)
model = f2d.model


def gaussian_shape(x0, y0, sigma):
    x = np.sqrt((xr-param.Lx*x0)**2+(yr-param.Ly*y0)**2)
    y = 2-2./(1+np.exp(-x/sigma))
    # y = cos(x/sigma*pi/2)
    # y[x>sigma]=0.
    return y



# 1/ set the stream function

state = model.var.get('tracer')
u = model.var.get('u')
v = model.var.get('v')
psi = model.var.get('psi')
xr, yr = grid.xr, grid.yr

if flow_config == 0:  # pure translation
    assert param.geometry == 'perio', '[ERROR]: param.geometry should be perio'
    angle = 0.  # with respect to the x axis
    speed = 1.
    u[:] = speed*np.cos(angle*np.pi/180.)
    v[:] = speed*np.sin(angle*np.pi/180.)

elif flow_config == 1:  # body rotation
    assert param.geometry != 'perio', '[ERROR]: param.geometry should not be perio'
    state[:] = 1.
    state[:] = state[:]*grid.msk

elif flow_config == 2:  # horizontal shear
    assert param.geometry == 'xchannel', '[ERROR]: param.geometry should be xchannel'
    state[:] = -1.
    state[:] = state[:]*grid.msk

elif flow_config == 3:  # single vortex
    sigma = grid.dx*5
    state[:] = 100.*gaussian_shape(0.5, 0.5, sigma)

elif flow_config == 4:  # quadrupole
    assert param.geometry == 'perio', '[ERROR]: param.geometry should be perio'
    sigma = 5*grid.dx
    state[:] = gaussian_shape(0.75, 0.75, sigma)  # single vortex
    state[:] += gaussian_shape(0.25, 0.25, sigma)  # single vortex
    state[:] -= gaussian_shape(0.25, 0.75, sigma)  # single vortex
    state[:] -= gaussian_shape(0.75, 0.25, sigma)  # single vortex


elif flow_config == 5:  # a complicated flow
    assert param.geometry in ['perio', 'closed'], '[ERROR]: param.geometry should be perio or closed'
    np.random.seed(42)
    for k in range(8):
        sigma = np.random.uniform(10)*grid.dx
        sign = np.random.randint(2)*2-1
        x0 = np.random.uniform()
        y0 = np.random.uniform()
        state[:] += sign*gaussian_shape(x0, y0, sigma)/np.sqrt(sigma)

else:
    raise ValueError("[ERROR]: undefined flow_config")

model.set_psi_from_tracer()



# 2/ set an initial tracer field

if tracer_config == 0:  # single localized patch
    state[:] = 1.2*gaussian_shape(0.3, 0.4, patch_width)

elif tracer_config == 1:  # tiles
    state[:] = (np.round(xr*6) % 2 + np.round(yr*6) % 2)/2

elif tracer_config == 2:  # a closed line
    x = xr/param.Lx
    y = yr/param.Ly
    msk = np.zeros_like(yr)
    msk[(x > 0.4) & (x < 0.7) & (y > 0.4) & (y < 0.7)] = 1
    z = np.roll(msk, -1, axis=1)+np.roll(msk, -1, axis=0) + \
        np.roll(msk, +1, axis=1)+np.roll(msk, +1, axis=0)-4*msk
    state[:] = 0.
    state[z > 0] = 4.

else:
    raise ValueError("[ERROR]: undefined tracer_config")


state[:] *= grid.msk

dummy = state*1.
grid.fill_halo(dummy)
state[:] = dummy

f2d.loop()

maxspeed = f2d.model.diags['maxspeed']

print(f"time step : {f2d.dt:.3e}")
print(f"max speed : {maxspeed:.3e}")
print(f"grid size : {grid.dx:.3e}")
print(f"CFL       : {param.cfl:.2f}")
