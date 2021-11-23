# test all about advection schemes
from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np
#from numpy import arange, exp, min, max, sqrt, nan, round, pi, cos, sin, roll, zeros_like, random


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
param.geometry = 'disc'  # 'square', 'perio', 'disc'

# time
param.tend = 10.
param.cfl = 1.5
param.adaptable_dt = True
param.dt = 1e4
param.dtmax=1e4
param.asselin_cst = 0.1

# *** discretization ***
param.order = 5  # 1,2,3,4,5 : order of the spatial discretization. Even
# cases are centered fluxes, odd cases are upwinded
# fluxes
param.timestepping = 'RK3_SSP'  # LFAM3, Heun, EF, LF, RK3_SSP, AB2, AB3.
# See core/timestepping.py to see the
# list, you may add your own in there
#
# table of CFL (from Lemarie et al 2015)
#
# their order=4 discretization is not the same as in Fluid2d
#
# |-------------+-------+-------+-------|
# |    order    |   2   |   3   |   4   |
# |-------------+-------+-------+-------|
# | LFRA nu=0.1 | 0.904 | 0.472 | 0.522 | use 'LF', the Robert-Asselin constant is set to 0.1
# | LFAM3       | 1.587 | 0.871 | 0.916 |
# | AB2 eps=0.1 | 0.503 | 0.554 | 0.29  | 'eps' is set to 0.1 in core/timestepping.py
# | AB3         | 0.724 | 0.397 | 0.418 |
# | RK3         |  1.73 | 1.626 |     1 |
# |-------------+-------+-------+-------|
#

#
# ********* PARAMETERS CONTROLLING THE FLOW AND THE TRACER ********
#
flow_config = 1  # controls the flow (0=translation, 1=body rotation, 2=shear, 3=vortex)
tracer_config = 0  # controls the tracer (0 = isolated patch, 1=tiles, 2=square contour)


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


def vortex(x0, y0, sigma):
    x = np.sqrt((xr-param.Lx*x0)**2+(yr-param.Ly*y0)**2)
    y = 2-2./(1+np.exp(-x/sigma))
    #y = cos(x/sigma*pi/2)
    # y[x>sigma]=0.
    return y


xr, yr = grid.xr, grid.yr
# 1/ set the stream function

state = model.var.get('tracer')
u = model.var.get('u')
v = model.var.get('v')
psi = model.var.get('psi')

if flow_config == 0:  # pure translation
    angle = 0.  # with respect to the x axis
    speed = 1.
    if param.geometry == 'perio':
        u[:] = speed*np.cos(angle*np.pi/180.)
        v[:] = speed*np.sin(angle*np.pi/180.)
    else:
        print('param.geometry should be perio')
        exit(0)

if flow_config == 1:  # body rotation
    if param.geometry != 'perio':
        state[:] = 1.
        state[:] = state[:]*grid.msk
        # determine the stream function from the vorticity
        model.set_psi_from_tracer()
    else:
        print('param.geometry should not be perio')
        exit(0)

if flow_config == 2:  # horizontal shear
    if param.geometry == 'xchannel':
        state[:] = -1.
        state[:] = state[:]*grid.msk
        # determine the stream function from the vorticity
        model.set_psi_from_tracer()
    else:
        print('param.geometry should be xchannel')
        exit(0)

if flow_config == 3:  # single vortex
    sigma = grid.dx*5
    state[:] = 100.*vortex(0.5, 0.5, sigma)
    # determine the stream function from the vorticity
    model.set_psi_from_tracer()


if flow_config == 4:  # quadrupole
    if param.geometry == 'perio':
        sigma = 5*grid.dx
        state[:] = vortex(0.75, 0.75, sigma)  # single vortex
        state[:] += vortex(0.25, 0.25, sigma)  # single vortex
        state[:] -= vortex(0.25, 0.75, sigma)  # single vortex
        state[:] -= vortex(0.75, 0.25, sigma)  # single vortex
        # determine the stream function from the vorticity
        model.set_psi_from_tracer()
    else:
        print('param.geometry should be perio')
        exit(0)

if flow_config == 5:  # a complicated flow

    for k in range(6):
        sigma = np.random.uniform(10)*grid.dx
        x0 = np.random.uniform()
        y0 = np.random.uniform()
        state[:] += vortex(x0, y0, sigma)  # single vortex
        # determine the stream function from the vorticity
        model.set_psi_from_tracer()


# normalize all flows with 'maxspeed=1'
f2d.model.diagnostics(f2d.model.var, 0.)
maxspeed = f2d.model.diags['maxspeed']
# u = u/maxspeed
# v = v/maxspeed
# psi = psi/maxspeed


# 2/ set an initial tracer field
sigma = 0.02*param.Lx


if tracer_config == 0:  # single localized patch
    state[:] = 1.2*vortex(0.3, 0.4, sigma)

if tracer_config == 1:  # tiles
    state[:] = (np.round(xr*6) % 2 + np.round(yr*6) % 2)/2

if tracer_config == 2:  # a closed line
    x = xr/param.Lx
    y = yr/param.Ly
    msk = np.zeros_like(yr)
    msk[(x > 0.4) & (x < 0.7) & (y > 0.4) & (y < 0.7)] = 1
    z = np.roll(msk, -1, axis=1)+np.roll(msk, -1, axis=0) + \
        np.roll(msk, +1, axis=1)+np.roll(msk, +1, axis=0)-4*msk
    state[:] = 0.
    state[z > 0] = 4.


state[:] *= grid.msk

dummy = state*1.
grid.fill_halo(dummy)
state[:] = dummy

f2d.loop()

print(f"pas de temps : {f2d.dt:.3e}")
print(f"vitesse max  : {maxspeed:.3e}")
print(f"grid size dx : {grid.dx:.3e}")
print(f"CFL          : {param.cfl:.2f}")
