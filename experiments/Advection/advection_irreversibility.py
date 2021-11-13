# irreversibility effects for an advection scheme (space+time)

import numpy as np
from param import Param
from grid import Grid
from fluid2d import Fluid2d
from plotting import Plotting
import matplotlib.pyplot as plt


param = Param('default')
param.modelname = 'advection'
param.expname = 'irrev_0'

# domain and resolution
param.nx = 64*2
param.ny = 64*2
param.npx = 1
param.npy = 1
param.Lx = 1.
param.Ly = 1.
param.geometry = 'perio'

# time
param.tend = 1.
param.cfl = 1.2
param.adaptable_dt = True
param.dt = .1
# param.dtmax=1e4

# *** discretization ***
param.order = 5  # 1,2,3,4,5 : order of the spatial discretization. Even
# cases are centered fluxes, odd cases are upwinded
# fluxes
param.timestepping = 'RK3_SSP'  # LFAM3, Heun, EF, LF, RK3_SSP, AB2,
# AB3. See core/timestepping.py to see the
# list, you may add your own in there

# output
param.var_to_save = ['tracer']
param.list_diag = ['mean', 'rms']
param.freq_his = 1
param.freq_diag = 1

# plot
param.plot_var = 'tracer'  # 'psi' or 'tracer'
param.freq_plot = 1
param.colorscheme = 'imposed'
param.cax = [-0.2, 2.2]
param.plot_interactive = True
# param.plotting_module='plotting_adv'

# physics
param.diffusion = False


# ********* PARAMETERS CONTROLLING THE FLOW AND THE TRACER ********
flow_config = 0  # controls the flow
tracer_config = 0  # controls the tracer


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
    angle = 20.
    speed = 0.1
    if param.geometry == 'perio':
        u[:] = speed*np.cos(angle*np.pi/180.)
        v[:] = speed*np.sin(angle*np.pi/180.)
    else:
        print('param.geometry should be perio')
        exit(0)

if flow_config == 1:  # body rotation
    state[:] = 1.
    state[:] = state[:]*grid.msk
    # determine the stream function from the vorticity
    model.set_psi_from_tracer()

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
        state *= 100.
        # determine the stream function from the vorticity
        model.set_psi_from_tracer()
    else:
        print('param.geometry should be perio')
        exit(0)


# normalize all flows with 'maxspeed=1'
f2d.model.diagnostics(f2d.model.var, 0.)
maxspeed = f2d.model.diags['maxspeed']
u = u/maxspeed
v = v/maxspeed
psi = psi/maxspeed


# 2/ set an initial tracer field
sigma = 0.05*param.Lx


if tracer_config == 0:  # single localized patch
    state[:] = 2*vortex(0.3, 0.4, sigma)

if tracer_config == 1:  # tiles
    state[:] = round(xr*4) % 2 + round(yr*4) % 2

if tracer_config == 2:  # a closed line
    x = xr/param.Lx
    y = yr/param.Ly
    msk = np.zeros_like(yr)
    msk[(x > 0.3) & (x < 0.7) & (y > 0.3) & (y < 0.7)] = 1
    z = np.roll(msk, -1, axis=1)+np.roll(msk, -1, axis=0) + \
        np.roll(msk, +1, axis=1)+np.roll(msk, +1, axis=0)-4*msk
    state[:] = 0.
    state[z > 0] = 4.


state[:] *= grid.msk


# ***** instead of looping in time we take control of the time steps
nh = param.nh
z2d = f2d.model.var.get(param.plot_var)[nh:-nh, nh:-nh]
msk = grid.msk[nh:-nh, nh:-nh]

t = 0.
var = f2d.model.var


# advect it for nite time step only
nite = 1


z0 = z2d.copy()  # copy initial state

f2d.loop()
# f2d.model.diagnostics(var,0.)
# f2d.set_dt()
# for k in range(nite):
#    f2d.model.step(t,f2d.dt)

# revert the velocity
var.state[2] *= -1.
var.state[3] *= -1.

f2d.t = 0
f2d.kt = 0
f2d.loop()
# for k in range(nite):
#    f2d.model.step(t,f2d.dt)

# substract the initial state
z2d = z2d-z0

print('everything was ok')

# look at the result


plt.ioff()

maxi = max(abs(z2d.flatten()))

fig = plt.figure()
im = plt.imshow(z2d, vmin=-maxi, vmax=maxi, origin='lower')
plt.title('Final state - Initial state')
plt.colorbar(im)
plt.show()
