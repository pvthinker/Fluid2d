from param import Param
from grid import Grid
from fluid2d import Fluid2d
from numpy import exp,sqrt,pi,cos,random,shape

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'rayleighbenard_0'

# domain and resolution
param.nx = 64*4
param.ny = param.nx/4
param.npx=1
param.Lx = 4.
param.Ly = param.Lx/4
param.geometry='xchannel'

#time
param.tend=1e3
param.cfl=1.
param.adaptable_dt=True
param.dt = 1.
param.dtmax=.2

# discretization
param.order=5

# output
param.plot_var='buoyancy'
param.var_to_save=['vorticity','buoyancy','v','psi']
param.list_diag=['ke','pe','energy','vorticity','enstrophy','brms']
param.freq_his=1
param.freq_diag=1

# plot
param.plot_interactive=True
param.plotting_module='plotting_rayleigh'
param.freq_plot=10
param.colorscheme='imposed'
param.cax=[-1,1]

# physics
param.gravity=1.
param.forcing = True
param.forcing_module='forcing_rayleigh'
param.diffusion=False
param.noslip = True

grid  = Grid(param)
# Prandtl number is Kvorticity / Kbuoyancy

prandtl = 0.1
param.Kdiff={}
param.Kdiff['vorticity']=0e-3*grid.dx
param.Kdiff['buoyancy']= 8e-2*grid.dx#param.Kdiff['vorticity'] / prandtl


f2d = Fluid2d(param,grid)
model = f2d.model

xr,yr = grid.xr,grid.yr
buoy = model.var.get('buoyancy')


def sigmoid(x):
    return 1/(1+exp(-(x-0.5)*param.nx/4))

def stratif():    
    b= 0.001*(sigmoid(yr/param.Ly + 0.1*cos(pi/4+2*xr/param.Lx*pi))-0.5)
    return b * grid.msk

# add noise to trigger the instability
noise = random.normal(size=shape(yr))*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area
grid.fill_halo(noise)

buoy[:]+=  1e-7*noise

model.set_psi_from_vorticity()

f2d.loop()
