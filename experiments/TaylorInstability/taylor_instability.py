from param import Param
from grid import Grid
from fluid2d import Fluid2d
from numpy import exp,sqrt,pi,cos,random,shape

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'taylor_instab_0'

# domain and resolution
param.nx = 64*4
param.ny = param.nx
param.npx=1
param.npy=1
param.Lx = 1.
param.Ly = param.Lx
param.geometry='xchannel'

#time
param.tend=10
param.cfl=.7
param.adaptable_dt=False
param.dt = 0.005
param.dtmax=.1

# discretization
param.order=5

# output
param.plot_var='buoyancy'
param.var_to_save=['vorticity','buoyancy','psi']
param.list_diag=['maxspeed','ke','pe','energy','vorticity','enstrophy','buoyancy','brms']
param.freq_his=1
param.freq_diag=1

# plot
param.plot_interactive=True
#param.plotting_module='plotting_rayleigh'
param.freq_plot=10
param.colorscheme='imposed'
param.cax=[-.6,.6]

# physics
param.gravity=1.
param.forcing = False
param.diffusion=False
param.noslip=False
param.enforce_momentum=False

grid  = Grid(param)

param.Kdiff=2e-3*grid.dx


f2d = Fluid2d(param,grid)
model = f2d.model

xr,yr = grid.xr,grid.yr
buoy = model.var.get('buoyancy')


def sigmoid(x,delta):
    return 1/(1+exp(-(x-0.5)/delta))

def stratif():    
    sigma = 3*grid.dx # width of the interface
    b= sigmoid(yr/param.Ly,sigma/param.Lx)
    return b * grid.msk

buoy[:]= (1-stratif() - 0.5)*grid.msk
# add noise to trigger the instability
noise = random.normal(size=shape(yr),scale=1.)*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area
grid.fill_halo(noise)

buoy +=  1e-3*noise

model.set_psi_from_vorticity()


f2d.loop()


