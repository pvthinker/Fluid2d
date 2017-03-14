from param import Param
from grid import Grid
from fluid2d import Fluid2d
from numpy import exp,sqrt,pi,cos,random,shape,tanh

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'gravcurrent_0'

# domain and resolution

param.ny = 128
param.nx = param.ny*2
param.npx=1
param.npy=1
param.Lx = 2.
param.Ly = 1.
param.geometry='square'

#time
param.tend=50
param.cfl=1.2
param.adaptable_dt=True
param.dt = .1
param.dtmax=.1

# discretization
param.order=5
#param.timestepping='LFAM3'

# output
param.plot_var='buoyancy'
param.var_to_save=['vorticity','buoyancy','psi']
param.list_diag=['maxspeed','ke','pe','energy','vorticity','enstrophy','buoyancy','brms']
param.freq_his=.2
param.freq_diag=.1

# plot
param.plot_interactive=True
#param.plotting_module='plotting_rayleigh'
param.freq_plot=10
param.colorscheme='imposed'
param.cax=[0,3]

# physics
param.gravity=1.
param.forcing = False
param.diffusion=True
param.noslip=False

grid  = Grid(param)
param.Kdiff=2e-3*grid.dx**2

xr,yr = grid.xr,grid.yr
# add a mask
hb  = exp( -(xr/param.Lx-0.5)**2 *50)*0.7
grid.msk[(yr<hb)]=0
grid.finalize_msk()

f2d = Fluid2d(param,grid)
model = f2d.model


buoy = model.var.get('buoyancy')


def sigmoid(x,delta):
    return 1/(1+exp(-(x-0.5)/delta))

def stratif():    
    sigma = 3*grid.dx # width of the interface
    b= sigmoid(xr/param.Lx,sigma/param.Lx)
    return (1-b) * grid.msk

# to have a lock exchange type flow
buoy[:]= (1-stratif() - 0.5)*grid.msk

buoy[:]= (yr *(1 + 1*(1+tanh( (xr-param.Lx/2)/0.1)))) * grid.msk

# add noise to trigger the instability
noise = random.normal(size=shape(yr),scale=1.)*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area
grid.fill_halo(noise)

buoy +=  1e-8*noise

model.set_psi_from_vorticity()


f2d.loop()


