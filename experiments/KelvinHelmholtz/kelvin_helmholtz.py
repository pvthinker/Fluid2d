from param import Param
from grid import Grid
from fluid2d import Fluid2d
from numpy import exp,sqrt,pi,cos,where,random,shape,tanh,cumsum,cosh
from restart import Restart
from island import Island

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'KH_0'

# domain and resolution
ratio = 2
param.ny = 64*2
param.nx =  param.ny*ratio
param.Ly = 1.  
param.Lx = param.Ly*ratio
param.npx=1
param.npy=1
param.geometry='xchannel'

# time
param.tend=600
param.cfl=0.9
param.adaptable_dt=True
param.dt = 1e-2
param.dtmax=1.

# discretization
param.order=5
param.timestepping='RK3_SSP'

# output
param.var_to_save=['vorticity','psi','u','v','buoyancy']
param.list_diag=['ke','pe','vorticity','enstrophy']
param.freq_plot=20
param.freq_his =1.
param.freq_diag=.1

#plot
param.plot_interactive=True
param.plot_var='v'
param.cax=[-5,5]
param.colorscheme='minmax'

# physics
param.forcing = False
param.noslip = False
param.diffusion=False
param.forcing = False

param.gravity = 1.

nh = param.nh

grid  = Grid(param)


f2d = Fluid2d(param,grid)
model = f2d.model

xr,yr = grid.xr,grid.yr
vor = model.var.get('vorticity')
buoy = model.var.get('buoyancy')

# control parameters of the experiment
N2 = 8. # Brunt Vaisala frequency squared
S2 = 40. # vertical shear squared

# linear stratification
buoy[:,:] = grid.yr0*grid.msk*N2


# we take the mean jet profile as a tanh
#U = tanh( grid.yr0 * sqrt(S2) ) 

vor[:,:] =  -sqrt(S2) / cosh( grid.yr0 * sqrt(S2) )**2.*grid.msk


# add noise to trigger the instability
noise = random.normal(size=shape(yr))*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area
grid.fill_halo(noise)

vor[:,:] += 1e-2*noise*grid.msk

vor *= grid.msk



model.set_psi_from_vorticity()


f2d.loop()

