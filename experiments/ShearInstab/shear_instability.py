from param import Param
from grid import Grid
from fluid2d import Fluid2d
from numpy import exp,sqrt,pi,cos,sin,abs,random,shape

param = Param('default.xml')
param.modelname = 'euler'
param.expname = 'shear_instab_0'

# domain and resolution
param.nx = 64*2
param.ny = 64
param.Lx = 2.
param.Ly = param.Lx/2
param.geometry='xchannel'

# time
param.tend=500.
param.cfl=1.
param.adaptable_dt=True
param.dt = 1.
param.dtmax=10.

# discretization
param.order=5

# output
param.var_to_save=['vorticity','psi','u']
param.list_diag=['ke','vorticity','enstrophy']
param.freq_plot=10
param.freq_his=5
param.freq_diag=1

#plot
param.plot_interactive=True
param.plot_var='psi'
param.cax=[-1,1]
param.colorscheme='minmax'

# physics
param.forcing = False
param.noslip = False
param.diffusion=False

grid  = Grid(param)
param.Kdiff=2e-4*grid.dx

# it's time to modify the mask and add obstacles  if you wish, 0 is land
#grid.msk[:55,35]=0
#grid.msk[:,:4]=0

f2d = Fluid2d(param,grid)
model = f2d.model

xr,yr = grid.xr,grid.yr
vor = model.var.get('vorticity')

def vortex(x0,y0,sigma):
    x = sqrt( (xr-param.Lx*x0)**2+(yr-param.Ly*y0)**2 )
    y = 4./(1+exp(x/(sigma/2)))
    #y = cos(x/sigma*pi/2)
    #y[x>sigma]=0.
    return y

# 2/ set an initial shear 
sigma = 0.05
yy= (yr/param.Ly -0.5)

# comment/uncomment choice A vs. B

# choice A/ corresponding to a gaussian shaped jet
#vor[:] = (yy/sigma)*exp( - (yy/sigma)**2/2 ) 

# choice B/ or a cosine shaped jet
vor[:] = sin(yy*pi/sigma)
vor[abs(yy/sigma)>1]=0.

# add noise to trigger the instability
noise = random.normal(size=shape(yr))*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area
grid.fill_halo(noise)

vor += noise  * 1e-3

model.set_psi_from_vorticity()

f2d.loop()
