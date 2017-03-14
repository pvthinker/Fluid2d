from param import Param
from grid import Grid
from fluid2d import Fluid2d
from numpy import exp,sqrt,pi,cos
from restart import Restart

param = Param('default.xml')
param.modelname = 'euler'
param.expname = 'DW_0'

# domain and resolution
param.nx = 64*4
param.ny = param.nx
param.Ly = param.Lx
param.npx=1
param.npy=1
param.geometry='square'

# time
param.tend=.2
param.cfl=1.2
param.adaptable_dt=True
param.dt = 1e-2
param.dtmax=1e-2

# discretization
param.order=5
param.timestepping='RK3_SSP'

# output
param.var_to_save=['vorticity','psi','v','tauw','wshear']
param.list_diag=['ke','vorticity','enstrophy']
param.freq_plot=10
param.freq_his =0.005
param.freq_diag=0.001

#plot
param.plot_interactive=True
param.plot_var='vorticity'
param.cax=[-2000,2000]
param.colorscheme='imposed'

# physics
param.forcing = False
param.noslip = True
param.diffusion=False

grid  = Grid(param)

param.Kdiff=5e-2*grid.dx

# it's time to modify the mask and add obstacles  if you wish, 0 is land
#grid.msk[:param.ny/4,:param.nx/4]=0
#grid.msk[:55,35]=0
#grid.msk[:,:4]=0
#grid.finalize()

f2d = Fluid2d(param,grid)
model = f2d.model

xr,yr = grid.xr,grid.yr
vor = model.var.get('vorticity')

def vortex(param,grid,x0,y0,sigma):
    xr,yr = grid.xr,grid.yr
    x = sqrt( (xr-param.Lx*x0)**2+(yr-param.Ly*y0)**2 )
    #y = 4./(1+exp(x/(sigma/2)))
    y = exp( -x**2/(2*sigma**2) )
    #y = cos(x/sigma*pi/2)
    #y[x>sigma]=0.
    return y

# 2/ set an initial tracer field
sigma = 0.04*param.Lx
vor[:]-= vortex( param,grid,0.3,0.52,sigma)
vor[:]+= vortex( param,grid,0.3,0.48,sigma)

# to have a chasing dipole (like smoke rings), uncomment below
#vor[:]+= vortex( 0.25,0.55,sigma)
#vor[:]+=-vortex( 0.25,0.45,sigma)
vor[:] = 2e3*vor*grid.msk



model.set_psi_from_vorticity()

f2d.loop()

