from param import Param
from grid import Grid
from fluid2d import Fluid2d
from numpy import exp,sqrt,pi,cos,where,random,shape
from restart import Restart
from island import Island

param = Param('default.xml')
param.modelname = 'euler'
param.expname = 'karman_0'

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
param.cfl=1.
param.adaptable_dt=True
param.dt = 1e-2
param.dtmax=1.

# discretization
param.order=5
param.timestepping='RK3_SSP'

# output
param.var_to_save=['vorticity','psi','dye','age']
param.list_diag=['ke','vorticity','enstrophy']
param.freq_plot=10
param.freq_his =1.
param.freq_diag=.1

#plot
param.plot_interactive=True
param.plot_var='vorticity'
param.cax=[-20,20]
param.colorscheme='minmax'

# physics
param.forcing = False
param.noslip = True
param.diffusion=False
param.forcing = False
param.isisland=True


param.additional_tracer=['dye','age']

nh = param.nh

grid  = Grid(param)

def disc(param,grid,x0,y0,sigma):
    xr,yr = grid.xr,grid.yr
    r = sqrt( (xr-param.Lx*x0)**2+(yr-param.Ly*y0)**2 )
    # add a fin
    #r[66:68,30:60]=0
    return where(r<=sigma)

def square(param,grid,x0,y0,sigma):
    xr,yr = grid.xr,grid.yr
    rx = abs(xr-param.Lx*x0)/sigma
    ry = abs(yr-param.Ly*y0)/sigma
    return where((rx<1) & (ry<1))

psi0 = -0.2 # this sets the inflow intensity
sigma = 0.03
idx=disc(param,grid,sigma*4,0.5,sigma)
grid.msk[idx]=0
grid.island.add(idx,psi0*.5)

if grid.j0==param.npy-1:
    msk = grid.msk.copy()*0
    msk[-nh:-1,:]=1
    idx = where(msk==1)
    grid.island.add(idx,psi0)

# tell the code to deactivate the no-slip along the outer walls
if grid.j0==0:
    grid.msknoslip[:nh,:]=1
if grid.j0==param.npy-1:
    grid.msknoslip[-nh:,:]=1



param.Kdiff=20e-3*grid.dx # regular street, with Heun and up3
param.Kdiff=10e-3*grid.dx # regular street, single vortices only
param.Kdiff=5e-3*grid.dx 
param.Kdiff=0 # start to send dipolewith up5, RK3, no dissipation, all good too


f2d = Fluid2d(param,grid)
model = f2d.model

xr,yr = grid.xr,grid.yr
vor = model.var.get('vorticity')

# add noise to trigger the instability
noise = random.normal(size=shape(yr))*grid.msk
grid.fill_halo(noise)
noise -= grid.domain_integration(noise)*grid.msk/grid.area

vor += 1e-7*noise*grid.msk

vor *= grid.msk



model.set_psi_from_vorticity()

f2d.loop()

