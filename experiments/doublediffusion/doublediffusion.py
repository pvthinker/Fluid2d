from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'boussinesqTS'
# if you want to not store on $HOME/data/fluid2d, uncomment
# param.datadir = '/workdir/myname/data/fluid2d'
param.expname = 'dbl_diff'

# domain and resolution
param.nx = 64*2

param.ny = param.nx*2
param.npx = 1
param.npy = 1
param.Lx = 1.
param.Ly = param.Lx*2
param.geometry = 'xchannel'

# time
param.tend = 100.
param.cfl = 1.
param.adaptable_dt = True
param.dt = 0.01
param.dtmax = 1e-2

# discretization
param.order = 5
param.relaxation = 'tridiagonal'

# output
param.plot_var = 'S'
param.var_to_save = ['vorticity', 'density', 'T', 'S', 'psi']
param.list_diag = 'all'
param.freq_his = 1.
param.freq_diag = 2e-1

# plot
param.plot_interactive = True
param.freq_plot = 10
param.colorscheme = 'minmax'
param.cax = np.asarray([-10., 10.])
param.generate_mp4 = True

# physics
param.gravity = 1.
param.forcing = True # set it to True to have r.h.s on prognostic variables
param.diffusion = True
param.noslip = False
param.alphaT = 1.
param.betaS = 1.
# also modify param.Kdiff with these new tracers
# param.additional_tracer = ['mytracer1', 'phyto']


grid = Grid(param)

xr, yr = grid.xr, grid.yr
Lx, Ly = param.Lx, param.Ly


K0 = 2e-2*grid.dx
param.Kdiff =  {'vorticity': K0*7., 'T': K0, 'S': K0/50.}
# add my own attributes (used in the forcing)
setattr(param, 'dTdz', 50.)
setattr(param, 'dSdz', 50.)

class Forcing:
    """ define the forcing """

    def __init__(self, param, grid):

        self.list_param = ['deltab', 'dt']
        param.copy(self, self.list_param)

        self.list_param = ['j0', 'npy', 'nh']
        grid.copy(self, self.list_param)

        # this is the z coordinate (from -param.Ly/2 to param.Ly/2)
        z = grid.yr0
        dz = grid.dy

        KT = param.Kdiff['T']
        KS = param.Kdiff['S']

        dTdz = param.dTdz
        dSdz = param.dSdz

        self.FluxT = KT*dTdz/dz
        self.FluxS = KS*dSdz/dz
        print('Temperature flux: %.2g' % self.FluxT)

    def add_forcing(self, x, t, dxdt, coef=1.):
        """ add the forcing term on
        - the temperature 'T' (x[6])
        - the salinity 'S' (x[7])
        """

        if self.j0 == self.npy-1:
            # add the forcing to the uppermost row
            # yet, inside the domain (the halo has
            # a width of nh points) and because
            # of MPI domain decomposition do that only
            # subdomains on the top row
            dxdt[6][-self.nh-1, :] += self.FluxT
            dxdt[7][-self.nh-1, :] += self.FluxS

        if self.j0 == 0:
            # same but for the bottom boundary condition
            dxdt[6][self.nh, :] -= self.FluxT
            dxdt[7][self.nh, :] -= self.FluxS





f2d = Fluid2d(param, grid)
model = f2d.model
if param.forcing:
    model.forc = Forcing(param, grid)

# let's define the initial state

x, z = grid.xr0, grid.yr0
# temp, salt and vor are 'pointers' to model arrays
temp = model.var.get('T')
salt = model.var.get('S')
vor = model.var.get('vorticity')

delta = 0.02*param.Lx
C = 50.

zz = z/delta

# assign the values to temp / the [:, :] is VERY important
temp[:, :] = param.dTdz*z
salt[:, :] = param.dSdz*z


RaT = 2*param.gravity*param.alphaT*C*param.Ly**3/(param.Kdiff['vorticity']*param.Kdiff['T'])
print('RaT = %.3g' % RaT)


# add noise to trigger the instability
np.random.seed(42)
noise = np.random.normal(size=np.shape(yr), scale=1.)*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area

grid.fill_halo(noise)

vor[:, :] = 1e-2*noise

model.set_density()
model.set_psi_from_vorticity()

f2d.loop()
