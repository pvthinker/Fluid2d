from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'boussinesq'
param.expname = 'HorConv_04'

# domain and resolution
param.nx = 64*4
param.ny = param.nx/4
param.npx = 1
param.Lx = 4.
param.Ly = 1.
param.geometry = 'closed'

# time
param.tend = 4000.
param.cfl = 0.4
param.adaptable_dt = True
param.dt = .1
param.dtmax = .1
param.exacthistime = False

# discretization
param.order = 5

# output
param.plot_var = 'buoyancy'
param.var_to_save = ['vorticity', 'buoyancy', 'v', 'psi']
param.list_diag = 'all'
param.freq_his = 20
param.freq_diag = 1

# plot
param.plot_interactive = True
param.freq_plot = 50
param.colorscheme = 'imposed'
param.cax = [-.5, .5]
param.generate_mp4 = True

# physics
param.gravity = 1.
param.forcing = True
#param.forcing_module = 'forcing_rayleigh'
param.diffusion = True
param.noslip = False

grid = Grid(param)
# Prandtl number is Kvorticity / Kbuoyancy

prandtl = 1.


deltab = 600  # this is directly the Rayleigh number is L=visco=diffus=1

param.deltab = deltab  # make it known to param

L = param.Ly
visco = .002*grid.dy
diffus = visco / prandtl

# Rayleigh number is
Ra = deltab * L ** 3 / (visco * diffus)
print('Rayleigh number is %4.1f' % Ra)

param.Kdiff = {}
param.Kdiff['vorticity'] = visco
param.Kdiff['buoyancy'] = diffus  # param.Kdiff['vorticity'] / prandtl

# time step is imposed by the diffusivity
# param.dt = 0.25*grid.dx**2 / max(visco, diffus)
# param.dtmax = param.dt
# print('dt = %f' % param.dt)


class Forcing:
    """ define the forcing """

    def __init__(self, param, grid):

        self.list_param = ['deltab', 'dt']
        param.copy(self, self.list_param)

        self.list_param = ['j0', 'npy', 'nh']
        grid.copy(self, self.list_param)

        yr = (grid.yr-param.Ly*.5)/param.Ly

        Q = 1e-2
        nh = param.nh
        self.forc = yr*0.

        type_forcing = 'horizontalconvection_v2'

        if type_forcing == 'coolroof':
            if grid.j0 == grid.npy-1:
                self.forc[-nh-1, :] = -Q
            if grid.j0 == 0:
                self.forc[nh, :] = +Q

        elif type_forcing == 'spreadcooling':
            # bottom heating
            if grid.j0 == 0:
                self.forc[nh, :] = +Q
            # uniform cooling
            self.forc -= Q/param.ny

        elif type_forcing == 'horizontalconvection':
            if grid.j0 == grid.npy-1:
                xr = grid.xr
                self.forc[xr < param.Lx/2] = Q
                self.forc[xr > param.Lx/2] = -Q
                self.forc[:-nh-1, :] = 0.

        elif type_forcing == 'horizontalconvection_v2':
            if grid.j0 == grid.npy-1:
                xr = grid.xr
                x0 = param.Lx*3/4
                coef = 1.#0.5/0.8
                self.forc[xr < x0] = Q*coef
                self.forc[xr > x0] = -Q*coef*3
                self.forc[:-nh-1, :] = 0.

        self.type_forcing = type_forcing

        self.forc *= grid.msk

        # transform the surface flux into a volume flux
        self.forc *= (1./grid.dx)

    def add_forcing(self, x, t, dxdt, coef=1.):
        """ add the forcing term on x[0]=the vorticity """
        if self.type_forcing == 'prescribed_buoyancy':
            #
            # for this experiment you need to use the fixed time step
            # whose value is controlled by max(visco,diffus)
            #
            if self.j0 == self.npy-1:
                dxdt[4][-self.nh-1, :] -= (x[4][-self.nh-1, :]
                                           + self.deltab*.5)/self.dt
            if self.j0 == 0:
                dxdt[4][self.nh, :] -= (x[4][self.nh, :]
                                        - self.deltab*.5)/self.dt

        else:
            dxdt[4] += self.forc
        dxdt[4] *= coef



f2d = Fluid2d(param, grid)
model = f2d.model

if param.forcing:
    model.forc = Forcing(param, grid)


xr, yr = grid.xr, grid.yr

buoy = model.var.get('buoyancy')
# add noise to trigger the instability
noise = np.random.normal(size=np.shape(yr))*grid.msk
noise -= grid.domain_integration(noise)*grid.msk/grid.area
grid.fill_halo(noise)

buoy += 1e-1*noise

model.set_psi_from_vorticity()

f2d.loop()
