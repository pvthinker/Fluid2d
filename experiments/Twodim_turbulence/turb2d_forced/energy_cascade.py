from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'euler'
param.expname = 'turb2d_forced_ss'

# domain and resolution
param.nx = 64*2
param.ny = param.nx
param.Ly = param.Lx
param.npx = 1
param.npy = 1
param.geometry = 'perio'

# time
param.tend = 500.
param.cfl = 1.2
param.adaptable_dt = True
param.dt = .05
param.dtmax = 10.

# discretization
param.order = 5

# output
param.var_to_save = ['vorticity', 'psi', 'tracer']
param.list_diag = 'all'
param.freq_plot = 5
param.freq_his = 20
param.freq_diag = 10
param.exacthistime = True

# plot
param.plot_interactive = True
param.plot_var = 'vorticity'
param.plot_psi = False
param.cmap = 'Spectral'
param.cax = np.array([-1, 1])*.3
param.colorscheme = 'imposed'
param.generate_mp4 = True

# physics
param.forcing = True
param.forcing_module = 'embedded'
param.tau = 1.
param.k0 = 0.5*param.nx
param.noslip = False
param.diffusion = False
param.decay = False
#
# add a passive tracer
param.additional_tracer = ['tracer']


class Forcing(Param):
    """ define the forcing """

    def __init__(self, param, grid):
        self.list_param = ['sizevar', 'tau', 'k0']
        param.copy(self, self.list_param)

        self.list_param = ['nh', 'msk', 'fill_halo', 'nx', 'ny',
                           'domain_integration', 'area']
        grid.copy(self, self.list_param)

        def set_x_and_k(n, L):
            k = ((n//2+np.arange(n)) % n) - n//2
            return (np.arange(n)+0.5)*L/n, 2*np.pi*k/L

        self.t0 = 0.

        _, kx = set_x_and_k(param.nx, np.pi)
        _, ky = set_x_and_k(param.ny, np.pi)
        kkx, kky = np.meshgrid(kx, ky)
        self.kk = np.sqrt(kkx**2 + kky**2)

        k0 = self.k0
        dk = .1

        # the forcing has a gaussian spectrum
        self.amplitude = 1e3*np.exp(-(self.kk-k0)**2/(2*dk**2))
        self.forc = np.zeros(self.sizevar)

    def add_forcing(self, x, t, dxdt):
        """ add the forcing term on x[0]=the vorticity """
        dt = t-self.t0
        self.t0 = t
        gamma = dt/self.tau
        nh = self.nh

        # white noise phase
        phase = np.random.normal(size=(self.ny, self.nx))*2*np.pi
        hnoise = self.amplitude*np.exp(1j*phase)

        forc = np.real(np.fft.ifft2(hnoise))
        # time filter to introduce some time correlation in the forcing
        self.forc[nh:-nh, nh:-nh] = gamma*forc + \
            self.forc[nh:-nh, nh:-nh]*(1-gamma)
        self.fill_halo(self.forc)
        # impose zero net vorticity (otherwise the solver does not converge)
        self.forc -= self.domain_integration(self.forc) / self.area

        # add the forcing to the vorticity tendency
        dxdt[0] += self.forc


grid = Grid(param)
param.Kdiff = 5e-4*grid.dx

f2d = Fluid2d(param, grid)
model = f2d.model
model.forc = Forcing(param, grid)

xr, yr = grid.xr, grid.yr
vor = model.var.get('vorticity')
trac = model.var.get('tracer')

# set an initial small scale random vorticity field (white noise)
np.random.seed(42)


def set_x_and_k(n, L):
    k = ((n//2+np.arange(n)) % n) - n//2
    return (np.arange(n)+0.5)*L/n, 2*np.pi*k/L


_, kx = set_x_and_k(param.nx, np.pi)
_, ky = set_x_and_k(param.ny, np.pi)
kkx, kky = np.meshgrid(kx, ky)
kk = np.sqrt(kkx**2 + kky**2)

k0 = param.k0
dk = 1

phase = np.random.normal(size=(param.ny, param.nx))*2*np.pi

hnoise = np.exp(-(kk-k0)**2/(2*dk**2))*np.exp(1j*phase)

noise = np.zeros_like(vor)
nh = grid.nh
noise[nh:-nh, nh:-nh] = 1e3*np.real(np.fft.ifft2(hnoise))

grid.fill_halo(noise)
noise -= grid.domain_integration(noise)*grid.msk/grid.area
vor[:] = noise*grid.msk

# initialize the passive tracer with square tiles (better idea?)
trac[:] = np.round(xr*6) % 2 + np.round(yr*6) % 2


model.set_psi_from_vorticity()
model.diagnostics(model.var, 0)
energy = model.diags['ke']
print('energy=', energy)
vor = vor*param.Lx/np.sqrt(energy)
vor = vor/1e3
model.set_psi_from_vorticity()
f2d.loop()
