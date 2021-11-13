from param import Param
from grid import Grid
from fluid2d import Fluid2d
import numpy as np

param = Param('default.xml')
param.modelname = 'euler'
param.expname = 'fd_00'

# domain and resolution
param.nx = 64*2
param.ny = param.nx
param.Ly = param.Lx
param.npx = 1
param.npy = 1
param.geometry = 'perio'

# time
param.tend = 1200.
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
param.freq_his = 10
param.freq_diag = 1
param.exacthistime = False

# plot
param.plot_interactive = True
param.plot_var = 'vorticity'
param.plot_psi = False
param.cmap = 'RdYlBu_r'  # 'Spectral'
# change the color scale if you change the magnitude of the  vorticity
param.cax = np.array([-1, 1])*.25
param.colorscheme = 'imposed'
param.generate_mp4 = True

# physics
param.forcing = False
param.noslip = False
param.diffusion = False
#
# add a passive tracer
param.additional_tracer = ['tracer']

grid = Grid(param)
param.Kdiff = 5e-4*grid.dx

f2d = Fluid2d(param, grid)
model = f2d.model


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

# k0 is the wave-number of the spectral peak
k0 = param.nx*0.48
dk = 1

phase = np.random.normal(size=(param.ny, param.nx))*2*np.pi
hnoise = np.exp(-(kk-k0)**2/(2*dk))*np.exp(1j*phase)

noise = np.zeros_like(vor)
nh = grid.nh
noise[nh:-nh, nh:-nh] = 1e3*np.real(np.fft.ifft2(hnoise))


grid.fill_halo(noise)

magnitude = 1.  # change the magnitude here
vor[:] = magnitude * noise  # also change param.cax, the color axis

# modifiy slightly the vorticity, on just one grid point!
# vor[10, 10] += 1.

# initialize the passive tracer with square tiles
# better idea? => do it, but watch out the syntax
# trac[:] = blabla (a 2D array)
trac[:] = np.round(xr*6) % 2 + np.round(yr*6) % 2


model.set_psi_from_vorticity()
model.diagnostics(model.var, 0)
energy = model.diags['ke']
print('energy=', energy)
# uncomment below to rescale the vorticity with an intrinsic time scale
# vor = vor*param.Lx/np.sqrt(energy)
model.set_psi_from_vorticity()
f2d.loop()
