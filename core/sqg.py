from operators import Operators
from variables import Var
from timescheme import Timescheme
import numpy as np
from importlib import import_module
import fortran_diag as fd
import fortran_operators as fo
import sys


class SQG(object):
    """ Surface Quasi-Geostrophic model

    The prognostic variable is 'pv'. It is the surface PV

    The streamfunction is computed in the Fourier space

    The geometry must be biperiodic

    This model does not support multicores (because of the fft)

    """

    def __init__(self, param, grid):
        self.list_param = ['forcing', 'diffusion', 'Kdiff',
                           'timestepping', 'ageostrophic',
                           'forcing_module', 'geometry']
        param.copy(self, self.list_param)

        # for diagnostics
        self.list_param = ['yr', 'nh', 'msk', 'area', 'mpitools']
        grid.copy(self, self.list_param)
        assert grid.geometry == "perio", "SQG imposes a biperiodic domain"

        assert param.ageostrophic == False, "Ageostrophic velocity not yet tested"

        # for variables
        param.varname_list = ['pv', 'psi', 'u', 'v', 'vorticity']

        if param.ageostrophic:
            # TODO: not yet tested for SQG
            param.varname_list += ['ua', 'va']

        param.sizevar = [grid.nyl, grid.nxl]
        self.var = Var(param)
        self.source = np.zeros(param.sizevar)

        self.ipv = self.var.varname_list.index('pv')
        self.ivor = self.var.varname_list.index('vorticity')
        self.ipsi = self.var.varname_list.index('psi')

        # for operators
        param.tracer_list = ['pv']

        if param.ageostrophic:
            param.tracer_list += ['ua', 'va']

        param.whosetspsi = ('pv')
        param.sqgoperator = True
        self.ope = Operators(param, grid)

        # for timescheme
        self.tscheme = Timescheme(param, self.var.state)
        self.dx0 = self.tscheme.dx0

        self.kt = 0

        if self.forcing:
            if self.forcing_module == 'embedded':
                print(
                    'Warning: check that you have indeed added the forcing to the model')
                print('Right below the line    : model = f2d.model')
                print('you should have the line: model.forc = Forcing(param, grid)')

                pass
            else:
                try:
                    f = import_module(self.forcing_module)

                except ImportError:
                    print('module %s for forcing cannot be found'
                          % self.forcing_module)
                    print('make sure file **%s.py** exists' %
                          self.forcing_module)
                    sys.exit(0)

                self.forc = f.Forcing(param, grid)

        self.diags = {}

        self.tscheme.set(self.dynamics, self.timestepping)

    def step(self, t, dt):
        """Integrate the model over one time step dt"""

        if self.ageostrophic:
            iu = self.var.varname_list.index('u')
            iv = self.var.varname_list.index('v')
            iua = self.var.varname_list.index('ua')
            iva = self.var.varname_list.index('va')

            self.var.state[iua][:, :] = self.var.state[iu].copy()
            self.var.state[iva][:, :] = self.var.state[iv].copy()

        self.dt = dt
        # 1/ integrate advection
        self.tscheme.forward(self.var.state, t, dt)

        self.set_psi_from_pv()

        if self.ageostrophic:
            """ diagnose the ageostrophic velocity (factor 1/f missing)

            du/dt + J(psi, u) = +f va + (pvback)*vg
            dv/dt + J(psi, v) = -f ua - (pvback)*ug

            implementation:

            - store u in 'ua'
            - step ua
            - compare update u (real u diagnosed from psi) with updated 'ua'
            - the difference is the ageostrophic velocity

            Yet to be done: implement the proper staggering, 'u' and
            'v' sit on edges while the advection scheme assumes cell
            centers quantities  """
            ug = self.var.state[iu]
            vg = self.var.state[iv]
            ua = -(vg - self.var.state[iva])/dt - self.pvback*ug
            va = +(ug - self.var.state[iua])/dt - self.pvback*vg
            self.var.state[iva][:, :] = va.copy()
            self.var.state[iua][:, :] = ua.copy()

    def dynamics(self, x, t, dxdt):
        """Compute the tendency terms + invert the streamfunction"""

        dxdt[:] = 0.
        self.ope.rhs_adv(x, t, dxdt)
        if (self.tscheme.kstage == self.tscheme.kforcing):
            if self.forcing:
                self.forc.add_forcing(x, t, dxdt)
            if self.diffusion:
                self.ope.rhs_diffusion(x, t, dxdt)

        else:
            self.ope.fourier_invert_vorticity(dxdt, flag='fast')

    def set_psi_from_pv(self):
        """Solve the elliptic equation for the streamfunction

        The unknown of the Helmholtz equation is the PV anomaly,
        not the full PV"""

        self.ope.fourier_invert_vorticity(self.var.state, flag='full')

    def diagnostics(self, var, t):
        """Integral diagnostics for the SQG model

        should provide at least 'maxspeed' (for cfl determination)

        WARNING: the diag for KE and APE are wrong. A SQG flow is 3D
        the integration should be carried over in the Fourier space

        The total PV should be conserved (in the absence of external
        forcing) and the enstrophy should be conserved as long as the
        flow is not developping filaments (grid-scale enstrophy). The
        numerical dissipation wipes out Filaments which yields
        enstrophy dissipation.

        """

        u = var.get('u')
        v = var.get('v')
        trac = var.get('pv')
        psi = var.get('psi')

        #fo.cornertocell(psi, self.ope.work)
        #psim, psi2 = fd.computesumandnorm(self.msk, self.ope.work, self.nh)
        #ape = 0.5 * psi2 / self.Rd**2

        ke, maxu = fd.computekemaxu(self.msk, u, v, self.nh)

        z, z2 = fd.computesumandnorm(self.msk, trac, self.nh)

        cst = self.mpitools.local_to_global([(maxu, 'max'),
                                             (ke, 'sum'),
                                             (z, 'sum'),
                                             (z2, 'sum')])

        self.diags['maxspeed'] = cst[0]
        self.diags['ke'] = cst[1] / self.area
        self.diags['pv'] = cst[2] / self.area
        self.diags['pv2'] = 0.5*cst[3] / self.area
        #self.diags['ape'] = cst[4] / self.area
        #self.diags['energy'] = (cst[1]+cst[4]) / self.area
