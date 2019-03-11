from operators import Operators
from variables import Var
from timescheme import Timescheme
import numpy as np
from importlib import import_module
import fortran_diag as fd


class QG2L(object):
    """Two-Layers Quasi-geostrophic model

    The prognostic variable is 'pv'. It is the full PV, which includes
    the planetary background PV 'pvback'. Full PV is materially
    conserved unlike the PV anamoly that has a beta*v source term

    """
    def __init__(self, param, grid):
        self.list_param = ['forcing', 'diffusion', 'Kdiff', 'noslip',
                           'timestepping', 'beta', 'Rd', 'ageostrophic',
                           'forcing_module']
        param.copy(self, self.list_param)

        # for diagnostics
        self.list_param = ['yr', 'nh', 'msk', 'area', 'mpitools']
        grid.copy(self, self.list_param)

        # for variables
        param.varname_list = ('pv', 'psi', 'u', 'v',
                              'pvanom', 'vorticity')

        if param.ageostrophic:
            param.varname_list += ('ua', 'va')
            print(param.varname_list)

        param.sizevar = [2, grid.nyl, grid.nxl]
        self.var = Var(param)
        self.source = np.zeros(param.sizevar)

        self.ipva = self.var.varname_list.index('pvanom')
        self.ipv = self.var.varname_list.index('pv')
        self.ivor = self.var.varname_list.index('vorticity')
        self.ipsi = self.var.varname_list.index('psi')

        # background pv
        self.pvback = self.beta*(grid.yr-grid.Ly*.5)*grid.msk

        # for operators
        param.tracer_list = ['pv']
        if param.ageostrophic:
            param.tracer_list += ['ua', 'va']

        param.whosetspsi = ('pvanom')
        param.qgoperator = True
        self.ope = Operators(param, grid)

        # for timescheme
        self.tscheme = Timescheme(param, self.var.state)
        self.dx0 = self.tscheme.dx0

        self.kt = 0

        if self.forcing:

            try:
                f = import_module(self.forcing_module)

            except:
                print('module %s for forcing cannot be found'
                      % self.forcing_module)
                print('make sure file **%s.py** exists' % self.forcing_module)
                exit(0)

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

        # 2/ integrate source
        if self.noslip:
            self.add_noslip(self.var.state)

        # 3/ diagnostic fields
        self.var.state[self.ipva] = self.var.state[self.ipv] - self.pvback
        self.var.state[self.ivor] = (self.var.state[self.ipva]
                                     - (self.Rd**-2)*self.var.state[self.ipsi])

        if self.ageostrophic:
            """ diagnose the ageostrophic velocity (factor 1/f missing)

            du/dt + J(psi, u) = -f va
            dv/dt + J(psi, v) = +f ua

            implementation:

            - store u in 'ua'
            - step ua
            - compare update u (real u diagnosed from psi) with updated 'ua'
            - the difference is the ageostrophic velocity

            Yet to be done: implement the proper staggering, 'u' and 'v' sit
            on edges while the advection scheme assumes cell centers quantities
            """
            ua = -(self.var.state[iv] - self.var.state[iva])/dt
            va = +(self.var.state[iu] - self.var.state[iua])/dt
            self.var.state[iva][:, :] = va.copy()
            self.var.state[iua][:, :] = ua.copy()

    def dynamics(self, x, t, dxdt):
        """Compute the tendency terms + invert the streamfunction"""

        dxdt[:] = 0.
        self.ope.rhs_adv_multilayers(x, t, dxdt)
        if (self.tscheme.kstage == self.tscheme.kforcing):
            if self.forcing:
                self.forc.add_forcing(x, t, dxdt)
            if self.diffusion:
                self.ope.rhs_diffusion(x, t, dxdt)
        # substract the background pv
        dxdt[self.ipva][:, :] = dxdt[self.ipv]  # -self.pvback*self.dt
#        self.ope.invert_vorticity(dxdt,flag='fast')
        self.ope.invert_vorticity(dxdt, flag='fast')
#        self.kt +=1

    def add_noslip(self, x):
        """Add the no-slip term """

        self.ope.rhs_noslip(x, self.source)
        self.ope.invert_vorticity(x, flag='fast', island=self.isisland)

    def add_backgroundpv(self):
        self.var.state[self.ipv] += self.pvback

    def set_psi_from_pv(self):
        """Solve the elliptic equation for the streamfunction

        The unknown of the Helmholtz equation is the PV anomaly,
        not the full PV"""

        state = self.var.state
        state[self.ipva] = state[self.ipv] - self.pvback
        self.ope.invert_vorticity(self.var.state, flag='full')

    def diagnostics(self, var, t):
        """ Integral diagnostics for the QG model

        should provide at least 'maxspeed' (for cfl determination) """

        u = var.get('u')
        v = var.get('v')
        trac = var.get('pv')

        ke, maxu = fd.computekemaxu(self.msk, u, v, self.nh)

        z, z2 = fd.computesumandnorm(self.msk, trac, self.nh)

        cst = self.mpitools.local_to_global([(maxu, 'max'),
                                             (ke, 'sum'),
                                             (z, 'sum'),
                                             (z2, 'sum')])

        self.diags['maxspeed'] = cst[0]
        self.diags['ke'] = cst[1] / self.area
        self.diags['pv'] = cst[2] / self.area
        self.diags['pv2'] = cst[3] / self.area
