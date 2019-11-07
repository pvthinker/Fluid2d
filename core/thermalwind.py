from param import Param
from operators import Operators
from variables import Var
from timescheme import Timescheme
from importlib import import_module
import fortran_diag as fd
from numpy import sqrt, sum, zeros, tanh


class Thermalwind(Param):
    """ Thermal wind model

    It provides the step(t,dt) function
    and 'var' containing the mode state
    """

    def __init__(self, param, grid):

        self.list_param = ['forcing', 'noslip', 'timestepping',
                           'forcing_module', 'additional_tracer',
                           'myrank',
                           'gravity', 'diffusion', 'Kdiff', 'f0']
        param.copy(self, self.list_param)

        # for potential energy
        self.list_param = ['xr', 'yr', 'nh',
                           'Lx', 'msk', 'area', 'mpitools', 'dx']
        grid.copy(self, self.list_param)

        # for variables
        param.varname_list = ['vorticity', 'psi',
                              'u', 'v', 'buoyancy', 'V', 'qE']
        param.tracer_list = ['vorticity', 'buoyancy', 'V']
        param.whosetspsi = ('vorticity')

        if hasattr(self, 'additional_tracer'):
            for k in range(len(self.additional_tracer)):
                trac = self.additional_tracer[k]
                param.varname_list.append(trac)
                param.tracer_list.append(trac)

        param.sizevar = [grid.nyl, grid.nxl]
        self.var = Var(param)

        # for operators
        self.ope = Operators(param, grid)

        # for timescheme
        self.tscheme = Timescheme(param, self.var.state)

        if self.forcing:
            if self.forcing_module == 'embedded':
                print('Warning: check that you have indeed added the forcing to the model')
                print('Right below the line    : model = f2d.model')
                print('you should have the line: model.forc = Forcing(param, grid)')

                pass
            else:
                try:
                    f = import_module(self.forcing_module)

                except ImportError:
                    print('module %s for forcing cannot be found'
                          % self.forcing_module)
                    print('make sure file **%s.py** exists' % self.forcing_module)
                    sys.exit(0)

                self.forc = f.Forcing(param, grid)


        self.diags = {}

    def step(self, t, dt):

        # 1/ integrate advection
        self.tscheme.set(self.dynamics, self.timestepping)
        self.tscheme.forward(self.var.state, t, dt)

        # 2/ integrate source
        if self.noslip:
            self.add_noslip(self.var.state)

        self.set_psi_from_vorticity()
        self.compute_pv()

    def compute_pv(self):
        qE = self.var.get('qE')
        v = self.var.get('V')
        b = self.var.get('buoyancy')

        # Ertel PV
        y = qE*0.
        y[1:-1, 1:-1] = self.ope.jacobian(v+self.f0*self.xr, b)
        y *= self.msk
        self.ope.fill_halo(y)
        qE[:, :] = y

    def dynamics(self, x, t, dxdt):
        self.ope.rhs_adv(x, t, dxdt)
        self.ope.rhs_thermalwind(x, t, dxdt)  # add the r.h.s. terms

        if (self.tscheme.kstage == self.tscheme.kforcing):
            if self.forcing:
                self.forc.add_forcing(x, t, dxdt)

            if self.diffusion:
                self.ope.rhs_diffusion(x, t, dxdt)

        else:
            self.ope.invert_vorticity(dxdt, flag='fast')

    def sources(self, x, t, dxdt):
        if self.noslip:
            self.ope.rhs_noslip(x, t, dxdt)
            #self.ope.invert_vorticity(x, flag='full')

    def set_psi_from_vorticity(self):
        self.ope.invert_vorticity(self.var.state)

    def diagnostics(self, var, t):
        """ should provide at least 'maxspeed' (for cfl determination) """

        nh = self.nh
        u = var.get('u')
        v = var.get('v')
        vort = var.get('vorticity')
        buoy = var.get('buoyancy')
        V = var.get('V')
        qE = var.get('qE')
        qEneg = qE.copy()
        qEneg[qE>0] = 0.
        
        ke, maxu = fd.computekemaxu(self.msk, u, v, self.nh)

        z, z2 = fd.computesumandnorm(self.msk, vort, self.nh)

        b, b2 = fd.computesumandnorm(self.msk, buoy, self.nh)

        vm, v2 = fd.computesumandnorm(self.msk, V, self.nh)

        q, q2 = fd.computesumandnorm(self.msk, qE, self.nh)
        qn, qn2 = fd.computesumandnorm(self.msk, qEneg, self.nh)

        # potential energy (minus sign because buoyancy is minus density)
        pe = fd.computesum(self.msk, buoy*self.yr, nh)
        pe = - self.gravity*pe

        cst = self.mpitools.local_to_global([
            (maxu, 'max'), (ke, 'sum'),
            (z, 'sum'), (z2, 'sum'),
            (pe, 'sum'), (b, 'sum'),
            (b2, 'sum'), (q, 'sum'),
            (q2, 'sum'), (qn, 'sum'),
            (qn2, 'sum'), (v2, 'sum')])

        a = self.area
        self.diags['maxspeed'] = cst[0]
        self.diags['ke'] = (cst[1]) / a
        self.diags['keV'] = (0.5*cst[11])/a
        self.diags['pe'] = cst[4] / a
        self.diags['energy'] = self.diags['ke'] + self.diags['pe'] + self.diags['keV']

        self.diags['vorticity'] = cst[2] / a
        self.diags['enstrophy'] = 0.5*cst[3] / a

        bm = cst[5] / a
        self.diags['buoyancy'] = bm
        self.diags['brms'] = sqrt(cst[6] / a - bm**2)

        pvm = cst[7] / a
        pvneg_mean = cst[9] / a
        self.diags['pv_mean'] = pvm
        self.diags['pv_std'] = sqrt(cst[8] / a - pvm**2)
        self.diags['pvneg_mean'] = pvneg_mean
        self.diags['pvneg_std'] = sqrt(cst[10] / a - pvneg_mean**2)
