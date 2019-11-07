import numpy as np
from operators import Operators
from variables import Var
from timescheme import Timescheme
from importlib import import_module
import fortran_diag as fd
import fortran_advection as fa

from timers import Timers


class Euler(object):
    """ Euler model

    It provides the step(t,dt) function
    and 'var' containing the mode state
    """

    def __init__(self, param, grid):

        self.list_param = ['forcing', 'noslip', 'timestepping',
                           'diffusion', 'Kdiff',
                           'forcing_module', 'additional_tracer',
                           'enforce_momentum', 'var_to_save',
                           'customized', 'custom_module', 'spongelayer']
        param.copy(self, self.list_param)

        # for diagnostics
        self.list_param = ['yr', 'nh', 'msk', 'area', 'mpitools',
                           'dx', 'xr', 'yr', 'r2',
                           'x0', 'y0', 'x2', 'y2',
                           'isisland', 'Lx', 'ny']
        grid.copy(self, self.list_param)

        # for variables
        param.varname_list = ['vorticity', 'psi', 'u', 'v',
                              'source']  # ,'tauw','wshear']
        param.tracer_list = ['vorticity']
        param.whosetspsi = ('vorticity')

        if 'tauw' in self.var_to_save:
            param.varname_list.append('tauw')

        if 'wshear' in self.var_to_save:
            param.varname_list.append('wshear')

        if hasattr(self, 'additional_tracer'):
            for k in range(len(self.additional_tracer)):
                trac = self.additional_tracer[k]
                param.varname_list.append(trac)
                param.tracer_list.append(trac)

        param.sizevar = [grid.nyl, grid.nxl]
        self.var = Var(param)
        self.work = np.zeros(param.sizevar)

        # for timing the code
        self.timers = Timers(param)

        # for operators
        self.ope = Operators(param, grid)

        # for timescheme
        self.tscheme = Timescheme(param, self.var.state)
        self.tscheme.set(self.dynamics, self.timestepping)

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
                    exit(0)
                self.forc = f.Forcing(param, grid)

        if self.spongelayer:
            # sponge layer zone [0 = full sponge, 1 = no sponge]
            self.spongemsk = (1-(1+np.tanh((self.xr - self.Lx)/0.1))*0.5)
            
        self.diags = {}

        if self.customized:
            try:
                f = import_module(self.custom_module)
                self.extrastep = f.Step(param, grid)
            except ImportError:
                print('module %s for forcing cannot be found'
                      % self.custom_module)
                print('make sure file **%s.py** exists'
                      % self.custom_module)
                exit(0)

    def step(self, t, dt):

        # 1/ integrate dynamics
        self.tscheme.forward(self.var.state, t, dt)

        # 2/ integrate source
        if self.noslip:
            self.add_noslip(self.var.state)
            source = self.var.get('source')
            source /= dt
            if 'tauw' in self.var_to_save:
                tauw = self.var.get('tauw')
                tauw[:, :] = source*self.dx*self.dx

            if 'wshear' in self.var_to_save:
                wshear = self.var.get('wshear')
                self.ope.wallshear(self.var.state, self.work)
                wshear[:, :] = self.work

        if self.customized:
            self.extrastep.do(self.var, t, dt)

        # # 3/ dye tracer
        if 'dye' in self.var.varname_list:
            i, jp, jm = 1, self.ny//2+3, self.ny//2-3

            dye = self.var.get('dye')
            dye[jp+self.nh, i] = 1
            dye[jm+self.nh, i] = -1

        if 'age' in self.var.varname_list:
            age = self.var.get('age')
            age += dt*self.msk
            age[:, self.nh] = 0.

        if self.spongelayer:
            # 4/ sponge layer
            w = self.var.get('vorticity')
            w *= self.spongemsk
            if 'dye' in self.var.varname_list:
                dye *= self.spongemsk
            if 'age' in self.var.varname_list:
                age *= self.spongemsk

        self.set_psi_from_vorticity()

    def dynamics(self, x, t, dxdt):
        """Compute the tendency terms + invert the streamfunction"""

        self.timers.tic('rhs_adv')
        self.ope.rhs_adv(x, t, dxdt)
        self.timers.toc('rhs_adv')

        if (self.tscheme.kstage == self.tscheme.kforcing):
            if self.forcing:
                self.forc.add_forcing(x, t, dxdt)

            if self.diffusion:
                self.ope.rhs_diffusion(x, t, dxdt)

            if self.diffusion or self.forcing:
                self.ope.invert_vorticity(dxdt, flag='fast')
        else:
            self.timers.tic('invert')
            self.ope.invert_vorticity(dxdt, flag='fast')
            #self.ope.invert_vorticity(dxdt)
            self.timers.toc('invert')

    def add_noslip(self, x):
        self.timers.tic('noslip')
        source = self.var.get('source')
        self.ope.rhs_noslip(x, source)
        self.timers.toc('noslip')

        # if not(self.enforce_momentum):
        self.timers.tic('invert')
        self.ope.invert_vorticity(x, flag='fast', island=self.isisland)
        self.timers.toc('invert')

    def set_psi_from_vorticity(self):
        self.ope.invert_vorticity(self.var.state, island=self.isisland)

    def diagnostics(self, var, t):
        """ should provide at least 'maxspeed' (for cfl determination) """
        self.timers.tic('diag')
        u = var.get('u')
        v = var.get('v')
        trac = var.get('vorticity')
        psi = var.get('psi')
        source = self.var.get('source')

        xr = self.xr
        yr = self.yr

        ke, maxu = fd.computekemaxu(self.msk, u, v, self.nh)
        # ke = fd.computekewithpsi(self.msk, trac, psi, self.nh)

        z, z2 = fd.computesumandnorm(self.msk, trac, self.nh)

        px = fd.computedotprod(self.msk, trac, xr, self.nh)
        py = fd.computedotprod(self.msk, trac, yr, self.nh)

        angmom = fd.computesum(self.msk, psi, self.nh)
        sce = fd.computedotprod(self.msk, trac, source, self.nh)

        cst = self.mpitools.local_to_global([(maxu, 'max'), (ke, 'sum'),
                                             (z, 'sum'), (z2, 'sum'),
                                             (px, 'sum'), (py, 'sum'),
                                             (angmom, 'sum'),
                                             (sce, 'sum')])

        self.diags['maxspeed'] = cst[0]
        self.diags['ke'] = cst[1] / self.area
        self.diags['vorticity'] = cst[2] / self.area
        self.diags['enstrophy'] = 0.5*cst[3] / self.area
        self.diags['px'] = cst[4] / self.area
        self.diags['py'] = cst[5] / self.area
        self.diags['angmom'] = cst[6] / self.area
        self.diags['source'] = cst[7] / self.area

        self.timers.toc('diag')
