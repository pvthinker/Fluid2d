from __future__ import print_function
from output import Output
import numpy as np
from importlib import import_module
import signal
from time import time as clock
import sys
from subprocess import call
import os
try:
    import fluxes as Flx
except:
    print("[ERROR] unable to import compiled module")
    print("[INFO]  you likely forgot to compile the Fortran modules")
    print("[INFO]  go in the main Fluid2d folder")
    print("> make")
    sys.exit()

class Fluid2d(object):
    def __init__(self, param, grid):

        # let's first check that no param is obviously incorrect
        param.checkall()
        
        # copy the launch script into 'expname.py' to allow
        # for experiment reproducibility
        launchscript = sys.argv[0]
        param.datadir = param.datadir.replace('~', os.getenv("HOME"))
        param.expdir = '%s/%s' % (param.datadir, param.expname)
        if param.myrank == 0:
            if os.path.isdir(param.expdir):
                pass
            else:
                os.makedirs(param.expdir)
            savedscript = '%s/%s.py' % (param.expdir, param.expname)
            outfile = '%s/output.txt' % param.expdir
            if os.path.exists(outfile):
                print('Warning: this experiment has already been ran, output.txt already exists')
                print('dummy.txt will be used instead')
                outfile = '%s/dummy.txt' % param.expdir
            
            sys.stdout = Logger(outfile)

            if os.path.exists(savedscript):
                print('Warning: the python script already exists in %s' % param.expdir)
                print('the script won''t be copied')
                pass
            else:
                self.savedscript = savedscript
                call(['cp', launchscript, savedscript])

        self.list_param = ['modelname', 'tend', 'adaptable_dt', 'dt', 'cfl',
                           'dtmax', 'myrank', 'nprint', 'exacthistime',
                           'rescaledtime', 'noslip', 'geometry',
                           'diag_fluxes', 'print_param',
                           'enforce_momentum', 'forcing', 'decay',
                           'plotting_module', 'freq_save', 'freq_plot',
                           'plot_interactive', 'nbproc',
                           'isisland', 'npx', 'npy', 'nx', 'ny']

        param.copy(self, self.list_param)

        self.dt0 = self.dt

        grid.finalize_msk()
        # print('momentum=',self.enforce_momentum)
        self.list_grid = ['dx', 'dy', 'nh', 'msk', 'xr0', 'yr0', 'x2', 'y2']
        grid.copy(self, self.list_grid)

        if param.modelname == 'euler':
            if self.geometry not in ['closed', 'disc']:
                self.enforce_momentum = False
            from euler import Euler
            self.model = Euler(param, grid)
        else:
            # not yet implemented in other models
            self.enforce_momentum = False

        if param.modelname == 'advection':
            from advection import Advection
            self.model = Advection(param, grid)

        if param.modelname == 'boussinesq':
            from boussinesq import Boussinesq
            self.model = Boussinesq(param, grid)

        if param.modelname == 'boussinesqTS':
            from boussinesqTS import BoussinesqTS
            self.model = BoussinesqTS(param, grid)

        if param.modelname == 'quasigeostrophic':
            from quasigeostrophic import QG
            self.model = QG(param, grid)

        if param.modelname == 'sqg':
            from sqg import SQG
            self.model = SQG(param, grid)

        if param.modelname == 'thermalwind':
            from thermalwind import Thermalwind
            self.model = Thermalwind(param, grid)

        if self.modelname == 'quasigeostrophic':
            self.enstrophyname = 'pv2'
        elif self.modelname == 'sqg':
            self.enstrophyname = 'pv2'
        else:
            self.enstrophyname = 'enstrophy'

        if self.isisland:
            grid.island.finalize(self.model.ope.mskp)
            self.model.ope.rhsp = grid.island.rhsp
            self.model.ope.psi = grid.island.psi
            # self.diag = Diag(param,grid)

        if self.diag_fluxes:
            self.flx = Flx.Fluxes(param, grid, self.model.ope)
            flxlist = self.flx.fullflx_list
        else:
            flxlist = None
            
        if self.plot_interactive:
            try:
                p = import_module(self.plotting_module)
                # print(self.plotting_module)
            except ImportError:
                print('problem with the interactive plotting')
                print('this might be due to a backend issue')
                print('try to rerun the code with')
                print('param.plot_interactive = False')
                exit(0)

            self.plotting = p.Plotting(param, grid,
                                       self.model.var,
                                       self.model.diags)

        self.tracer_list = param.tracer_list
        # here's a shortcut to the model state
        self.state = self.model.var.state

        self.t = 0.
        self.kt = 0

        self.output = Output(param, grid, self.model.diags, flxlist=flxlist)
        self.print_config(param, start=True)

    def print_config(self, param, start=True):
        outputfiles = [self.output.hisfile, self.output.diagfile]
        if self.diag_fluxes:
            outputfiles += [self.output.flxfile]
        if self.plot_interactive:
            if hasattr(self.plotting, 'mp4file'):
                outputfiles += [self.plotting.mp4file]

        if hasattr(self, 'savedscript'):
            outputfiles += [self.savedscript]

        if self.myrank == 0:
            if start:
                print('-'*50)
                print(' Fluid2d summary:')
                print('-'*50)
                print('  - model equations: %s' % self.modelname)
                print('  - grid size: %i x %i' % (self.nx, self.ny))
                print('  - integration time: %.2f' % self.tend)
                print('  - advection schemes applied to:')
                for trac in self.tracer_list:
                    print('    - %s' % trac)
            else:
                print(' Output files:')
                print('-'*50)
                for f in outputfiles:
                    print('  - %s' % f)
                print('-'*50)
                print(' You may recover all the experiment parameters by doing')
                print(' ncdump -h %s' % self.output.hisfile)
                print('')
                print(' or browse the history file by doing')
                print(' ncview %s' % self.output.hisfile)
                print('-'*50)                                

        if self.print_param and (self.myrank == 0) and start:
            print('-'*50)
            print(' Extended list of all parameters')
            print('-'*50)
            param.printvalues()

    def loop(self, joinhis=True, keepplotalive=False):
        """Fluid2d time loop"""
        if self.myrank == 0:
            print('-'*50)
            print(' Starting the time loop')
            if self.decay and (self.modelname == 'euler'):
                print('   in this simulation the kinetic energy is expected')
                print('   to be decreasing, the code will issue a warning msg')
                print('   during the integration if ke increases')
            print('-'*50)

        self.model.diagnostics(self.model.var, self.t)
        self.model.diags['dkedt'] = 0.
        self.model.diags['dvdt'] = 0.
        data = {'his': self.model.var.state,
                'diag': self.model.diags}
        if self.diag_fluxes:
            data['flx'] = self.flx.flx
            # it is important to set dt (in the case of 'adaptable_dt')
            # because otherwise this first diagnostic can go crazy
            # and for an unknown reason this can spoil the whole file
            self.set_dt(self.kt)
            self.flx.diag_fluxes(self.model.var.state, self.t, self.dt)

        self.output.do(data, self.t, self.kt)

        if self.plot_interactive:
            if hasattr(self.plotting, 'fig'):
                pass
            else:
                self.plotting.create_fig(self.t)

        nh = self.nh
        self.z2d = self.model.var.state[0][nh:-nh, nh:-nh]

        def signal_handler(signal, frame):
            if self.myrank == 0:
                print('\n hit ctrl-C, stopping', end='')
            self.stop = True

        signal.signal(signal.SIGINT, signal_handler)

        self.stop = False

        kt0 = 0
        t0 = clock()
        reduce = 0
        while (self.t < self.tend and not(self.stop)):

            self.set_dt(self.kt)

            # reduce the time step to arrive right at the desired
            # next history time when param.exacthistime == True (default)
            # this occurs when param.adaptable_dt == True
            #
            # to avoid having suddenly a very small time step, which breaks
            # the smoothness of the time integration, the time step is
            # adjusted 8 time steps ahead of time.
            if self.exacthistime and self.adaptable_dt:
                if (self.t+8*self.dt > self.output.tnexthis) and (reduce == 0):
                    reduce = 8
                if (reduce > 0):
                    self.dt = (self.output.tnexthis-self.t)/(reduce*0.95)
                    reduce -= 1
                if (self.t+self.dt > self.output.tnexthis):
                    dt = self.dt
                    reduce = 0
                    self.dt = self.output.tnexthis-self.t

            self.model.step(self.t, self.dt)

            if self.rescaledtime == 'enstrophy':
                self.t += self.dt * np.sqrt(self.model.diags['enstrophy'])

            else:
                self.t += self.dt

            self.kt += 1
            ke_old = self.model.diags['ke']
            ens_old = self.model.diags[self.enstrophyname]
            self.model.diagnostics(self.model.var, self.t)

            if self.enforce_momentum:
                self.enforce_zero_momentum()
                self.model.diagnostics(self.model.var, self.t)

            ke = self.model.diags['ke']
            ens = self.model.diags[self.enstrophyname]
            dkedt = (ke-ke_old)/self.dt
            dvdt = (ens-ens_old)/self.dt
            self.model.diags['dkedt'] = dkedt
            self.model.diags['dvdt'] = dvdt

            if ((ke > ke_old) and (self.myrank == 0)
                    and (self.decay) and (self.modelname == 'euler')):
                dke = (ke-ke_old)/ke

                warning = '\033[0;32;40m' + 'WARNING dlog(ke)'
                white = '\033[0m'
                print('\rkt=%-4i %s = %.2g%s' %
                      (self.kt, warning, dke, white), end='')

            if self.diag_fluxes and (self.t >= self.output.tnexthis):
                # this is a costly diagnostic, do it only before writing it
                self.flx.diag_fluxes(self.model.var.state, self.t, self.dt)
            
            self.output.do(data, self.t, self.kt)

            if self.dt == self.dtmax:
                flag = '*'
            else:
                flag = ''

            if (self.myrank == 0) and (self.kt % self.nprint == 0)\
               or (self.t >= self.tend):
                print('\rkt=%-4i / t=%-7.3f %s / dt=%-7.3f ' %
                      (self.kt, self.t, flag, self.dt), end='')

            if (self.kt % self.freq_plot == 0) and self.plot_interactive:
                self.plotting.update_fig(self.t, self.dt, self.kt)

            # check for blow-up
            if self.model.diags['maxspeed'] > 1e3:
                self.stop = True
                if self.myrank == 0:
                    print()
                    print('max|u| > 1000, blow-up detected, stopping')

        if self.myrank == 0:
            print('\ndone')

        if self.plot_interactive and not(keepplotalive):
            self.plotting.finalize()

        if hasattr(self.model, 'timers'):
            if (self.myrank == 0):
                print('-'*50)
                print(' A few model performances metrics')
                print('-'*50)
                self.model.timers._print()

            if self.myrank == 0:
                dt = clock()-t0
                nkt = (self.kt-kt0)
                gamma = dt*self.npx*self.npy/(nkt*self.nx*self.ny)
                print()
                print('  - Wall  time      : %f s' % dt)
                print('  - Nb of iterations: %i' % nkt)
                print('  - Time per ite    : %5.3f s' % (dt/nkt))
                print('  - Rescaled time   : %5.3e s (per ite, per dof)' % gamma)
                print('-'*50)

        if (self.myrank == 0) and (joinhis):
            self.output.dump_diag()

            if self.nbproc <= 64:
                self.output.join()

            else:
                print('too many cores will be idle ' +
                      'join the history files manually')
        self.print_config(None, start=False)

    def enforce_zero_momentum(self):
        if self.enforce_momentum:
            vor = self.model.var.get('vorticity')
            px = self.model.diags['px']/self.x2
            py = self.model.diags['py']/self.y2
            vor[:] -= (self.xr0 * px + self.yr0 * py)
            x = self.model.var.state
            self.model.ope.invert_vorticity(x, flag='fast')

    def relaunch(self, duration=0):
        t, dt, kt, tnextdiag, tnexthis = self.restart.read(self.model.var)
        self.t = t
        self.dt = dt
        self.kt = kt
        self.output.tnextdiag = tnextdiag
        self.output.tnexthis = tnexthis

        self.tend += duration

        print('t = %f / kt = %i / nextdiag = %f / nexthist =%f' %
              (t, kt, tnextdiag, tnexthis))

        self.model.set_psi_from_vorticity()
        self.loop()

    def set_dt(self, kt):

        if ((self.adaptable_dt) & (self.model.diags['maxspeed'] != 0)):
            dt = self.cfl * min(self.dx, self.dy) / \
                self.model.diags['maxspeed']
            # filter in time, change the filter_coef
            # if dt > 0.8*self.dt:
            #     self.filter_coef = 0.05
            # else:
            #     # time step decreases too fast
            #     # don't filter, otherwise the model blows up
            #     self.filter_coef = 1.
            # self.dt = (1.-self.filter_coef)*self.dt + self.filter_coef*dt
            self.dt = dt
            if self.dt > self.dtmax:
                self.dt = self.dtmax

        else:
            self.dt = self.dt0

        # transfer this information to the advection scheme
        self.model.ope.cst[3] = self.model.diags['maxspeed']


class Logger(object):
    def __init__(self, logfile):
        self.terminal = sys.stdout
        self.log = open(logfile, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass
