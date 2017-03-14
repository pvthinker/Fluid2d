from __future__ import print_function
from param import Param
from output import Output
from numpy import nan,sqrt,abs
from importlib import import_module
import signal


class Fluid2d(Param):
    def __init__(self,param,grid):

        self.list_param=['modelname','tend','adaptable_dt','dt','cfl','dtmax','myrank','nprint',
                         'rescaledtime','noslip','geometry','enforce_momentum','forcing',
                         'plotting_module','freq_save','freq_plot','plot_interactive','nbproc',
                         'isisland']

        param.copy(self,self.list_param)

        grid.finalize_msk()
        #print('momentum=',self.enforce_momentum)
        self.list_grid=['dx','nh','msk','xr0','yr0','x2','y2']
        grid.copy(self,self.list_grid)


        if param.modelname=='euler':
            if not(self.geometry in ['square','disc']):
                self.enforce_momentum = False
            from euler import Euler
            self.model = Euler(param,grid)
        else:
            # not yet implemented in other models
            self.enforce_momentum = False
                

        if param.modelname=='advection':
            from advection import Advection
            self.model = Advection(param,grid)

        if param.modelname=='boussinesq':
            from boussinesq import Boussinesq
            self.model = Boussinesq(param,grid)

        if param.modelname=='quasigeostrophic':
            from quasigeostrophic import QG
            self.model = QG(param,grid)


        if self.isisland:
            grid.island.finalize(self.model.ope.mskp)
            self.model.ope.rhsp=grid.island.rhsp
            self.model.ope.psi=grid.island.psi
#        self.diag = Diag(param,grid)

        if self.plot_interactive:
            try:
                p = import_module(self.plotting_module)
            except:
                print('module %s for plotting cannot be found'%self.plotting_module)
                print('make sure file **%s.py** exists'%self.plotting_module)
                exit(0)

            self.plotting = p.Plotting(param,grid,self.model.var,self.model.diags)
        

        # here's a shortcut to the model state
        self.state = self.model.var.state
                
        self.t  = 0.
        self.kt = 0
        
        self.output=Output(param,grid,self.model.diags)



    def loop(self,joinhis=True):
        if self.myrank==0:
            print('starting the time loop')

        self.model.diagnostics(self.model.var,self.t)
        self.model.diags['dkedt'] = 0.
        self.output.do(self.model.var.state,self.model.diags,self.t,self.kt)
        
        if self.plot_interactive:
            self.plotting.create_fig(self.t)

        nh=self.nh
        self.z2d = self.model.var.state[0][nh:-nh,nh:-nh]
        
        
        def signal_handler(signal, frame):
            if self.myrank==0:
                print('\nstop requested',end='')
            self.stop=True

        signal.signal(signal.SIGINT, signal_handler)

        self.stop=False
        while (self.t<self.tend and not(self.stop)):        
            
            self.set_dt()
            self.model.step(self.t,self.dt)
            if self.rescaledtime=='enstrophy':
                self.t += self.dt * sqrt( self.model.diags['enstrophy'] )
            else:
                self.t  += self.dt
            self.kt += 1
            ke_old = self.model.diags['ke']
            self.model.diagnostics(self.model.var,self.t)


            self.enforce_zero_momentum()
            if self.enforce_momentum:
                self.model.diagnostics(self.model.var,self.t)

            ke = self.model.diags['ke']
            dkedt = -(ke -ke_old)/self.dt
            self.model.diags['dkedt'] = dkedt
#            if(ke>ke_old) and self.myrank==0:
#                print('kt = %i / ke is increasing dke/dt = %.2g'%(self.kt,-dkedt))
            
            self.output.do(self.model.var.state,self.model.diags,self.t,self.kt)
                
            flag=''
            if self.dt==self.dtmax:
                flag='*'
            if (self.myrank==0) and (self.kt%self.nprint==0):
                print('\rkt = %-4i  /  t=%-7.3f %s'%(self.kt,self.t,flag),end='')
            if (self.kt%self.freq_plot ==0) and self.plot_interactive:
                self.plotting.update_fig(self.t,self.dt,self.kt)
        if self.myrank==0:
            print('\ndone')
            
        if self.plot_interactive:
            self.plotting.finalize()


        if (self.myrank==0) and (joinhis):
            self.output.dump_diag()
            if self.nbproc<=8:
                self.output.join()
            else:
                print('to many cores will be idle, join the history files manually')

    def enforce_zero_momentum(self):
        if self.enforce_momentum:
            vor = self.model.var.get('vorticity')
            px = self.model.diags['px']/self.x2
            py = self.model.diags['py']/self.y2
            vor[:] -= ( self.xr0 * px + self.yr0 * py)
            x = self.model.var.state
            self.model.ope.invert_vorticity(x,flag='fast')
        

    def relaunch(self,duration=0):
        t,dt,kt,tnextdiag,tnexthis = self.restart.read(self.model.var)
        self.t = t
        self.dt = dt
        self.kt = kt
        self.output.tnextdiag = tnextdiag
        self.output.tnexthis  = tnexthis
               
        self.tend += duration

        print('t = %f / kt = %i / nextdiag = %f / nexthist =%f'%(t,kt,tnextdiag,tnexthis))

        self.model.set_psi_from_vorticity()
        self.loop()


    def set_dt(self):

        if ((self.adaptable_dt) & (self.model.diags['maxspeed'] != 0)):
            dt = self.cfl * self.dx / self.model.diags['maxspeed']
            # filter in time, change the filter_coef
            self.filter_coef=1.
            self.dt = (1.-self.filter_coef)*self.dt + self.filter_coef*dt
            if self.dt>self.dtmax:
                self.dt=self.dtmax
        else:
            pass

        # transfer this information to the advection scheme
        self.model.ope.cst[2]=self.model.diags['maxspeed']



