from param import Param
from operators import Operators
from variables import Var
from timescheme import Timescheme
from numpy import zeros,sqrt,sum
from importlib import import_module
from fortran_diag import *


class QG(Param):
    """ Quasi-geostrophic model 

    It provides the step(t,dt) function
    and 'var' containing the mode state
    """
    def __init__(self,param,grid):
        self.list_param=['forcing','diffusion','Kdiff','noslip','timestepping','beta','Rd',
                         'forcing_module']
        param.copy(self,self.list_param)

        # for diagnostics
        self.list_param=['yr','nh','msk','area','mpitools']
        grid.copy(self,self.list_param)

        # for variables
        param.varname_list=('pv','psi','u','v','pvanom','vorticity')
        param.sizevar     =[grid.nyl,grid.nxl]
        self.var = Var(param)
        self.source = zeros(param.sizevar)

        self.ipva = self.var.varname_list.index('pvanom')
        self.ipv  = self.var.varname_list.index('pv')
        self.ivor = self.var.varname_list.index('vorticity')
        self.ipsi = self.var.varname_list.index('psi')

        # background pv
        self.pvback = self.beta*(grid.yr-grid.Ly*.5)*grid.msk

        # for operators
        param.tracer_list=['pv']
        param.whosetspsi=('pvanom')
        param.qgoperator=True
        self.ope = Operators(param,grid)

        # for timescheme
        self.tscheme = Timescheme(param,self.var.state)
        self.dx0 = self.tscheme.dx0

        self.kt = 0

        if self.forcing:

            try:
                f = import_module(self.forcing_module)
            except:
                print('module %s for forcing cannot be found'%self.forcing_module)
                print('make sure file **%s.py** exists'%self.forcing_module)
                exit(0)

            self.forc = f.Forcing(param,grid)

        self.diags={}

        self.tscheme.set(self.dynamics, self.timestepping)


    def step(self,t,dt):
        self.dt = dt
        # 1/ integrate advection
        self.tscheme.forward(self.var.state,t,dt)

        # 2/ integrate source
        if self.noslip:
            self.add_noslip(self.var.state)

        # 3/ diagnostic fields
        self.var.state[self.ipva] = self.var.state[self.ipv] - self.pvback
        self.var.state[self.ivor] = self.var.state[self.ipva] + self.Rd**2*self.var.state[self.ipsi]


    def dynamics(self,x,t,dxdt):
        dxdt[:]=0.
        self.ope.rhs_adv(x,t,dxdt)
        if (self.tscheme.kstage==self.tscheme.kforcing):
            if self.forcing:
                self.forc.add_forcing(x,t,dxdt)
            if self.diffusion:
                self.ope.rhs_diffusion(x,t,dxdt)
        # substract the background pv
        dxdt[self.ipva][:,:] = dxdt[self.ipv] #- self.pvback*self.dt
#        self.ope.invert_vorticity(dxdt,flag='fast')
        self.ope.invert_vorticity(dxdt,flag='fast')
#        self.kt +=1
            
    def add_noslip(self,x):
        self.ope.rhs_noslip(x,self.source)        
        self.ope.invert_vorticity(x,flag='fast')
            

    def add_backgroundpv(self):
        self.var.state[self.ipv][:,:] = self.var.state[self.ipv] + self.pvback

    def set_psi_from_pv(self):
        state = self.var.state
        state[self.ipva] = state[self.ipv] - self.pvback
        self.ope.invert_vorticity(self.var.state,flag='full')


    def diagnostics(self,var,t):
        """ should provide at least 'maxspeed' (for cfl determination) """

        nh = self.nh
        u = var.get('u')
        v = var.get('v')
        trac = var.get('pv')

        ke,maxu = computekemaxu(self.msk,u,v,self.nh)
            
        z,z2    = computesumandnorm(self.msk,trac,self.nh)

        cst = self.mpitools.local_to_global( [(maxu,'max'),(ke,'sum'),(z,'sum'),(z2,'sum')] )

        self.diags['maxspeed'] = cst[0]            
        self.diags['ke']       = cst[1] / self.area
        self.diags['pv']       = cst[2] / self.area
        self.diags['pv2']      = cst[3] / self.area

