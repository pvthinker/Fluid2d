from param import Param
from operators import Operators
from variables import Var
from timescheme import Timescheme
from numpy import zeros,sqrt,sum,shape,tanh
from importlib import import_module
from fortran_diag import *


class Euler(Param):
    """ Euler model 

    It provides the step(t,dt) function
    and 'var' containing the mode state
    """
    def __init__(self,param,grid):

        self.list_param=['forcing','noslip','timestepping','diffusion','Kdiff',
                         'forcing_module','additional_tracer','enforce_momentum','var_to_save',
                         'customized','custom_module']
        param.copy(self,self.list_param)

        # for diagnostics
        self.list_param=['yr','nh','msk','area','mpitools','dx','xr','yr','r2','x0','y0','x2','y2','isisland','Lx']
        grid.copy(self,self.list_param)

        # for variables
        param.varname_list=['vorticity','psi','u','v']#,'tauw','wshear']
        param.tracer_list=['vorticity']
        param.whosetspsi=('vorticity')

        if 'tauw' in self.var_to_save:
            param.varname_list.append('tauw')

        if 'wshear' in self.var_to_save:
            param.varname_list.append('wshear')

        if hasattr(self,'additional_tracer'):
            for k in range(len(self.additional_tracer)):
                trac=self.additional_tracer[k]
                param.varname_list.append(trac)
                param.tracer_list.append(trac)
                print('Tracers are :',param.tracer_list)

        param.sizevar     =[grid.nyl,grid.nxl]
        self.var = Var(param)
        self.source = zeros(param.sizevar)


        # for operators
        self.ope = Operators(param,grid)

        # for timescheme
        self.tscheme = Timescheme(param,self.var.state)
        self.tscheme.set(self.advection, self.timestepping)

        if self.forcing:

            try:
                f = import_module(self.forcing_module)
            except:
                print('module %s for forcing cannot be found'%self.forcing_module)
                print('make sure file **%s.py** exists'%self.forcing_module)
                exit(0)

            self.forc = f.Forcing(param,grid)

        self.diags={}

        if self.customized:
            try:
                f = import_module(self.custom_module)
                self.extrastep = f.step(param,grid)
            except:
                print('module %s for forcing cannot be found'%self.custom_module)
                print('make sure file **%s.py** exists'%self.custom_module)
                exit(0)


    def step(self,t,dt):

        # 1/ integrate advection
        self.tscheme.forward(self.var.state,t,dt)

        # 2/ integrate source
        if self.noslip:
            self.add_noslip(self.var.state)
        
            if 'tauw' in self.var_to_save:
                tauw = self.var.get('tauw')
                tauw[:,:] = self.source/dt*self.dx*self.dx

            if 'wshear' in self.var_to_save:
                wshear = self.var.get('wshear')
                self.ope.wallshear(self.var.state,self.source) 
                wshear[:,:]=self.source            

        if self.customized:
            self.extrastep.do(self.var,t,dt)

        # # 3/ dye tracer
        # i,jp,jm=1,67,61
        
        # dye = self.var.get('dye')
        # dye[jp+self.nh,i]=1
        # dye[jm+self.nh,i]=-1

        # age = self.var.get('age')
        # age += dt*self.msk
        # age[:,self.nh]=0.


        # # 4/ sponge layer
        # damping = ( 1- (1+tanh( (self.xr - self.Lx)/0.1))*0.5)
        # w = self.var.get('vorticity')
        # w *= damping
        # dye *= damping
        # age *= damping
        

        
    def advection(self,x,t,dxdt):
        self.ope.rhs_adv(x,t,dxdt)
        if (self.tscheme.kstage==self.tscheme.kforcing):
            if self.forcing:
                self.forc.add_forcing(x,t,dxdt)
            if self.diffusion:
                self.ope.rhs_diffusion(x,t,dxdt)

        self.ope.invert_vorticity(dxdt,flag='fast')
    

    def add_noslip(self,x):
        self.ope.rhs_noslip(x,self.source)        
        if not(self.enforce_momentum):
            self.ope.invert_vorticity(x,flag='fast',island=self.isisland)
            
    def set_psi_from_vorticity(self):
        self.ope.invert_vorticity(self.var.state,island=self.isisland)


    def diagnostics(self,var,t):
        """ should provide at least 'maxspeed' (for cfl determination) """

        nh = self.nh
        u = var.get('u')
        v = var.get('v')
        trac = var.get('vorticity')
        
        xr = self.xr
        yr = self.yr
        r2 = self.r2

        ke,maxu = computekemaxu(self.msk,u,v,self.nh)            

        z,z2    = computesumandnorm(self.msk,trac,self.nh)

        px = computedotprod(self.msk,trac,xr,self.nh)
        py = computedotprod(self.msk,trac,yr,self.nh)
        angmom = computedotprod(self.msk,trac,r2,self.nh)

        cst = self.mpitools.local_to_global( [(maxu,'max'),(ke,'sum'),(z,'sum'),(z2,'sum'),
                                          (px,'sum'),(py,'sum'),(angmom,'sum')] )

        self.diags['maxspeed'] = cst[0]
        self.diags['ke']       = cst[1] / self.area
        self.diags['vorticity']= cst[2] / self.area
        self.diags['enstrophy']= cst[3] / self.area
        self.diags['px']       = cst[4] / self.area
        self.diags['py']       = cst[5] / self.area
        self.diags['angmom']   = cst[6] / self.area


        
