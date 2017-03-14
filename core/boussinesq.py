from param import Param
from operators import Operators
from variables import Var
from timescheme import Timescheme
from importlib import import_module
from fortran_diag import *
from numpy import sqrt,sum,zeros,tanh

class Boussinesq(Param):
    """ Boussinesq model 

    It provides the step(t,dt) function
    and 'var' containing the mode state
    """
    def __init__(self,param,grid):

        self.list_param=['forcing','noslip','timestepping','diffusion','Kdiff',
                         'forcing_module','gravity','isisland',
                         'customized','custom_module','additional_tracer']
        param.copy(self,self.list_param)

        # for potential energy
        self.list_param=['xr','yr','nh','Lx','msk','area','mpitools']
        grid.copy(self,self.list_param)

        # for variables
        param.varname_list=['vorticity','psi','u','v','buoyancy']
        param.tracer_list=['vorticity','buoyancy']
        param.whosetspsi=('vorticity')

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
        self.tscheme.set(self.dynamics, self.timestepping)

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
                print(f)
                self.extrastep = f.Step(param,grid)
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

        if self.customized:
            self.extrastep.do(self.var,t,dt)

    def dynamics(self,x,t,dxdt):
        self.ope.rhs_adv(x,t,dxdt)
        self.ope.rhs_torque(x,t,dxdt) # db/dx is a source term for the vorticity
        if (self.tscheme.kstage==self.tscheme.kforcing):
            if self.forcing:
                self.forc.add_forcing(x,t,dxdt)
            if self.diffusion:
                self.ope.rhs_diffusion(x,t,dxdt)

        self.ope.invert_vorticity(dxdt,flag='fast')


                
    def add_noslip(self,x):
        self.ope.rhs_noslip(x,self.source)        
        self.ope.invert_vorticity(x,flag='fast',island=self.isisland)
            
    def set_psi_from_vorticity(self):
        self.ope.invert_vorticity(self.var.state,island=self.isisland)

    def diagnostics(self,var,t):
        """ should provide at least 'maxspeed' (for cfl determination) """

        nh = self.nh
        u = var.get('u')
        v = var.get('v')
        vort = var.get('vorticity')
        buoy = var.get('buoyancy')

        ke,maxu = computekemaxu(self.msk,u,v,self.nh)
            
        z,z2    = computesumandnorm(self.msk,vort,self.nh)

        b,b2    = computesumandnorm(self.msk,buoy,self.nh)


        # potential energy (minus sign because buoyancy is minus density)
        pe = computesum(self.msk,buoy*self.yr,nh)
        pe = - self.gravity*pe 

        cst = self.mpitools.local_to_global( [(maxu,'max'),(ke,'sum'),(z,'sum'),(z2,'sum'),
                                          (pe,'sum'),(b,'sum'),(b2,'sum')] )

        
        self.diags['maxspeed'] = cst[0]
        self.diags['ke']       = cst[1] / self.area
        self.diags['pe']       = cst[4] / self.area
        self.diags['energy']   = (cst[1]+cst[4]) / self.area
        self.diags['vorticity']= cst[2] / self.area
        self.diags['enstrophy']= cst[3] / self.area
        self.diags['buoyancy'] = cst[5] / self.area
        self.diags['brms']= sqrt( cst[6] / self.area -(cst[5] / self.area)**2 )
