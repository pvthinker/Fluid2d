from param import Param
from operators import Operators
from variables import Var
from timescheme import Timescheme
from fortran_diag import *
from numpy import sqrt

class Advection(Param):
    """ Advection model 

    It provides the step(t,dt) function
    and 'var' containing the mode state
    """
    def __init__(self,param,grid):

        self.list_param=['timestepping','diffusion','Kdiff']
        param.copy(self,self.list_param)

        self.list_grid=['msk','nh','area','mpitools']
        grid.copy(self,self.list_grid)

        # for variables
        param.varname_list=('tracer','psi','u','v','vorticity')
        param.sizevar     =[grid.nyl,grid.nxl]
        self.var = Var(param)

        # for operators
        param.tracer_list=['tracer']
        param.whosetspsi=('tracer')
        self.ope = Operators(param,grid)

        # for timescheme
        self.tscheme = Timescheme(param,self.var.state)
        self.tscheme.set(self.advection, self.timestepping)        

        self.diags={}


    def step(self,t,dt):
        self.tscheme.forward(self.var.state,t,dt)
        self.diagnostics(self.var,t)
        
    def advection(self,x,t,dxdt):
        self.ope.rhs_adv(x,t,dxdt)

        if (self.tscheme.kstage==self.tscheme.kforcing):
            if self.diffusion:
                self.ope.rhs_diffusion(x,t,dxdt)

    def set_psi_from_tracer(self):
        self.ope.invert_vorticity(self.var.state)

                
    def diagnostics(self,var,t):
        """ should provide at least 'maxspeed' (for cfl determination) """

        if t==0.:
            # since the flow is prescribed, compute it only at t=0
            u = var.get('u')
            v = var.get('v')
            ke,maxu = computekemaxu(self.msk,u,v,self.nh)
            
            maxu = self.mpitools.local_to_global( [(maxu,'sum')] )

            self.diags['maxspeed'] = maxu
            
        trac = var.get('tracer')

        z,z2    = computesumandnorm(self.msk,trac,self.nh)

        cst = self.mpitools.local_to_global( [(z,'sum'),(z2,'sum')] )

        z  = cst[0] / self.area
        z2 = cst[1] / self.area

        self.diags['mean']=z
        self.diags['rms'] = sqrt( z2 )
