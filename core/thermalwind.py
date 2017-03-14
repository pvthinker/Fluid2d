from param import Param
from operators import Operators
from variables import Var
from timescheme import Timescheme

class Thermalwind(Param):
    """ Thermal wind model 

    It provides the step(t,dt) function
    and 'var' containing the mode state
    """
    def __init__(self,param,grid):

        self.list_param=['forcing','noslip','timestepping']
        param.copy(self,self.list_param)

        # for variables
        param.varname_list=('vorticity','psi','u','v','buoyancy','V')
        param.sizevar     =[grid.nyl,grid.nxl]
        self.var = Var(param)

        # for operators
        param.tracer_list=['vorticity','buoyancy','V']
        param.whosetspsi=('vorticity')
        self.ope = Operators(param,grid)

        # for timescheme
        self.tscheme = Timescheme(param,self.var.state)
        

    def step(self,t,dt):

        # 1/ integrate advection
        self.tscheme.set(self.dynamics, self.timestepping)
        self.tscheme.forward(self.var.state,t,dt)

        # 2/ integrate source
        if self.forcing or self.noslip:
            self.tscheme.set(self.sources, 'EF')
            self.tscheme.forward(self.var.state,t,dt)

        
    def dynamics(self,x,t,dxdt):
        self.ope.rhs_adv(x,t,dxdt)
        self.ope.rhs_thermalwind(x,t,dxdt) # add the r.h.s. terms
        self.ope.invert_vorticity(dxdt,flag='fast')
    

    def sources(self,x,t,dxdt):
        if self.forcing:
            self.ope.rhs_forcing(x,t,dxdt)
            self.ope.invert_vorticity(dxdt,flag='fast')
        if self.noslip:
            self.ope.rhs_noslip(x,t,dxdt)
            self.ope.invert_vorticity(x,flag='full')
            
    def set_psi_from_vorticity(self):
        self.ope.invert_vorticity(self.var.state)


