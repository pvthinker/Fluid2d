import numpy as np
from param import Param
from variables import Var


class Timescheme(object):
    """ Catalog of time schemes

    It provides a forward() method that uses a generic rhs function """

    def __init__(self, param, x):
        self.list_param = ['timestepping']
        param.copy(self, self.list_param)

        self.timeschemelist = {'EF': self.EulerForward,
                               'LF': self.LeapFrog,
                               'Heun': self.Heun,
                               'AB2': self.AB2,
                               'AB3': self.AB3,
                               'LFAM3': self.LFAM3,
                               'RK3_SSP': self.RK3_SSP,
                               'RK3': self.RK3,
                               'RK4_LS': self.RK4_LS}

        self.asselin_cst = 0.1
        self.ab2_epsilon = 0.1

        # coef to apply on terms that are evaluated only at the
        # last stage of a multi-stages scheme
        coefondt = {'EF': 1., 'LF': 1., 'Heun': 2.,
                    'AB2': 1., 'AB3': 1., 'LFAM3': 1.,
                    'RK3_SSP': 1.5, 'RK3': 1., 'RK4_LS': 1.}
        self.dtcoef = coefondt[self.timestepping]

        # internal arrays
        self.x = np.zeros_like(x)
        self.dx0 = np.zeros_like(x)
        self.dx1 = np.zeros_like(x)

        if self.timestepping in ['RK3_SSP', 'AB3', 'RK3']:
            self.dx2 = np.zeros_like(x)

        if self.timestepping in ['LF', 'LFAM3']:
            self.xb = np.zeros_like(x)

        self.first = True
        self.second = True

        # the rhs and the timestepping has to be provided later
        # via the self.set()
        def f(**kwargs):
            print('define a rhs and a timestepping before calling forward()')
            print('I''m stopping now !!!')
            exit(0)
        # self.forward will be overwritten by self.set()
        self.forward = f

        self.kstage = 0

    # ----------------------------------------
    def set(self, rhs, timestepping):
        """ assign the 'rhs' & the 'timestepping' to the 'forward' method """
        self.rhs = rhs
        self.forward = self.timeschemelist[timestepping]

        self.kforcing = 0

        if self.timestepping in ['RK4_LS']:
            self.kforcing = 3

        elif self.timestepping in ['RK3_SSP']:
            self.kforcing = 2

        elif self.timestepping in ['Heun', 'LFAM3']:
            self.kforcing = 1

    # ----------------------------------------
    def EulerForward(self, x, t, dt, **kwargs):
        self.rhs(x, t, self.dx0)
        x += dt * self.dx0

    # ----------------------------------------
    def AB2(self, x, t, dt):
        self.rhs(x, t, self.dx0)
        if self.first:
            x += dt * self.dx0
            self.first = False
        else:
            x += ((1.5+self.ab2_epsilon)*dt) * self.dx0\
                - ((0.5+self.ab2_epsilon)*dt)*self.dx1
        self.dx1[:] = self.dx0

    # ----------------------------------------
    def AB3(self, x, t, dt, **kwargs):
        self.rhs(x, t, self.dx0)
        if self.first:
            x += dt * self.dx0
            self.first = False
        elif self.second:
            x += (1.5*dt) * self.dx0 - (0.5*dt)*self.dx1
            self.second = False
        else:
            x += (23*dt/12.) * self.dx0\
                - (16*dt/12.)*self.dx1+(5*dt/12.)*self.dx2
        self.dx2[:] = self.dx1
        self.dx1[:] = self.dx0

    # ----------------------------------------
    def LeapFrog(self, x, t, dt, **kwargs):
        self.x[:] = x
        self.rhs(x, t, self.dx0)
        if self.first:
            x += dt * self.dx0
            self.first = False
        else:
            x[:] = self.xb + (2*dt) * self.dx0
            # Asselin's filter
            self.x += self.asselin_cst*(x+self.xb-2*self.x)
        self.xb[:] = self.x

    # ----------------------------------------
    def LFAM3(self, x, t, dt, **kwargs):
        self.x[:] = x
        self.kstage = 0
        self.rhs(x, t, self.dx0)
        if self.first:
            x += dt * self.dx0
            self.first = False
        else:
            # LF x is at n+1
            x[:] = self.xb + (2*dt) * self.dx0
            # AM3 gives x at n+1/2
            x[:] = (1./12.)*(5.*x + 8.*self.x-self.xb)
            # do the corrector step at n+1/2
            self.kstage = 1
            self.rhs(x, t+dt*.5, self.dx0)
            x[:] = self.x + dt*self.dx0
        self.xb[:] = self.x

    # ----------------------------------------
    def Heun(self, x, t, dt, **kwargs):
        # predictor
        self.kstage = 0
        self.rhs(x, t, self.dx0)
        self.x = x + dt * self.dx0
        # corrector
        self.kstage = 1
        self.rhs(self.x, t+dt, self.dx1)
        x += (0.5*dt)*(self.dx0+self.dx1)

    # ----------------------------------------
    def RK3(self, x, t, dt, **kwargs):
        # "classical" RK3
        self.kstage = 0
        self.rhs(x, t, self.dx0)
        self.x = x + (dt/3.) * self.dx0

        self.kstage = 1
        self.rhs(self.x, t+dt/3., self.dx1)
        self.x = x + (0.5*dt)*self.dx1

        self.kstage = 2
        self.rhs(self.x, t+0.5*dt, self.dx2)
        x += dt*self.dx2

    # ----------------------------------------
    def RK3_SSP(self, x, t, dt, **kwargs):
        # Stably Strongly Preserving RK3
        # convex combination of rhs
        self.kstage = 0
        self.rhs(x, t, self.dx0)
        self.x = x + dt * self.dx0

        self.kstage = 1
        self.rhs(self.x, t+dt, self.dx1)
        self.x = x + (0.25*dt)*(self.dx0+self.dx1)

        self.kstage = 2
        self.rhs(self.x, t+0.5*dt, self.dx2)
        x += (dt/6.)*(self.dx0+self.dx1+4*self.dx2)

    # ----------------------------------------
    def RK4_LS(self, x, t, dt, **kwargs):
        # low-storage RK4
        # 4th order for linear waves
        # 2nd order for nonlinear dynamics
        self.kstage = 0
        self.rhs(x, t, self.dx0)
        self.x = x + (0.25*dt) * self.dx0

        self.kstage = 1
        self.rhs(self.x, t+dt*0.25, self.dx0)
        self.x = x + (dt/3.)*self.dx0

        self.kstage = 2
        self.rhs(self.x, t+dt/3., self.dx0)
        self.x = x + (dt/2.)*self.dx0

        self.kstage = 3
        self.rhs(self.x, t+0.5*dt, self.dx0)
        x += dt*self.dx0


if __name__ == "__main__":

    param = Param()
    param.timestepping = 'AB2'
    param.varname_list = ['vort', 'psi', 'u']
    param.sizevar = [4, 2]

    var = Var(param)

    tscheme = Timescheme(param, var.state)

    print(tscheme.x)
    tscheme.forward()
