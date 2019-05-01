import fortran_fluxes as fa
import fortran_operators as fo
import numpy as np
from timescheme import Timescheme
import copy
import sys

class Fluxes(object):

    def __init__(self, param, grid, ope):
        self.list_param = ['timestepping', 'varname_list',
                           'tracer_list', 'order', 'aparab',
                           'sizevar', 'flux_splitting_method']

        p = copy.deepcopy(param)
        newvariables = ['uc', 'vc', 'psic']
        #p.tracer_list += newvariables
        #p.varname_list += newvariables
        flx_list = []
        for v in p.tracer_list:
            for d in ['x', 'y']:
                flx_list += ['flx_%s_%s' % (v, d)]

        self.flx_list = flx_list
        self.nvarstate = len(p.varname_list)
        p.varname_list  += flx_list
        
        fullflx_list = []
        for r in ['rev', 'irr']:
            for v in p.tracer_list:
                for d in ['x', 'y']:
                    fullflx_list += ['%s_%s_%s' % (r, d, v)]

        self.fullflx_list = fullflx_list
        p.copy(self, self.list_param)

        self.list_grid = ['nh', 'dx', 'dy', 'msk']
        grid.copy(self, self.list_grid)

        self.ope = ope

        varsize = [len(self.varname_list)]+self.sizevar
        self.x = np.zeros(varsize)
        self.xwork = np.zeros(varsize)
        self.xwork2 = np.zeros(varsize)
        
        # tracsize = [len(flx_list)]+self.sizevar
        # # work array to store the fluxes after the forward phase and the backward phase
        # self.workflx = np.zeros(tracsize)
        # self.forward = np.zeros(tracsize)

        tracsize = [len(self.fullflx_list)]+self.sizevar
        # full stack of fluxes (including reversible and irreversible)
        self.flx = np.zeros(tracsize)
                
        # select the proper flux discretization
        if self.order % 2 == 0:
            self.fortran_adv = fa.adv_centered
        else:
            self.fortran_adv = fa.adv_upwind

        self.cst = np.zeros(5,)
        self.cst[0] = grid.dx
        self.cst[1] = grid.dy
        self.cst[2] = 0.05
        self.cst[3] = 0  # umax
        self.cst[4] = self.aparab
        
        # for timescheme
        self.tscheme = Timescheme(p, self.x)
        self.tscheme.set(self.advection, self.timestepping)

        # controls the flux splitting method
        # 0 = min/max
        # 1 = parabolic
        list_fs_method = ['minmax', 'parabolic']
        if self.flux_splitting_method in list_fs_method:
            self.fs_method = list_fs_method.index(
                self.flux_splitting_method)
        else:
            print('Warning: %s does not exist' % self.flux_splitting_method)
            print('replaced with the default: parabolic')
            self.fs_method = list_fs_method.index('parabolic')
        
        
        
    def diag_fluxes(self, x, t, dt):
        """advance x from one time step, then reverse the time and advance
        another time step. Diagnose the reversible and the
        irreversible part of the fluxes

        """
        iu = self.varname_list.index('u')
        iv = self.varname_list.index('v')
        ip = self.varname_list.index('psi')

        self.x[:self.nvarstate, :, :] = x
        # set fluxes entries to zero
        self.x[self.nvarstate:, :, :] = 0.       
        self.tscheme.forward(self.x, t, dt*1e-6)
        self.xwork[:, :, :] = self.x

        # reverse velocity
        # self.x[iu] *= -1
        # self.x[iv] *= -1
        # self.x[ip] *= -1
        self.x *= -1
        # set fluxes entries to zero
        self.x[self.nvarstate:, :, :] = 0.    
        self.tscheme.forward(self.x, t+dt, -dt*1e-6)
        self.xwork2[:, :, :] = self.x
        
        # # reverse velocity
        # self.x[iu] *= -1
        # self.x[iv] *= -1
        # self.x[ip] *= -1
        # # set fluxes entries to zero
        # self.x[self.nvarstate:, :, :] = 0.        
        # self.tscheme.forward(self.x, t-dt, dt)
        # self.x[self.nvarstate:, :, :] += self.xwork[self.nvarstate:, :, :]

        # don't forget to divide by 'dt' to get the tendency
        cff = 0.5/(dt*1e-6)
        nflx = len(self.flx_list)
        for k in range(nflx):
            l = self.nvarstate+k
            j = nflx+k
            self.flx[k, :, :] = cff*(self.xwork[l]-self.x[l])
            self.flx[j, :, :] = cff*(self.xwork[l]+self.x[l])
            # self.flx[k, :, :] = cff*(self.x[l]+self.xwork2[l])
            # self.flx[j, :, :] = cff*(self.x[l]-self.xwork2[l])
            
        
    def advection(self, x, t, dxdt):
        self.rhs_adv(x, t, dxdt)
        self.ope.invert_vorticity(dxdt, flag='fast')

    def rhs_adv(self, x, t, dxdt):
        """ compute -div(u*tracer) using finite volume flux discretization
        the flux is computed at edge cells using p-th order interpolation
        for p even, the flux is centered
        for p odd, the flux is upwinded (more points on the upwind side) """
        iu = self.varname_list.index('u')
        iv = self.varname_list.index('v')
        u = x[iu]
        v = x[iv]

        itrac = 0
        for trac in self.tracer_list:
            ik = self.varname_list.index(trac)
            y = dxdt[ik]
            ix = self.nvarstate+itrac*2
            iy = self.nvarstate+itrac*2+1            
            xf = dxdt[ix]
            yf = dxdt[iy]
            self.fortran_adv(self.msk, x[ik], y, u, v,
                             xf, yf,
                             self.cst, self.nh,
                             self.fs_method,
                             self.order)
            self.ope.fill_halo(y)
            # for an unknown reason dxdt[ik] is
            # not updated by the Fortran routine
            # it should be done manually
            # (this yields an excessive data movement)
            dxdt[ik][:, :] = y
            dxdt[ix, :, :] = xf
            dxdt[iy, :, :] = yf
            itrac += 1
