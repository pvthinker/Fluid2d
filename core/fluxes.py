import copy
import numpy as np
import fortran_fluxes as fa
from timescheme import Timescheme

class Fluxes():
    """Class to diagnostic the irreversible part of advective fluxes

    The development status is
    - 'euler' : ok
    - 'boussinesq' : ok
    - 'thermalwind' : ok
    - 'qg' : not done

    Note that for 'thermalwind' we should change sign of both V and f0,
    but we don't, because the equations are still invariant under time
    reversal if we don't do it. The rationale is that the implementation
    is easier.
"""
    def __init__(self, param, grid, ope):
        self.list_param = ['timestepping', 'varname_list',
                           'tracer_list', 'order', 'aparab',
                           'sizevar', 'flux_splitting_method',
                           'modelname']

        p = copy.deepcopy(param)
        # we add two new variables in the state vector: uc and vc
        # the velocities at cell centers
        # their fluxes are used to estimate the KE dissipation
        newvariables = ['uc', 'vc']
        p.tracer_list += newvariables
        p.varname_list += newvariables
        flx_list = []
        for v in p.tracer_list:
            for d in ['x', 'y']:
                flx_list += ['flx_%s_%s' % (v, d)]

        self.flx_list = flx_list
        self.nvarstate = len(p.varname_list)
        p.varname_list += flx_list

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

        # extended sate vector to host the model state vector and (uc,vc)
        self.xe = np.zeros([self.nvarstate]+self.sizevar)
        self.xwork = np.zeros(varsize)

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
        #p.timestepping = 'EF'
        self.tscheme = Timescheme(p, self.x)
        self.tscheme.set(self.advection, p.timestepping)

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
        iw = self.varname_list.index('vorticity')

        # set locally dt to a very small value to make the time scheme
        # dissipation vanishing (hard coded so far)
        #
        # dt = dt*1e-6
        self.xe[:self.nvarstate-2, :, :] = x
        #
        # we add two more variables 'uc' and 'vc', the cell centered velocities
        # we store them in the extended state self.xe
        #
        # i) uc
        y = 0.5*(x[iu]+np.roll(x[iu], 1, axis=1))
        self.ope.fill_halo(y)
        self.xe[self.nvarstate-2, :, :] = y
        # ii) vc
        y = 0.5*(x[iv]+np.roll(x[iv], 1, axis=0))
        self.ope.fill_halo(y)
        self.xe[self.nvarstate-1, :, :] = y

        # we copy this extended state into self.x
        self.x[:self.nvarstate, :, :] = self.xe

        # set fluxes entries to zero
        self.x[self.nvarstate:, :, :] = 0.
        self.tscheme.forward(self.x, t, dt)
        self.xwork[:, :, :] = self.x

        # reverse velocity
        self.x[:self.nvarstate, :, :] = self.xe
        self.x[iu] *= -1
        self.x[iv] *= -1
        self.x[ip] *= -1
        self.x[iw] *= -1
        self.x[self.nvarstate-2] *= -1 # uc
        self.x[self.nvarstate-1] *= -1 # vc

        # set fluxes entries to zero
        self.x[self.nvarstate:, :, :] = 0.
        self.tscheme.forward(self.x, t+dt, -dt)

        # don't forget to divide by 'dt' to get the tendency
        cff = 0.5/(dt)
        nflx = len(self.flx_list)
        # this loop below should be made more readable
        # k sweeps across all reversible fluxes
        # (two components per advected quantity)
        for k in range(nflx):
            # l points to the flux index in the self.x array
            l = self.nvarstate+k
            # j points to the associated irreversible flux
            j = nflx+k
            if (k < 2) or (k >= (nflx-4)):
                # k = 0, 1 are the vorticity fluxes
                # k = nflx-4, nflx-3 are the 'uc' fluxes
                # k = nflx-2, nflx-1 are the 'vc' fluxes
                sign = -1 #for vorticity, 'uc' and 'vc'
            else:
                sign = 1 # for buoyancy and tracer (if present)

            self.flx[k, :, :] = cff*(self.xwork[l]+sign*self.x[l])
            self.flx[j, :, :] = cff*(self.xwork[l]-sign*self.x[l])


    def advection(self, x, t, dxdt):
        self.rhs_adv(x, t, dxdt)
        if self.modelname == 'boussinesq':
            self.ope.rhs_torque(x, t, dxdt)
        elif self.modelname == 'thermalwind':
            self.ope.rhs_thermalwind(x, t, dxdt)
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
        # don't forget to adjust the maxspeed
        self.cst[3] = self.ope.cst[3]
        for itrac, trac in enumerate(self.tracer_list):
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
            self.ope.fill_halo(xf)
            self.ope.fill_halo(yf)
            # for an unknown reason dxdt[ik] is
            # not updated by the Fortran routine
            # it should be done manually
            # (this yields an excessive data movement)
            dxdt[ik, :, :] = y
            dxdt[ix, :, :] = xf
            dxdt[iy, :, :] = yf
