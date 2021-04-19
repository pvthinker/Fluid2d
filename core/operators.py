from param import Param
from numpy import zeros, zeros_like, roll
from gmg.hierarchy import Gmg
import fortran_advection as fa
import fortran_operators as fo
import fortran_diag as fd
from fourier import Fourier


class Operators(Param):

    def __init__(self, param, grid):
        self.list_param = ['varname_list', 'tracer_list',
                           'whosetspsi', 'mpi', 'npx', 'npy',
                           'nh', 'gravity', 'f0', 'beta', 'Rd',
                           'qgoperator', 'order', 'Kdiff', 'diffusion',
                           'enforce_momentum', 'isisland', 'aparab',
                           'flux_splitting_method', 'hydroepsilon',
                           'myrank', 'geometry', 'sqgoperator']

        param.copy(self, self.list_param)

        self.list_grid = ['msk', 'nxl', 'nyl', 'dx', 'dy',
                          'bcarea', 'mpitools', 'msknoslip',
                          'mskbc', 'domain_integration',
                          'nh', 'xr0', 'yr0', 'i0', 'j0', 'area']

        grid.copy(self, self.list_grid)
        self.first_time = True

        # internal work array for the inversion
        self.work = zeros((self.nyl, self.nxl))
        self.work2 = zeros((self.nyl, self.nxl))

        pp = {'np': param.npx, 'mp': param.npy, 'nh': param.nh,
              'n': param.nx//param.npx, 'm': param.ny//param.npy,
              'omega': 8./9., 'npmpmax': 1, 'verbose': False,
              'dx': grid.dx, 'dy': grid.dy, 'n1': 32, 'n0': 4,
              'method': 'deep', 'nagglo': 2,
              'hydroepsilon': param.hydroepsilon,
              'relaxation': param.relaxation}

        # load the multigrid solver
        #
        # WARNING: the multigrid needs the mask at cell corners!!!
        #         not at cell centers
        mskr = self.msk*1.

        # this piece is a bit awkward: to initialize gmg, we need
        # a mask with a halo properly filled but the fill_halo method
        # belongs to gmg. We have a circular definition.
        # the trick: define a dummy gmg first a msk=1 everywhere
        # then grab the fill_halo method and redefine once again the
        # multigrid, this time with the proper mask
        # self.gmg = Gmg(pp,mskr)
        # borrow the fill_halo from the multigrid
        # self.fill_halo = self.gmg.grid[0].halo.fill

        fo.celltocorner(mskr, self.work)
        # self.fill_halo(self.work)

        # del self.gmg
        # del self.fill_halo

        self.work[self.work < 1.] = 0.
        self.mskp = self.msk*0
        self.mskp[self.work == 1.] = 1
        pp['verbose'] = True
        if self.myrank == 0:
            print('-'*50)
            print(' Multigrid hierarchy')
            print('-'*50)

        if hasattr(self, 'qgoperator'):
            pp['qgoperator'] = True
            pp['Rd'] = self.Rd
            self.gmg = Gmg(pp, self.work)
        else:
            self.gmg = Gmg(pp, self.work)
        if hasattr(self, 'sqgoperator'):
            self.fourier = Fourier(param, grid)

        # borrow the fill_halo from the multigrid
        self.fill_halo = self.gmg.grid[0].halo.fill
        grid.fill_halo = self.gmg.grid[0].halo.fill

        self.blwidth = param.Lx*0.05

        # tentative for a regularized no-slip source term
        coef = 0.*zeros_like(self.mskp)
        coef[1:, 1:] = (self.mskp[:-1, 1:]+self.mskp[:-1, :-1]
                        + self.mskp[1:, 1:]+self.mskp[1:, :-1])
        # nbpsibc is the number of land psi-points surrounding a fluid cell
        self.nbpsibc = (4.-coef)*self.msk
        self.nbpsibc[self.nbpsibc > 0] = 1.

        self.set_boundary_msk()

        self.cst = zeros(5,)
        # select the proper flux discretization
        if self.order % 2 == 0:
            self.fortran_adv = fa.adv_centered
            self.cst[0] = grid.dx
            self.cst[1] = grid.dy
            self.cst[2] = 0.05
            self.cst[3] = 0  # umax
            self.cst[4] = 0  # unused
            # should be updated at each timestep
            # self.cst[3]=param.umax

        else:
            self.fortran_adv = fa.adv_upwind
            self.cst[0] = grid.dx
            self.cst[1] = grid.dy
            self.cst[2] = 0.05
            self.cst[3] = 0  # umax
            self.cst[4] = self.aparab
            # should be updated at each timestep
            # self.cst[3]=param.umax

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

        # these coefficients below are used for the thermalwind model
        coef = 0.*zeros_like(self.msk)
        coef[1:-1, 1:-1] = self.msk[1:-1, 2:]+self.msk[1:-1, 0:-2]
        coef[coef < 2] = 0.
        coef[coef == 2] = 0.5
        self.fill_halo(coef)
        self.coefb = coef*1.

        coef = 0.*zeros_like(self.msk)
        coef[1:-1, 1:-1] = self.msk[2:, 1:-1]+self.msk[0:-2, 1:-1]
        coef[coef < 2] = 0.
        coef[coef == 2] = 0.5
        self.fill_halo(coef)
        self.coefV = coef*1.

        if type(self.Kdiff) != dict:
            K = self.Kdiff
            self.Kdiff = {}
            for trac in self.tracer_list:
                self.Kdiff[trac] = K
        if self.diffusion:
            print('diffusion coefficients')
            print('  => ', self.Kdiff)

    def set_boundary_msk(self):
        """ for the no slip boundary source term """
        # nh = self.nh
        msk = self.msknoslip
        z = (roll(msk, -1, axis=1)+roll(msk, -1, axis=0)
             + roll(msk, +1, axis=1)+roll(msk, +1, axis=0)-4*msk)
        z = z*msk
        self.mskbc = self.msk*0
        self.mskbc[z < 0] = 1
        # the halo will be filled later in operator.py
        # when fill_halo will become available
        # we can now fix the boundary mask

        # to go with the new definition for the source term
        # self.mskbc = self.nbpsibc.copy()

        self.mskbc *= self.msknoslip
        self.fill_halo(self.mskbc)

        # idx in the 2D array where boundary terms are computed
        # used for storage and i/o
        # self.idxbc = where(self.mskbc==1)

        self.bcarea = self.domain_integration(self.mskbc)
        self.x2bc = self.domain_integration((self.xr0)**2
                                            * self.mskbc*self.msknoslip)
        self.y2bc = self.domain_integration((self.yr0)**2
                                            * self.mskbc*self.msknoslip)

        return

        def smooth(msk, msk0, k):
            y = (+roll(msk, -1, axis=1)+roll(msk, -1, axis=0)
                 + roll(msk, +1, axis=1)+roll(msk, +1, axis=0))
            z = msk*1.
            z[y > 0] = k+1
            z[msk > 0] = msk[msk > 0]
            z[msk0 == 0] = 0
            self.fill_halo(z)
            return z

        z0 = 1-msk*1.
        nk = int(round((self.blwidth/self.dx)))
        nk = 1
        for k in range(nk):
            z = z0*1.
            z0 = smooth(z, msk, k)
        z0[z0 == 0] = nk
        z0 = z0/nk
        z0 *= msk
        z0 = (1.-z0)**(nk/2.)
        z0[msk == 0] = 1
        self.cv = (roll(z0, -1, axis=1)+z0)*.5
        self.cu = (roll(z0, -1, axis=0)+z0)*.5
#        self.mskbc = z0
#        self.bcarea = self.domain_integration(z0)

    def rhs_adv(self, x, t, dxdt):
        """ compute -div(u*tracer) using finite volume flux discretization
        the flux is computed at edge cells using p-th order interpolation
        for p even, the flux is centered
        for p odd, the flux is upwinded (more points on the upwind side) """
        iu = self.varname_list.index('u')
        iv = self.varname_list.index('v')
        u = x[iu]
        v = x[iv]

        for trac in self.tracer_list:
            ik = self.varname_list.index(trac)
            y = dxdt[ik]
            self.fortran_adv(self.msk, x[ik], y, u, v,
                             self.cst, self.nh,
                             self.fs_method,
                             self.order)
            self.fill_halo(y)
            # for an unknown reason dxdt[ik] is
            # not updated by the Fortran routine
            # it should be done manually
            # (this yields an excessive data movement)
            dxdt[ik][:, :] = y

    def wallshear(self, x, shear):
        # ip = self.varname_list.index('psi')

        # meansource = fo.computewallshear(self.msk, x[ip],
        #                                  shear, self.dx, self.nh)
        return

    def rhs_noslip(self, x, source):
        """ add the vorticity source term along the boundary to enforce
        zero tangential velocity (=no-slip) """

        ip = self.varname_list.index('psi')
        # iu = self.varname_list.index('u')
        # iv = self.varname_list.index('v')
        iw = self.varname_list.index(self.whosetspsi)

        fo.cornertocell(x[ip], self.work)

        meansource = fo.computenoslipsourceterm(
            self.msknoslip, x[ip], self.work, self.dx, self.dy, self.nh)

        # K = self.dx*self.dy * 0.25
        # self.work2[:, :] = self.work[:, :]
        # for kt in range(2):
        #     fo.add_diffusion(self.msk, self.work, self.dx,
        #                      self.nh, K, self.work2)
        #     self.fill_halo(self.work2)
        #     fo.add_diffusion(self.msk, self.work2, self.dx,
        #                      self.nh, K, self.work)
        #     self.fill_halo(self.work)

        # # self.work =  self.work/(self.dx**2)*self.mskbc

        source[:, :] = self.work

        # this step is SUPER important to ensure GLOBAL vorticity conservation
        meansource = self.domain_integration(source) / self.bcarea

        source -= meansource*self.mskbc

        if self.enforce_momentum:
            xr = self.xr0
            yr = self.yr0
            # this step ensures the zero momentum
            px = fd.computedotprod(self.msk, source, xr, self.nh)
            py = fd.computedotprod(self.msk, source, yr, self.nh)
            cst = self.mpitools.local_to_global([(px, 'sum'), (py, 'sum')])

            px, py = cst[0]/self.x2bc, cst[1]/self.y2bc
            source -= (px*xr+py*yr)*self.mskbc

        self.fill_halo(source)
        x[iw] -= source

    def rhs_diffusion(self, x, t, dxdt, coef=1.):
        """ add a diffusion term on the tracer variables """

        for trac in self.tracer_list:
            ik = self.varname_list.index(trac)
            y = dxdt[ik]
            fo.add_diffusion(self.msk, x[ik], self.dx, self.nh,
                             coef*self.Kdiff[trac], y)
            self.fill_halo(y)
            dxdt[ik] = y

    def rhs_torque(self, x, t, dxdt):
        """ compute g*db/dx for the Boussinesq model """
        ib = self.varname_list.index('buoyancy')
        iw = self.varname_list.index('vorticity')

        y = dxdt[iw]
        b = x[ib]
        #y[1:-1, 1:-1] += self.gravity*self.diffx(b)
        y *= self.msk
        fo.add_torque(self.msk, b, self.dx, self.nh, self.gravity, y)
        self.fill_halo(y)
        dxdt[iw][:, :] = y

    def rhs_torque_density(self, x, t, dxdt):
        """ compute g*db/dx for the Boussinesq model """
        ib = self.varname_list.index('density')
        iw = self.varname_list.index('vorticity')

        y = dxdt[iw]
        b = x[ib]
        #y[1:-1, 1:-1] += self.gravity*self.diffx(b)
        # y *= self.msk
        # trick: use -gravity to account that density is opposite to buoyancy
        fo.add_torque(self.msk, b, self.dx, self.nh, -self.gravity, y)
        self.fill_halo(y)
        dxdt[iw][:, :] = y

    def diffx(self, x):
        nh = self.nh
        if self.i0 == self.npx-1:
            x[:, -nh] = 2*x[:, -nh-1]-x[:, -nh-2]
        if self.i0 == 0:
            x[:, nh-1] = 2*x[:, nh]-x[:, nh+1]
        return 0.5*(x[1:-1, 2:]-x[1:-1, :-2])/self.dx

    def diff1x(self, x):
        nh = self.nh
        if self.i0 == self.npx-1:
            x[:, -nh] = 2*x[:, -nh-1]-x[:, -nh-2]
        if self.i0 == 0:
            x[:, nh-1] = 2*x[:, nh]-x[:, nh+1]
        return (x[:, 1:]-x[:, :-1])/self.dx

    def diffz(self, x):
        nh = self.nh
        if self.j0 == self.npy-1:
            x[-nh, :] = 2*x[-nh-1, :]-x[-nh-2, :]
        if self.j0 == 0:
            x[nh-1, :] = 2*x[nh, :]-x[nh+1, :]
        return 0.5*(x[2:, 1:-1]-x[:-2, 1:-1])/self.dy

    def jacobian(self, x, y):
        return self.diffx(x)*self.diffz(y)-self.diffz(x)*self.diffx(y)

    def rhs_thermalwind(self, x, t, dxdt):

        iu = self.varname_list.index('u')
        ib = self.varname_list.index('buoyancy')
        iw = self.varname_list.index('vorticity')
        iV = self.varname_list.index('V')

        nh = self.nh

        # add the themal wind balance
        # g*db/dx + f0*dV/dz
        # to domega/dt
        b = x[ib]
        V = x[iV]
        dw = dxdt[iw]
        y = self.work
        # dw[1:-1, 1:-1] += self.diffx(b)*self.gravity
        # dw[1:-1, 1:-1] -= self.diffz(V)*self.f0

        y[1:-1, 1:-1] = self.diffx(b)*self.gravity
        y[1:-1, 1:-1] -= self.diffz(V)*self.f0

        # dw[1:-1, 1:-1] += self.coefb[1:-1, 1:-1]*self.diffx(b)*self.gravity
        # dw[1:-1, 1:-1] -= self.coefV[1:-1, 1:-1]*self.diffz(V)*self.f0
        # dw[:, :nh+1] = 0
        # dw[:, -nh-1:] = 0
        y *= self.msk
        self.fill_halo(y)
        dw[:, :] += y

        u = x[iu]
        # dxdt[iV][:, 1:] -= 0.5*self.f0*(u[:, :-1]+u[:, 1:])
        # dxdt[iV] *= self.msk
        # self.fill_halo(dxdt[iV])
        y[:, 1:] = - 0.5*self.f0*(u[:, :-1]+u[:, 1:])
        y *= self.msk
        self.fill_halo(y)
        dxdt[iV][:, :] += y

    def fourier_invert_vorticity(self, x, flag='full'):
        """ invert using Fourier transform """
        iu = self.varname_list.index('u')
        iv = self.varname_list.index('v')
        ip = self.varname_list.index('psi')
        ivor = self.varname_list.index('vorticity')
        ipv = self.varname_list.index('pv')

        u = x[iu]
        v = x[iv]
        psi = x[ip]
        pv = x[ipv]
        vor = x[ivor]

        self.fourier.invert(pv, psi, vor)
        self.fill_halo(psi)
        self.fill_halo(vor)

        self.first_time = False

        # compute (u,v) @ U,V points from psi @ cell corner
        fo.computeorthogradient(self.msk, psi, self.dx, self.dy, self.nh, u, v)
        x[iu] = u
        x[iv] = v

    def invert_vorticity(self, x, flag='full', island=False):
        """ compute psi from vorticity (or 'whosetspsi' in general)

        this routine interpolates the vorticity from cell centers to
        cell corners (where psi is defined)

        it then solves div*grad psi = omega with psi=0 along the boundary
        (Dirichlet condition) using a multigrid

        the non-divergent velocity is computed from psi"""
        iu = self.varname_list.index('u')
        iv = self.varname_list.index('v')
        ip = self.varname_list.index('psi')
        iw = self.varname_list.index(self.whosetspsi)

        u = x[iu]
        v = x[iv]
        psi = x[ip]

        fo.celltocorner(x[iw], self.work)
        #fo.celltocornerbicubic(x[iw], self.work)
        if island:
            # correcting RHS for islands
            self.work[:, :] -= self.rhsp

        if flag == 'fast':
            ite, res = self.gmg.twoVcycle(psi,
                                          self.work,
                                          {'maxite': 1,
                                           'tol': 1e-6,
                                           'verbose': True})
            # ite, res = self.gmg.solve(psi, self.work,
            #                           {'maxite': 2,
            #                            'tol': 1e-8,
            #                            'verbose': False})

        else:
            # compute to machine accuracy
            if self.first_time:
                verbose = True
            else:
                verbose = False
            if (self.myrank == 0) and verbose:
                print('-'*50)
                print(' Convergence of the vorticity inversion')
                print('    the residual should decrease by several orders')
                print('    of magnitude otherwise something is wrong')
                print('-'*50)

            ite, res = self.gmg.solve(psi, self.work,
                                      {'maxite': 4,
                                       'tol': 1e-11,
                                       'verbose': verbose})
            if self.geometry == 'perio':
                # make sure psi has zero mean (to avoid the drift)
                psim = self.domain_integration(psi) / self.area
                psi -= psim

        # don't apply the fill_halo on it
        # [because fill_halo, as it is, is applying periodic BC]
        psi = psi*self.mskp
        if island:
            # we set psi on the boundary values by adding
            # self.psi (defined in island module)
            # before that line, psi=0 along all boundaries
            psi += self.psi
            # it should be added only if we invert for the total psi
            # it should not be added if we compute the increment of psi

        self.first_time = False

        # compute (u,v) @ U,V points from psi @ cell corner
        fo.computeorthogradient(self.msk, psi, self.dx, self.dy, self.nh, u, v)
        # self.fill_halo(u)
        # self.fill_halo(v)
        x[iu] = u
        x[iv] = v
        x[ip] = psi
