###########################################################
#
# functions relative to multigrid
#
############################################################
try:
    from mpi4py import MPI
except:
    print('- mpi4py : not found, please install it')
    exit(0)

from numpy import zeros, arange, NaN
from gmg.level import *
#from analyticals import *
#from plotutils import *
from time import time

#from scipy.signal import convolve2d


class Gmg(object):

    def __init__(self, param, mskp):

        info = Gridinfo(param)
        self.info = info
        self.nlevs = len(info)

        self.npre = 1
        self.npost = 1
        self.nvcyc = 1
        self.ndeepest = 16
        self.nglo = param['n']*param['np']
        self.mglo = param['m']*param['mp']
        self.myrank = MPI.COMM_WORLD.Get_rank()
        self.setup_hierarchy(param, mskp)

    def setup_hierarchy(self, param, mskp):
        """ main initialisation part of the multigrid """
        # define the grids hierarchy
        self.grid = []
        for lev in range(self.nlevs):
            self.grid.append(Grid(self.info[lev]))

        g = self.grid[0]
        msk = mskp.astype(int8)
#        msk=zeros( (g.mv,g.nv),dtype=int8)+1
#        name = param['name']
#        msk = set_fine_msk(param,self.myrank,name)

#        self.smooth_msk(msk)

        self.x = []
        self.b = []
        self.r = []
        # populate them
        for lev in range(self.nlevs):
            if lev == 0:
                self.grid[lev].set_msk_laplacian(prevlev=0, mskf=msk)
            else:
                self.grid[lev].set_msk_laplacian(prevlev=self.grid[lev-1])
            mv = self.grid[lev].mv
            nv = self.grid[lev].nv
            # allocate memory for internal arrays
            # x is the field
            # b is the RHS
            # r is the residual r=b-A*x
            self.x.append(zeros((mv, nv)))
            self.b.append(zeros((mv, nv)))
            self.r.append(zeros((mv, nv)))

#            self.grid[lev].halo.init_fill_optimized(self.x[lev],'x')

        # here we truncate  A[:,:,0:8] into A[:,:0:4] on all grids
        # to take advantage of the fact that A is a symmetric matrix
        for lev in range(self.nlevs):
            self.grid[lev].ninetofive()

        # transform the Laplacian into a Helmholtz operator if requiered
        if 'qgoperator' in param.keys():
            for lev in range(self.nlevs):
                #print('****************Helmholtz operator')
                msk = self.grid[lev].A[:, :, 4]
                self.grid[lev].A[msk != 0., 4] -= 1./(param['Rd']**2)

# ----------------------------------------
    def smooth_msk(self, msk):

        coef = msk*1.
        weight = ones((3, 3))/9.
        for k in range(15):
            coef = convolve2d(coef, weight, mode='same', boundary='wrap')
            self.grid[0].halo.fill(coef)
        msk[coef <= 0.5] = 0
        msk[coef > 0.5] = 1

# ----------------------------------------
    def Vcycle(self, lev1):
        deepest_lev = self.nlevs-1
        for lev in range(lev1, deepest_lev):
            if lev > lev1:
                self.x[lev][:, :] = 0.
            if self.grid[lev].flag != 'peak':
                self.grid[lev].smooth(self.x[lev], self.b[lev], self.npre)
                self.grid[lev].residual(self.x[lev], self.b[lev], self.r[lev])
                # the residual will be the RHS on the next grid
                finetocoarse(self.grid[lev], self.grid[lev+1],
                             self.r[lev], self.b[lev+1])
            else:
                finetocoarse(self.grid[lev], self.grid[lev+1],
                             self.b[lev], self.b[lev+1])

        # at the coarsest level we simply smooth, exact solution is not requiered
        lev = deepest_lev
        self.x[lev][:, :] = 0.
        self.grid[lev].smooth(self.x[lev], self.b[lev], self.ndeepest)

        for lev in range(deepest_lev-1, lev1-1, -1):
            if self.grid[lev].flag == 'peak':
                coarsetofine(self.grid[lev+1], self.grid[lev],
                             self.x[lev+1], self.x[lev])
            else:
                coarsetofine(self.grid[lev+1], self.grid[lev],
                             self.x[lev+1], self.r[lev])
                self.x[lev] += self.r[lev]
                self.grid[lev].smooth(self.x[lev], self.b[lev], self.npost)


# ----------------------------------------

    def Fcycle(self, lev1):
        # push the RHS down on all grids
        deepest_lev = self.nlevs-1
        for lev in range(lev1, deepest_lev):
            finetocoarse(self.grid[lev], self.grid[lev+1],
                         self.b[lev], self.b[lev+1])

        # start solving the pb on the coarsest grid
        lev = deepest_lev
        self.x[lev][:, :] = 0.
        self.grid[lev].smooth(self.x[lev], self.b[lev], self.ndeepest)

        # pull the solution up
        for lev in range(deepest_lev-1, lev1-1, -1):
            coarsetofine(self.grid[lev+1], self.grid[lev],
                         self.x[lev+1], self.x[lev])
            if self.grid[lev].flag != 'peak':
                # improve the solution by applying one or two Vcycle (down from this level)
                for k in range(self.nvcyc):
                    # ncyc ==2 corresponds to a W-cycle
                    self.Vcycle(lev)

# ----------------------------------------
    def solve(self, x, b, param):
        """ driver to solve A*x=b using the Fcycle
        x is the first guess, b the right hand side"""

        g = self.grid[0]
        g.residual(x, b, self.b[0])
        normeb = g.norm(b)
        if normeb > 0:
            res0 = (g.norm(self.b[0])/normeb)
        else:
            return 0, 0.
        res = res0
        nite = 0
        nite_diverge = 0
        # improve the solution until one of this condition is met
        while (nite < param['maxite']) and (res0 > param['tol']):
            self.Fcycle(0)
            x += self.x[0]
            g.residual(x, b, self.b[0])

            res = (g.norm(self.b[0])/normeb)
            conv = res0/res

            res0 = res
            nite += 1
            if (self.myrank == 0) and (param['verbose']):
                print(' ite = {} / res = {:.2e} / conv = {:8.4f}'.format(nite, res, conv))

            if (conv < 1):
                nite_diverge += 1

            if (nite_diverge > 4):
                if (self.myrank == 0):
                    print('solver is not converging')
                exit(0)
            # zglo=combine_global(g,self.b[0])
            # if myrank==0:
            #     plot2d(zglo,'zglobal')
        return nite, res

# ----------------------------------------
    def oneFcycle(self, x, b, param):
        """ driver to solve A*x=b using one Fcycle
        x is the first guess, b the right hand side"""

        g = self.grid[0]
        g.residual(x, b, self.b[0])
        self.x[0][:, :] = 0.
        self.Fcycle(0)
        x += self.x[0]
        return 1, 0.

# ----------------------------------------
    def twoVcycle(self, x, b, param):
        """ driver to solve A*x=b using two Vcycles
        x is the first guess, b the right hand side"""

        g = self.grid[0]
        self.x[0][:, :] = x
        self.b[0][:, :] = b
        for k in range(2):
            g.residual(self.x[0], self.b[0], self.r[0])
            self.Vcycle(0)
        x[:, :] = self.x[0]
        return 1, 0.

# ----------------------------------------
    def Vsolve(self, x, b, param):
        """ driver to solve A*x=b using V-cycles
        x is the first guess, b the right hand side"""

        g = self.grid[0]
        self.x[0][:, :] = x
        self.b[0][:, :] = b
        g.residual(x, b, self.r[0])
        normeb = g.norm(b)
        res0 = (g.norm(self.r[0])/normeb)

        nite = 0
        # improve the solution until one of this condition is met
        while (nite < param['maxite']) and (res0 > param['tol']):
            self.Vcycle(0)
            g.residual(self.x[0], self.b[0], self.r[0])

            res = (g.norm(self.r[0])/normeb)
            conv = res0/res

            res0 = res
            nite += 1
            if (self.myrank == 0) and (param['verbose']):
                print(' ite = %i / res = %g / conv = %g' % (nite, res, conv))

            #zglo = combine_global(g, self.r[0])
#            if myrank==0:
#                plot2d(zglo,'zglobal')

        x[:, :] = self.x[0]
        return nite, res

# ----------------------------------------
    def coarsen(self, x, y):
        # x is a 2D field on lev=0
        # y is a grid structure
        y[0][:, :] = x[:, :]
        for lev in range(self.nlevs-1):
            finetocoarse(self.grid[lev], self.grid[lev+1], y[lev], y[lev+1])

# ----------------------------------------
    def innerproduct(self, x, y):
        # x and y are two grid structure
        # y is a grid structure
        z = zeros((self.nlevs,))
        for lev in range(self.nlevs):
            z[lev] = (4.**lev)*self.grid[lev].inner(x[lev], y[lev])
        return z


# ----------------------------------------


    def generate_random_large_scale_field(self):
        """ mg is a gmg object """
        from scipy import random

        # level on which the random field is generated
        # this controls the scale of the field
        # current limitation: this has to be done on a level for which
        # each node has the whole domain
        flev = self.nlevs-1  # min(6,self.nlevs-1)

        # gaussian noise parameters
        mu = 0.
        sigma = 1.

        # generate it on rank==0 then broadcast it
        # so that every rank has the same field
        if (self.myrank == 0):
            random.seed(1)
            forc = random.normal(
                mu, sigma, (self.grid[flev].mv, self.grid[flev].nv))
            forc = forc*self.grid[flev].msk
        else:
            forc = None
        forc = MPI.COMM_WORLD.bcast(forc, root=0)

        # interpolate it on the finest grid
        self.x[flev][:, :] = forc
        for lev in range(flev-1, -1, -1):
            if self.grid[lev].flag == 'peak':
                coarsetofine(self.grid[lev+1], self.grid[lev],
                             self.x[lev+1], self.x[lev])
            else:
                coarsetofine(self.grid[lev+1], self.grid[lev],
                             self.x[lev+1], self.r[lev])
                self.x[lev] += self.r[lev]
                # 15oct changed ,1 to ,2
                self.grid[lev].smooth(self.x[lev], self.b[lev], 2)

        return self.x[0]

# ----------------------------------------
    def generate_random_small_scale_field(self):
        """ mg is a gmg object """
        from scipy import random

        # level on which the random field is generated
        # this controls the scale of the field
        # current limitation: this has to be done on a level for which
        # each node has the whole domain
        lev = 0  # min(6,self.nlevs-1)

        # gaussian noise parameters
        mu = 0.
        sigma = 1.

        # generate it on rank==0 then broadcast it
        # so that every rank has the same field

        forc = random.normal(mu, sigma, (self.grid[lev].mv, self.grid[lev].nv))

        # smooth twice on the finest grid
        self.x[lev][:, :] = forc
        self.grid[lev].smooth(self.x[lev], self.b[lev],
                              2)  # 15oct changed ,1 to ,2

        return self.x[0]


# ----------------------------------------


    def print_timers(self):
        c1 = 0
        c2 = 0
        c3 = 0
        c4 = 0
        print('----------------------------------------------------------')
        for lev in range(self.nlevs):
            tt = self.grid[lev].time
            cc = self.grid[lev].ncalls
            comp = tt['smooth']+tt['res']+tt['sum'] + \
                tt['norm']+tt['interpolate']+tt['restrict']
            comp = tt['smooth']
            comm = tt['halo']+tt['reduce']+tt['split']+tt['gather']
            halo = tt['halo']
            barr = tt['barrier']
            c1 += comp
            c2 += comm
            c3 += halo
            c4 += barr
            print('%2i / comp =%5.3f / comm =%5.3f / halo =%5.3f / barrier =%5.3f / total=%5.3f' %
                  (lev, comp, comm, halo, barr, comp+comm))

        print('----------------------------------------------------------')
        print('total      %5.3f /       %5.3f /       %5.3f /       %5.3f /       %5.3f' %
              (c1, c2, c3, c4, c1+c2))
        print('----------------------------------------------------------')
        for lev in range(self.nlevs):
            cc = self.grid[lev].ncalls
            print('%2i / nsmooth = %4i / nres = %4i/ nhalos = %4i' %
                  (lev, cc['smooth'], cc['res'], cc['halo']))
        print('----------------------------------------------------------')

        for lev in range(self.nlevs):
            tt = self.grid[lev].time['halo']
            cc = self.grid[lev].ncalls['halo']
            if cc == 0:
                cc = 1
            print('%2i / halo = %5.3e' % (lev, tt/cc))
        print('--------------------------------------------')

    def benchmark_comm(self):
        nite = 1000
        comm1 = zeros((self.nlevs, 1))
        comm2 = zeros((self.nlevs, 1))
        for lev in range(self.nlevs):
            if self.grid[lev].flag != 'peak':
                x = self.x[lev]
                b = self.b[lev]
                t0 = time.time()
                for kt in range(nite):
                    self.grid[lev].halo.fill(x)
                    self.grid[lev].halo.fill(b)
                t1 = time.time()
                comm1[lev] = (t1-t0)/2

        for kt in range(nite):
            for lev in range(self.nlevs):
                if self.grid[lev].flag != 'peak':
                    x = self.x[lev]
                    b = self.b[lev]
                    t0 = time.time()
                    self.grid[lev].halo.fill(x)
                    self.grid[lev].halo.fill(b)
                    t1 = time.time()
                    comm2[lev] += (t1-t0)/2

        for lev in range(self.nlevs):
            if self.myrank == 0:
                print('%2i / halo = %5.3e / %5.3e' %
                      (lev, comm1[lev]/nite, comm2[lev]/nite))


# ----------------------------------------
if __name__ == '__main__':
    #    from plotutils import plot2d
    from numpy import zeros, arange, sqrt, exp, NaN, log
    import time

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nbofproc = MPI.COMM_WORLD.Get_size()

    # this is where we decide the subdomains decomposition
    n0 = int(sqrt(nbofproc))
    if n0**2 == nbofproc:
        # if nbofproc is a power of 4, then square decomposition
        np, mp = n0, n0
    else:
        # if nbofproc is of the form 2 *(4**p), then
        n0 = int(sqrt(nbofproc/2))
        if n0**2 == nbofproc/2:
            np, mp = n0, 2*n0  # twice more procs in the 'y' direction
        else:
            # not a power of 2 => error
            if myrank == 0:
                print('pb with number of procs, should be 2**p')
            exit()

    nh = 3  # width of the halo (2 or 3)

    n = 512*2/np
    m = 512*2/mp

    n = 512*4/np
    m = 512*4/mp

    dx = 1.

    name = 'headlands'
#    name='perio'

    Lx = n*np*dx
    Ly = m*mp*dx

    omega = .93  # 8./9#2./3.#0.83

    npmpmax = 1  # np*mp # number of subdomains at the coarsest resolution
    n1 = 8  # resolution to perform the gathering
    n0 = 8  # coarsest resolution ever

    # define the parameters
    param = {'np': np, 'mp': mp, 'nh': nh, 'n': n, 'm': m, 'npmpmax': npmpmax, 'name': name, 'omega': omega,
             'dx': dx, 'n1': n1, 'n0': n0, 'method': 'deep', 'nagglo': 2}

    gmg = Gmg(param)

    # for k in range(len(gmg.grid)):
    #     g=gmg.grid[k]
    #     zglo=combine_global(g,g.msk*1.)
    #     if myrank==0:
    #         plot2d(zglo,'zglobal')

    g = gmg.grid[0]
    mv = g.mv
    nv = g.nv
    b = zeros((mv, nv))
    x = zeros((mv, nv))

    # use the global coordinates
    xp, yp = get_global_coords(param, myrank)

    # to set up a RHS, here we pick a gaussian
    r2 = ((xp+0.3*Lx)**2+(yp-0.3*Ly)**2)
    b = exp(-r2/(2*(Lx*.4)**2))

    # don't forget to *mask* the RHS
    b = b*g.msk

    b = gmg.generate_random_large_scale_field()

    area = g.sum(b*0+1)
    sumb = g.sum(b)
    b = b-sumb/area
    # and to fill its halo
    g.halo.fill(b)

    # reset the timer of grid[0] to 0 ...
    g.set_timers_tozero()
    g.time['halo'] = 0

    x0 = x.copy()
    for k in range(10):
        x = x0.copy()
#        nite,res=gmg.solve(x,b,{'maxite':20,'tol':1e-18})

    for lev in range(gmg.nlevs):
        gmg.grid[lev].set_timers_tozero()

    x = x0.copy()
    # timing the resolution
    t0 = time.time()
#    nite,res=gmg.Vsolve(x,b,{'maxite':20,'tol':1e-8})
    nite, res = gmg.solve(x, b, {'maxite': 20, 'tol': 1e-8})
    t1 = time.time()

    # conclusion
    if myrank == 0:
        walltime = (t1-t0)
        rescaledtime = -(walltime/nite)*(g.np*g.mp) / \
            (gmg.nglo*gmg.mglo)/log(res)
        print('walltime (per ite)= %8.3e / rescaled time = %8.3e / average conv rate (per ite)= %4.2f digits / on %i x %i using %i x %i procs' %
              (walltime/nite, rescaledtime, -log(res)/nite/log(10), gmg.nglo, gmg.nglo, g.np, g.mp))

        print('TOTAL TIME : %5.3f' % (walltime))
        gmg.print_timers()

    gmg.benchmark_comm()
#    exit(0)


#    if myrank==0:
#        plot_times(gmg)
