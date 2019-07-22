############################################################
#
# functions relative to levels
#
############################################################
try:
    from mpi4py import MPI
    # print('- mpi4py : found')
except:
    print('- mpi4py : not found, please install it')
    exit(0)

from numpy import zeros, arange, sqrt, array, meshgrid, ones, concatenate, unique, random, pi, int8, size, isnan, exp, NaN, asarray
from time import time
# from fortran_multigrid import *
from gmg.subdomains import *
from gmg.halo import *
import sys
# from misc_mpi import *
# from plotutils import plot2d
# from analyticals import *


def Gridinfo(param):
    """ return the gridinfo vector
    defining a hierarchy of grids
    """
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    np = param['np']
    mp = param['mp']
    npmpmax = param['npmpmax']
    nh = param['nh']
    n = param['n']
    m = param['m']
    dx = param['dx']
    dy = param['dy']
    n1 = param['n1']  # coarsest grid
    n0 = param['n0']  # grid size for the valley (to form a peak)
    method = param['method']
    nagglo = param['nagglo']
    omega = param['omega']
    verbose = param['verbose']
    relaxation = param['relaxation']
    if 'hydroepsilon' in param.keys():
        hydroepsilon = param['hydroepsilon']
    else:
        hydroepsilon = 1
#    if param.has_key('verbose'):
#        verbose=param['verbose']

    np0 = np
    mp0 = mp

    ix = 1
    iy = 1

    gridinfo = []
    lev = 0
    ok = False
    while not(ok):
        ix0, iy0 = 1, 1
        if lev > 0:
            # by default : coarsen
            n = n//2
            m = m//2
            dx = dx*2
            dy = dy*2
            flag = 'regular'
            if((n < n1)or(m < n1))and(np*mp > 1):  # but sometimes gather
                dx = dx/2
                dy = dy/2
                flag = 'peak'
                if np > 1:
                    np = np//2
                    n = n*4
                    ix0 = 2
                else:
                    n = n*2
                if mp > 1:
                    mp = mp//2
                    m = m*4
                    iy0 = 2
                else:
                    m = m*2
        else:
            flag = 'finest'

        if(np*mp == npmpmax)and((n <= n0)or(m <= n0)):
            #        if( (n==n0)or(m==n0) ):
            ok = True
            flag = 'coarsest'
        if 'qgoperator' in param.keys():
            gridinfo.append({'n': n, 'm': m, 'np': np, 'mp': mp, 'np0': np0, 'mp0': mp0, 'omega': omega,
                             'qgoperator': True, 'Rd': param['Rd'],
                             'nh': nh, 'ix': ix, 'iy': iy, 'dx': dx, 'dy': dy, 'flag': flag, 'lev': lev,
                             'relaxation': relaxation})
        else:
            gridinfo.append({'n': n, 'm': m, 'np': np, 'mp': mp, 'np0': np0, 'mp0': mp0, 'omega': omega,
                             'qgoperator': False,
                             'nh': nh, 'ix': ix, 'iy': iy, 'dx': dx, 'dy': dy, 'flag': flag, 'lev': lev,
                             'hydroepsilon': hydroepsilon, 'relaxation': relaxation})

        if myrank == 0 and verbose:
            print("Level %2i: %5ix%5i / Procs %2ix%2i  / ix,iy=%i,%i / dx=%3.3f / Flag= %s"
                  % (lev, n, m, np, mp, ix, iy, dx, flag))
        lev += 1
        ix = ix*ix0
        iy = iy*iy0

        if lev > 20:
            if myrank == 0:
                print('too many levels, something is fishy: abort')
            exit()

    return gridinfo


class Grid(object):
    import time

    # flag,myrank,np,mp,incr,nh,n,m,comm,gather,msk0,A):
    def __init__(self, gridinfo):

        comm = MPI.COMM_WORLD
        myrank = comm.Get_rank()

        n = gridinfo['n']
        m = gridinfo['m']
        nh = gridinfo['nh']
        np = gridinfo['np']
        mp = gridinfo['mp']
        np0 = gridinfo['np0']
        mp0 = gridinfo['mp0']
        ix = gridinfo['ix']
        iy = gridinfo['iy']
        if 'hydroepsilon' in gridinfo.keys():
            self.hydroepsilon = gridinfo['hydroepsilon']
        else:
            self.hydroepsilon = 1.

        self.dx = gridinfo['dx']
        self.dy = gridinfo['dy']
        self.nh = nh
        self.np = np
        self.mp = mp
        self.np0 = np0
        self.mp0 = mp0
        self.lev = gridinfo['lev']
        self.relaxation = gridinfo['relaxation']

        aspect_ratio = self.hydroepsilon*self.dy/self.dx
        if (aspect_ratio <= .2):
            if self.mp0 > 1:
                if (myrank == 0):
                    raise ValueError(
                        'Small aspect ratio experiment requires param.npy = 1')
            else:
                if (myrank == 0) and (self.lev == 0):
                    print('Small aspect ratio detected')
                    print('Activating the tridiagonal solver for the relaxation')
                self.relaxation = 'tridiagonal'
        if (myrank == 0) and (self.lev == 0):
            print('Multigrid relaxation: %s' % self.relaxation)

        self.myrank = myrank

        # for timing purpose
        self.nbsmooth = 1
        self.nbresidual = 1
        self.nbfinetocoarse = 1
        self.nbcoarsetofine = 1

        nv = n+2*nh
        mv = m+2*nh

        self.flag = gridinfo['flag']
        self.n = n
        self.m = m
        self.nv = nv
        self.mv = mv

        self.qgoperator = gridinfo['qgoperator']
        if self.qgoperator:
            self.Rd = gridinfo['Rd']

        # finer resolution
        self.nf = n*2+2*nh
        self.mf = m*2+2*nh

        self.omega = gridinfo['omega']

        self.set_timers_tozero()

        # method for halo fill
        if np*mp == 1:
            method = 0  # fill halo with itself
#            method = -1 # fill halo with 0s
        else:
            method = 2  # use MPI communications and Fortran routines

        if self.flag == 'peak':
            family = set_family(myrank, np0, mp0, ix, iy)
            reducj = family.shape[0]
            reduci = family.shape[1]
            # coarser resolution
            self.nc = n//reduci+2*nh
            self.mc = m//reducj+2*nh
            self.subdom = Subdomains(
                nh, self.nc, self.mc, comm, family, method=method)

        neighbours = set_neighbours(myrank, np0, mp0, ix, iy)
        self.halo = Halo(nh, nv, mv, comm, neighbours,
                         method=method)  # check method

        self.msk = ones((mv, nv), dtype=int8)  # default mask full of ones

    def set_timers_tozero(self):
        # time in sec
        self.time = {'smooth': 0, 'res': 0, 'sum': 0, 'norm': 0, 'interpolate': 0,
                     'restrict': 0, 'halo': 0, 'reduce': 0, 'split': 0, 'gather': 0, 'barrier': 0}
        # number of calls
        self.ncalls = {'smooth': 0, 'res': 0, 'sum': 0, 'norm': 0, 'interpolate': 0,
                       'restrict': 0, 'halo': 0, 'reduce': 0, 'split': 0, 'gather': 0, 'barrier': 0}

    def set_msk_laplacian(self, prevlev=0, mskf=0):
        """ set the good mask and the laplacian using previous level """
        if prevlev == 0:
            self.msk = mskf
            self.compute_A_atfinest()
        else:
            mskf = prevlev.msk*1.0
            coef = ones((self.mv, self.nv))
            finetocoarse(prevlev, self, mskf, coef)
            self.msk[coef <= 0.5] = 0
            np = self.np
            mp = self.mp
            myrank = self.myrank
            iloc = myrank % np
            jloc = (myrank//np) % mp
            nh = self.nh
            self.compute_A_atcoarser(prevlev, self.msk)

        # # coefficient for the Jacobi iteration, A[:,:,4] is the main diagonal
        # val=abs(self.A[:,:,4]).max()
        # self.coef = MPI.COMM_WORLD.allreduce(val, op=MPI.MAX)
        # if self.coef!=0.:
        #     self.coef=1./self.coef
        # else:
        #     if self.myrank==0:
        #         print('matrix diagonal is zero')
        #         print('fix something!')
        #     exit()

        # buffer for 'smoothertwice' the Fortran subroutine
        self.yo = zeros((3, self.nv))

# ----------------------------------------

    def compute_A_atfinest(self):
        """ define the Laplacian on the finest grid """
        a, b, c = -4., 1., 0.
        bx = self.dy/self.dx*self.hydroepsilon
        by = self.dx/self.dy
        a = -2*(bx+by)
        stencil = array([[c, by, c], [bx, a, bx], [c, by, c]]).astype(float)

        # 9 points 2D Laplacian stencil [the factor 1/2 is weird but
        # in practice it should be here].

        if self.dx == self.dy:
            if self.myrank == 0:
                print('  => use the nine points stencil')
            a, b, c = -6./2, 1./2, 0.5/2
            stencil = array([[c, b, c], [b, a, b], [c, b, c]]).astype(float)
        else:
            if self.myrank == 0:
                print('  => use the five points stencil')

        Afinest = zeros((self.mv, self.nv, 9))
        I, J = meshgrid(arange(self.nv-2)+1, arange(self.mv-2)+1)
        I = asarray(I, dtype=int)
        J = asarray(J, dtype=int)

        coef = 1./(self.dx*self.dy)

        for k in range(9):
            Acoef = zeros((self.mv, self.nv))
            i, j = (k % 3)-1, (k//3)-1
            # check out how the mask is included here
            Acoef[J, I] = stencil[j+1, i+1]*coef * \
                self.msk[J+j, I+i]*self.msk[J, I]
#            if (k==4) and (self.qgoperator): # diagonal term => add the stretching for the QG model
#                Acoef[J,I] = Acoef[J,I] - self.msk[J,I]/(self.Rd**2)
#                print('Helmholtz operator')
            self.halo.fill(Acoef)
            Afinest[:, :, k] = Acoef
#        if self.qgoperator:
#            Afinest=-Afinest
        self.stencil = stencil
        self.A = Afinest

    def compute_A_atcoarser(self, prevlev, msk):
        """ compute the Laplacian on grid 'self' using previous grid 'prevlev' """
        self.stencil = prevlev.stencil
        if self.flag == 'peak':
            # the matrix ain't different from the previous level
            # it's simply spans a larger domain
            A = zeros((self.mv, self.nv, 9))
            for k in range(9):
                coef = prevlev.A[:, :, k].copy()
                coef2 = zeros((self.mv, self.nv))
                self.subdom.gather(coef, coef2)
                A[:, :, k] = coef2.copy()

        else:
            # let's coarsen the matrix, done in Fortran
            # we want
            # Acoarse = R * Afine * I
            A = coarsenmatrix(prevlev.A, prevlev.msk, msk, prevlev.nh)
            # don't forget to fill the halo
            for k in range(9):
                i, j = k % 3, k//3
                coef = A[:, :, k].squeeze()
                self.halo.fill(coef)
                A[:, :, k] = coef

        self.A = A
        self.stencil = prevlev.stencil

    def ninetofive(self):
        """ in this little function we take advantage that the matrix is symmetric
        so we discard the 5,6,7,8 diagonals """
        self.A = self.A[:, :, :5]


# ----------------------------------------

    def smooth(self, x, b, nite):
        if self.relaxation == 'tridiagonal':
            nite *= 3
        for l in range(self.nbsmooth):
            for k in range(nite):
                t0 = time()
#                MPI.COMM_WORLD.Barrier()
                # we apply the smoothing twice
                if self.relaxation == 'tridiagonal':
                    smoothtridiag(self.msk, self.A, x, b)
                else:
                    smoothtwicewitha(self.msk, self.A, x,
                                     b, self.omega, self.yo)

#                smoothoncewitha(self.msk,self.A,x,b,self.omega,self.yo)
#                smoothoncewitha(self.msk,self.A,x,b,self.omega,self.yo)
                t1 = time()
                self.time['smooth'] += t1-t0
                self.ncalls['smooth'] += 1

#                MPI.COMM_WORLD.Barrier()
                t0 = time()
                self.time['barrier'] += t0-t1
                self.ncalls['barrier'] += 1

                self.halo.fill(x)
                t1 = time()
                self.time['halo'] += t1-t0
                self.ncalls['halo'] += 1

    def residual(self, x, b, r):
        for l in range(self.nbresidual):
            t0 = time()
#            MPI.COMM_WORLD.Barrier()
            computeresidualwitha(self.msk, self.A, x, b, r)
            t1 = time()
            self.time['res'] += t1-t0
            self.ncalls['res'] += 1

#            MPI.COMM_WORLD.Barrier()
            t0 = time()
            self.time['barrier'] += t0-t1
            self.ncalls['barrier'] += 1

            self.halo.fill(r)
            t1 = time()
            self.time['halo'] += t1-t0
            self.ncalls['halo'] += 1

    def norm(self, x):
        """ norm = sum(x*x) """
        nbduplicates = (self.np0*self.mp0)/(self.np*self.mp)
        # computenorm is done in Fortran
        self.typenorm = 'l2'
        t0 = time()
#        MPI.COMM_WORLD.Barrier()
        if self.typenorm == 'l2':
            local_sum = computenorm(self.msk, x, self.nh)
            t1 = time()
            self.time['norm'] += t1-t0
            self.ncalls['norm'] += 1
            z = MPI.COMM_WORLD.allreduce(local_sum, op=MPI.SUM) / nbduplicates
            z = sqrt(z)
            t0 = time()
            self.time['reduce'] += t0-t1
            self.ncalls['reduce'] += 1

        if self.typenorm == 'inf':
            local_z = computemax(self.msk, x, self.nh)
            t1 = time()
            self.time['norm'] += t1-t0
            self.ncalls['norm'] += 1
            z = MPI.COMM_WORLD.allreduce(local_z, op=MPI.MAX)
            t0 = time()
            self.time['reduce'] += t0-t1
            self.ncalls['reduce'] += 1
        return z

    def inner(self, x, y):
        """ inner = sum(x*y) """
        nbduplicates = (self.np0*self.mp0)/(self.np*self.mp)
        # computenorm is done in Fortran

        local_sum = computeinner(self.msk, x, y, self.nh)
        z = MPI.COMM_WORLD.allreduce(local_sum, op=MPI.SUM) / nbduplicates

        return z

    def sum(self, x):
        """  sum(x) """
        t0 = time()
#        MPI.COMM_WORLD.Barrier()
        nbduplicates = (self.np0*self.mp0)/(self.np*self.mp)
        local_sum = computesum(self.msk, x, self.nh)
        t1 = time()
        self.time['sum'] += t1-t0
        self.ncalls['sum'] += 1
        global_sum = array(MPI.COMM_WORLD.allgather(local_sum))
        t0 = time()
        self.time['reduce'] += t0-t1
        self.ncalls['reduce'] += 1

        return global_sum.sum() / nbduplicates

# ----------------------------------------


def coarsetofine(coarse, fine, x, y):
    """ input = x is on coarse / output = y is on fine """
    # these are function outside the object Grid
    # because it works on two different grids ...
    for k in range(coarse.nbcoarsetofine):
        if coarse.flag == 'peak':
            t0 = time()
            coarse.subdom.split(x, y)
            t1 = time()
            coarse.time['split'] += t1-t0
            coarse.ncalls['split'] += 1
        else:
            t0 = time()
            interpolate(fine.msk, coarse.msk, x, fine.nh, y)
            t1 = time()
            fine.time['interpolate'] += t1-t0
            fine.ncalls['interpolate'] += 1
            if coarse.nh <= 2:
                fine.halo.fill(y)
                t0 = time()
                fine.time['halo'] += t0-t1
                fine.ncalls['halo'] += 1
            # else there is enough points in the halo to skip filling!
#        MPI.COMM_WORLD.Barrier()


def finetocoarse(fine, coarse, x, y):
    for k in range(fine.nbfinetocoarse):
        if coarse.flag == 'peak':
            t0 = time()
            coarse.subdom.gather(x, y)
            t1 = time()
            coarse.time['gather'] += t1-t0
            coarse.ncalls['gather'] += 1
        else:
            t0 = time()
            restrict(coarse.msk, x, coarse.nh, y)
            t1 = time()
            coarse.time['restrict'] += t1-t0
            coarse.ncalls['restrict'] += 1

#            MPI.COMM_WORLD.Barrier()
            t0 = time()
            coarse.time['barrier'] += t0-t1
            coarse.ncalls['barrier'] += 1

            coarse.halo.fill(y)
            t1 = time()
            coarse.time['halo'] += t1-t0
            coarse.ncalls['halo'] += 1

# ----------------------------------------


def combine_global(grid, x):
    """ return the global x array living on 'grid' """

    comm = MPI.COMM_WORLD

    nh = grid.nh
    nv = grid.nv
    mv = grid.mv
    N = nv*mv
    np = grid.np0
    mp = grid.mp0
    n = grid.n
    m = grid.m

    di, dj = n, m
    di, dj = nv+2, mv+2

    nv0 = di*np
    mv0 = dj*mp

    sizes = ones(np*mp)*N
    offsets = arange(np*mp).reshape((mp, np)).ravel()*N

    buff_loc = x.ravel().copy()
    buff_loc[grid.msk.ravel() == 0] = NaN
    buff_glo = zeros(N*np*mp)

    comm.Allgatherv(buff_loc, [buff_glo, sizes, offsets, MPI.DOUBLE])

    I, J = meshgrid(arange(nv), arange(mv))

    xglo = zeros((mv0, nv0))
    for j in range(mp):
        for i in range(np):
            k = i+j*np
            xglo[J+j*dj, I+i*di] = buff_glo[k*N:(k+1)*N].reshape((mv, nv))

    return xglo

# ----------------------------------------
# to check how the matrix looks like on the various grids
# interesting!


def combinematrix(grid):
    """ return the matrix (9,n,m) as a global (3n,3m) array """

    nh = grid.nh
    nv = grid.nv
    mv = grid.mv

    I, J = meshgrid(arange(nv), arange(mv))

    xglo = zeros((mv*3, nv*3))
    for k in range(9):
        i, j = k % 3, k//3
        coef = grid.A[:, :, k].reshape((mv, nv))
        xglo[J+j*mv, I+i*nv] = coef
    return xglo


# ----------------------------------------
# for debugging
def checkfornan(x):
    return size(isnan(x.ravel()).nonzero()[0])


def checkstencil(grid):
    for k in range(9):
        i, j = k % 3, k//3
        dx = grid.dx
        c = grid.A[:, :, k]*(dx**2)
        print(c[4, 4])

        c = c[c != grid.stencil[j, i]]
        nb = size(c)


# ----------------------------------------
if __name__ == '__main__':
    from plotutils import plot2d

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nbofproc = MPI.COMM_WORLD.Get_size()

    np = int(sqrt(nbofproc))
    mp = np

#    mp = 2
#    np = nbofproc/mp

    n = 64
    m = 64
    nh = 3

    lev = 0
    ix = 1
    iy = 1
    gather = 0
    A = 0

    # level=Level('finest',myrank,np,mp,incr,nh,n,m,comm,gather,msk0,A)

    # if myrank==0:
    #     print('----------------------------------------')
    #     print(level.A)
    #     print('----------------------------------------')
    #     print(level.msk1)

    n = 128//np
    m = 128//mp

    dx = 1.  # /(n*np)
    # define the parameters
    param = {'np': np, 'mp': mp, 'nh': nh, 'n': n, 'm': m, 'omega': 0.85,
             'dx': dx, 'dy': dx, 'n1': 32, 'n0': 8, 'method': 'deep', 'nagglo': 2}

    # setup the gridinfo vector
    info = Gridinfo(param)
    nlevs = len(info)

    # define the grids hierarchy
    grid = []
    for lev in range(nlevs):
        grid.append(Grid(info[lev]))

    msk = zeros((m+2*nh, n+2*nh), dtype=int8)+1

    name = 'circular'
#    name = 'perio'
    msk = set_fine_msk(param, myrank, name)

    if myrank == 0:
        print('----------------------------------------')
        print(' domain is %s' % (name))
    for lev in range(nlevs-1):
        comm.Barrier()
        x = zeros((grid[lev].mv, grid[lev].nv))+myrank
        xp, yp = meshgrid(arange(grid[lev].nv)*1., arange(grid[lev].mv)*1.)
        x = xp
        grid[lev].halo.fill(x)
        if myrank == 0:
            print(grid[lev].nv, grid[lev].mv, grid[lev-1].nv, grid[lev-1].mv)
            sys.stdout.write('sending x from lev %i to %i: ' % (lev, lev+1))
        y = zeros((grid[lev+1].mv, grid[lev+1].nv))
        finetocoarse(grid[lev], grid[lev+1], x, y)
        try:
            y = zeros((grid[lev+1].mv, grid[lev+1].nv))
            finetocoarse(grid[lev], grid[lev+1], x, y)
            if myrank == 0:
                sys.stdout.write(' ok\n')
            # xglo=combine_global(grid[lev],x)

            # if myrank==0:
            #     l=lev
            #     plot2d(xglo,'y / lev=%i / %s / n*m=%ix%i'%(l,grid[l].flag,grid[l].n,grid[l].m))
            # del xglo
        except:
            if myrank == 0:
                print('failed')
            exit(0)
        # if myrank==0:
        #     plot2d(y,'y / lev=%i'%(lev))

    if myrank == 0:
        print('----------------------------------------')
    for lev in range(nlevs-1, 0, -1):
        comm.Barrier()
        x = zeros((grid[lev].mv, grid[lev].nv))+myrank
        xp, yp = meshgrid(arange(grid[lev].nv)*1., arange(grid[lev].mv)*1.)
        x = yp
        grid[lev].halo.fill(x)

        if myrank == 0:
            print(grid[lev].nv, grid[lev].mv, grid[lev-1].nv, grid[lev-1].mv)
            sys.stdout.write('sending x from lev %i to %i: ' % (lev, lev-1))
        comm.Barrier()
        y = zeros((grid[lev-1].mv, grid[lev-1].nv))

        coarsetofine(grid[lev], grid[lev-1], x, y)

        try:
            y = zeros((grid[lev-1].mv, grid[lev-1].nv))
            coarsetofine(grid[lev], grid[lev-1], x, y)
            if myrank == 0:
                sys.stdout.write(' ok\n')

        except:
            if myrank == 0:
                print('failed')
            exit(0)
#         if myrank==0:
# #            z[0,0]=-1.
# #            print(y)
#             plot2d(y,'y / lev=%i'%(lev-1))
        # xglo=combine_global(grid[lev],x)
        # if myrank==0:
        #     plot2d(xglo,'x')

    comm.Barrier()

    # populate them
    for lev in range(nlevs):
        comm.Barrier()
        if myrank == 0:
            print('---------------- level %i ------------------------' % (lev))
        if lev == 0:
            grid[lev].set_msk_laplacian(prevlev=0, mskf=msk)
        else:
            grid[lev].set_msk_laplacian(prevlev=grid[lev-1])
        comm.Barrier()
        checkstencil(grid[lev])
        if myrank == 0:
            print('matrix on grid %i is ok' % (lev))

        A = grid[lev].A
        if myrank == 0:
            print('lev %i / dx = %f' % (lev, grid[lev].dx))
            print(unique(A))

        # check for NaNs
        nbnan = checkfornan(A)
        if nbnan != 0:
            print('found NaNs in A!')
            exit()

        if myrank == 0:
            xglo = combinematrix(grid[lev])
        #     try:
#            plot2d(xglo,'matrix / lev=%i / proc=%i'%(lev,myrank))
        #     except:
        #         print('pb with plot matrix / lev=%i'%(lev))
        #         print(xglo.shape)

    for lev in range(nlevs):
        grid[lev].ninetofive()

    # transform the Laplacian into a Helmholtz operator if requiered
    if self.qgoperator:
        for lev in range(nlevs):
            print('****************Helmholtz operator')
            grid[lev].A[:, :, 4] -= 1./(self.Rd**2)

    # test soundness of laplacian

    # test smoothing
    for lev in range(nlevs):
        g = grid[lev]
        if myrank == 0:
            print('----------------------------------------')
            print('lev %i' % (lev))

        if g.flag != 'peak':
            nv = g.nv
            mv = g.mv
            nh = g.nh

            area = g.sum(g.msk*1.)
#            b=random.uniform(-1,1,(nv,mv))*g.msk
            b = zeros((mv, nv))
            xp, yp = meshgrid(arange(nv), arange(mv))
#            xp,yp=get_global_coords(param,myrank)
#            b=xp
            x = zeros((mv, nv))
            #            x[:,:]=xp.astype(float)*g.msk
            #            x=sin(2*pi*xp/(nv-2*nh))
            i0 = n/(2*2**lev)
            j0 = m/(2*2**lev)
            if myrank == 0:
                b[j0, i0/3] = 1.
                b[j0, i0*2/3] = -1.

            r2 = ((xp*1./nv)-0.35)**2+((yp*1./mv)-0.45)**2
            r3 = ((xp*1./nv)-0.45)**2+((yp*1./mv)-0.55)**2
            b = exp(-r2*50)-exp(-r3*50)
            b = (b-g.sum(b)/area)*g.msk
            b = b*g.msk
            g.halo.fill(b)
            r = zeros((mv, nv))

            g.residual(b, x, r)
#            zglo=combine_global(g,r)
            print('sum r=', g.sum(r))
#            if myrank==0:
#                plot2d(zglo,'r')

            normeb = g.norm(b)
            if myrank == 0:
                print('    normeb=%g' % (normeb))
            x0 = x.copy()
            for k in range(14):
                comm.Barrier()
                g.smooth(x, b, g.n/2)
#                print('sum x=',g.sum(x))
#                if myrank==0:
#                    print(x)
                comm.Barrier()
                nbnan = checkfornan(x)
                if nbnan != 0:
                    print('found NaNs in x!')
                    print(x)
                    exit()
                nbnan = checkfornan(r)
                if nbnan != 0:
                    print('found NaNs in r!')
                    exit()

                g.residual(x, b, r)
                zglo = combine_global(g, r)
                if (myrank == 0) and (k % 3) == 1:
                    plot2d(zglo, 'r /lev=%i / ite=%i' % (lev, k))

                res = g.norm(r)/normeb
                x0 = x.copy()
                if myrank == 0:
                    print('    res=%g' % (res))

            if lev > 0:
                r1 = zeros((grid[lev-1].mv, grid[lev-1].nv))
                coarsetofine(grid[lev], grid[lev-1], r, r1)
                zglo = combine_global(grid[lev-1], r1)
#                zglo=combine_global(g,x)
                if myrank == 0:
                    plot2d(zglo, 'r')
        else:
            if myrank == 0:
                print('  no smoothing ')
