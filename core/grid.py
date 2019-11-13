import numpy as np
from param import Param
try:
    from mpi4py import MPI
except:
    MPI = 0
from mpitools import Mpitools


class Grid(Param):

    def __init__(self, param):
        if type(MPI) == int:
            param.myrank = 0
        else:
            param.myrank = MPI.COMM_WORLD.Get_rank()

        param.nbproc = param.npx*param.npy

        param.nx = int(param.nx)
        param.ny = int(param.ny)
        self.list_param = ['nx', 'ny', 'npx', 'npy',
                           'Lx', 'Ly', 'myrank', 'nh',
                           'geometry', 'mpi',
                           'enforce_momentum', 'isisland', 'hydroepsilon']
        param.copy(self, self.list_param)
        self.debug = False

        self.mpitools = Mpitools(param)

        self.nxl = self.nx//self.npx+2*self.nh
        self.nyl = self.ny//self.npy+2*self.nh
        # local position of the rank in the matrix of subdomains
        self.i0 = self.myrank % self.npx
        self.j0 = (self.myrank // self.npx) % self.npy
        # grid size
        self.dx = self.Lx/self.nx
        self.dy = self.Ly/self.ny

        if self.dx != self.dy:
            if(self.myrank == 0):
                print('dx and dy are different')
                print('the model does not allow it')
                print('model is not yet fully validated')
            # exit(0)

        ishift = self.i0*self.nx//self.npx
        jshift = self.j0*self.ny//self.npy
        self.x1d = (np.arange(self.nxl)+0.5-self.nh+ishift)*self.dx
        self.y1d = (np.arange(self.nyl)+0.5-self.nh+jshift)*self.dy

        self.xr, self.yr = np.meshgrid(self.x1d, self.y1d)

        self.set_msk()

        if self.isisland:
            from island import Island
            self.island = Island(param, self)

        #
        # the boundary mask
        # will be done later (during operator initialization)
        # when the fill_halo will be available
        #
        # self.set_boundary_msk()

    def set_msk(self):
        """Set the mask array to default value

        according to param.geometry. The mask can be later modified by
        the user

        """
        npx = self.npx
        npy = self.npy
        nxl = self.nxl
        nyl = self.nyl
        nh = self.nh

        self.msk = np.ones((nyl, nxl), dtype=np.int8)

        if self.geometry == 'perio':
            pass

        if self.geometry in ['xperio', 'xchannel', 'closed', 'disc']:
            if self.j0 == npy-1:  # northern boundary
                self.msk[-nh:, :] = 0
            if self.j0 == 0:  # southern boundary
                self.msk[:nh, :] = 0

        if self.geometry in ['yperio', 'ychannel', 'closed', 'disc']:
            if self.i0 == npx-1:  # east
                self.msk[:, -nh:] = 0
            if self.i0 == 0:  # west
                self.msk[:, :nh] = 0

        if self.geometry == 'disc':
            r = np.sqrt((self.xr/self.Lx-0.5)**2 + (self.yr/self.Ly-0.5)**2)
            self.msk[r >= 0.5] = 0

        self.msknoslip = self.msk.copy()
        self.finalize_msk()

    def finalize_msk(self):
        """ Define grid quantities that depend on the mask

        should be after the mask array has been touched"""
        msk = self.msk
        #self.msknoslip = self.msk.copy()-self.msknoslip

        self.area = self.domain_integration(msk)

        # compute r2 for angular momentum
        x0 = self.domain_integration(self.xr*msk) / self.area
        y0 = self.domain_integration(self.yr*msk) / self.area

        self.xr0 = (self.xr - x0)*msk
        self.yr0 = (self.yr - y0)*msk

        x2 = self.domain_integration((self.xr0)**2*msk) / self.area
        y2 = self.domain_integration((self.yr0)**2*msk) / self.area
        self.x0 = x0
        self.y0 = y0
        self.x2 = x2
        self.y2 = y2
        if (self.myrank == 0) and (self.debug):
            print('domain barycenter is at (x0,y0)=(%g,%g)' % (x0, y0))
            print('domain has %i interior points' % self.area)
            print('ratio fluid/solid is  %g ' % (self.area/self.nx/self.ny))
        self.r2 = self.xr0**2 + self.yr0**2

    def domain_integration(self, z2d):
        """Define the domain integral function on grid cells"""

        nh = self.nh
        integral = np.sum(z2d[nh:-nh, nh:-nh])*1.

        integral = self.mpitools.local_to_global([(integral, 'sum')])

        return integral


if __name__ == "__main__":

    param = Param()
    param.myrank = 0
    param.npx = 2
    param.npy = 2
    param.nx = 10
    param.ny = 10
    param.geometry = 'closed'
    grid = Grid(param)

    print("myrank is :", param.myrank)
    print(grid.xr[0, :])
    print(grid.yr[:, 0])
    print(grid.msk)
    print('my coordinates in the subdomains matrix (%i,%i)'
          % (grid.j0, grid.i0))
    print('global domain area        =%f' % grid.area)
    print('global boundary perimeter =%f' % grid.bcarea)
