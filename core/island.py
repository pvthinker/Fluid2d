from numpy import ones, arange, meshgrid, sum, int8, sqrt, roll, zeros
from param import Param
from fortran_operators import celltocorner


class Island(Param):
    def __init__(self, param, grid):

        self.list_param = ['nxl', 'nyl', 'dx', 'dy']
        grid.copy(self, self.list_param)

        self.rhsp = zeros((self.nyl, self.nxl))
        self.psi = zeros((self.nyl, self.nxl))
        self.nbisland = 0
        self.data = []

    def add(self, idx, psi0):
        self.data.append({'idx': idx, 'psi0': psi0})
        self.nbisland += 1

    def finalize(self, mskp_model):
        print('found %i islands' % self.nbisland)
        mskp = zeros((self.nyl, self.nxl), dtype=int8)
        work = zeros((self.nyl, self.nxl))
        mskr = zeros((self.nyl, self.nxl))
        for k in range(self.nbisland):
            idx = self.data[k]['idx']
            psi0 = self.data[k]['psi0']
            mskr[:, :] = 1.
            mskp[:, :] = 0
            mskr[idx] = 0.
            celltocorner(mskr, work)
            mskp[work == 1] = 1
            mskp = 1-mskp

            vort = (roll(mskp, -1, axis=1)+roll(mskp, -1, axis=0)
                    + roll(mskp, +1, axis=1)+roll(mskp, +1, axis=0))

            z = (vort)*psi0/(self.dx*self.dy)  # *(1-mskp)
            self.rhsp[vort > 0] = z[vort > 0]
            self.psi[mskp == 1] = psi0
#            print(self.psi[:,10])
        print('island are ok')
