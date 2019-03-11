import numpy as np
import fortran_diag as fd


class Diag(object):
    def __init__(self, param, grid):
        self.list_param = ['npx', 'npy', 'mpi', 'list_diag']
        param.copy(self, self.list_param)

        self.list_grid = ['msk', 'nh', 'area']
        grid.copy(self, self.list_grid)

        if self.mpi == 1:
            from mpi4py import MPI

    def integrals(self, state):
        u = state.get('u')
        v = state.get('v')
        vort = state.get('vorticity')

#        u2,umax = computenormmaxu(self.msk,u,self.nh)
#        v2,vmax = computenormmaxv(self.msk,v,self.nh)

        ke, maxu = fd.computekemaxu(self.msk, u, v, self.nh)
        # ke = computekewithpsi(self.msk,psi,vort,self.nh)


        z, z2 = fd.computesumandnorm(self.msk, vort, self.nh)

#        u2 = u2 / self.area
#        v2 = v2 / self.area
        ke = ke / self.area
        z  = z / self.area
        z2 = z2 / self.area

        if self.mpi == -1:  # TODO
            cst = np.zeros((6,))
            cst[0] = u2
            cst[1] = v2
            cst[2] = z
            cst[3] = z2
            cst[4] = umax
            cst[5] = vmax

            cst_glo = np.array(MPI.COMM_WORLD.allgather(cst))

            u2, v2, z, z2, umax, vmax = (
                0., 0., 0., 0., 0., 0.)

            for k in range(self.npx*self.npy):
                u2 += cst_glo[k][0]
                v2 += cst_glo[k][1]
                z += cst_glo[k][2]
                z2 += cst_glo[k][3]
                umax = max(umax, cst_glo[k][4])
                vmax = max(vmax, cst_glo[k][5])

        # todo: add potential energy for the Boussinesq model
        self.ke = ke   # 0.5*(u2+v2)
        self.vorticity = z
        self.enstrophy = 0.5*z2
        # self.umax = umax
        # self.vmax = vmax
        self.maxspeed = maxu  # max([umax,vmax])

    def forwardbackward(self, t, dt):
        """advect by one time step then reverse time and integrate back to
        start point. The difference between the two states is due to
        the irreversibility of the flow """
        # WORK IN PROGRESS

        # backup the initial state
        x = self.var.state.copy()

        self.step(t, dt)
        self.step(t+dt, -dt)

        return self.var.state-x
