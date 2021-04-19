import numpy as np

def set_x_and_k(n, L):
    k = ((n//2+np.arange(n)) % n) - n//2
    return (np.arange(n)+0.5)*L/n, 2*np.pi*k/L

class Fourier(object):
    def __init__(self, param, grid):
        dx  =grid.dx
        dy = grid.dy
        self.nx = param.nx
        self.ny = param.ny
        self.Lx = param.Lx
        self.Ly = param.Ly
        self.nh = param.nh

        self.x, self.kx = set_x_and_k(self.nx, self.Lx)
        self.y, self.ky = set_x_and_k(self.ny, self.Ly)

        self.xx, self.yy = np.meshgrid(self.x, self.y)
        self.kxx, self.kyy = np.meshgrid(self.kx, self.ky)
        self.ktot = np.sqrt(self.kxx**2+self.kyy**2)
        self.ktot[0,0] = 1. # <- to avoid division by zero
        # shift in Fourier space to move
        # psi from cell centers to cell corners
        shift = np.exp(1j*(self.kxx*dx*0.5+self.kyy*dy*0.5))
        self.pv2psi = -(1/self.ktot)*shift
        self.pv2vor = self.ktot
        self.pv2psi[0,0] = 0.
        self.ktot[0,0] = 0.

    def invert(self, pv, psi, vor):
        nh = self.nh
        hpv = np.fft.fft2(pv[nh:-nh,nh:-nh])
        hpsi = hpv*self.pv2psi
        hvor = hpv*self.pv2vor
        psi[nh:-nh,nh:-nh] = np.real(np.fft.ifft2(hpsi))
        # Take the opportunity to compute the vorticity
        #
        # Because vorticity is a diagnostic variable
        # it should be computed only for output purpose
        # but in this version, I found it easier to
        # compute the vorticity at all calls
        vor[nh:-nh,nh:-nh] = np.real(np.fft.ifft2(hvor))
