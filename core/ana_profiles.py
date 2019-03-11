import numpy as np


def vortex(xr, yr, Lx, Ly,
           x0, y0, sigma,
           vortex_type, ratio=1):

    # ratio controls the ellipticity, ratio=1 is a disc
    x = np.sqrt((xr-Lx*x0)**2+(yr-Ly*y0)**2*ratio**2)

    y = x.copy()*0.

    if vortex_type in ('gaussian', 'cosine', 'step'):
        if vortex_type == 'gaussian':
            y = np.exp(-x**2/(sigma**2))

        if vortex_type == 'cosine':
            y = np.cos(x/sigma*np.pi/2)
            y[x > sigma] = 0.

        if vortex_type == 'step':
            y[x <= sigma] = 1.
    else:
        print('this kind of vortex (%s) is not defined' % vortex_type)

    return y
