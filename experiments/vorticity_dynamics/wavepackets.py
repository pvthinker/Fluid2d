import numpy as np

def square(xr, yr, sigma):
    phi = np.zeros_like(xr)
    s2 = sigma/2.
    phi[(xr >= -s2) & (xr <= s2) & (yr >= -s2) & (yr <= s2)] = 1.
    return phi

def triangle(xr, yr, sigma):
    j = np.exp(1j*2*np.pi/3)
    z = xr + 1j*yr
    z = (z/sigma)

    p1x = np.real(j**0)
    p1y = np.imag(j**0)
    p2x = np.real(j**1)
    p2y = np.imag(j**1)
    p0x = np.real(j**2)
    p0y = np.imag(j**2)
    px = np.real(z)
    py = np.imag(z)

    area = 0.5*(-p1y*p2x + p0y*(-p1x + p2x) + p0x*(p1y - p2y) + p1x*p2y)
    s = 1/(2*area)*(p0y*p2x - p0x*p2y + (p2y - p0y)*px + (p0x - p2x)*py)
    t = 1/(2*area)*(p0x*p1y - p0y*p1x + (p0y - p1y)*px + (p1x - p0x)*py)

    phi = np.zeros_like(xr)
    phi[(t > 0) & (s > 0) & (1-t-s > 0)] = 1.

    return phi
