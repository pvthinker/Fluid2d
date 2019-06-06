# Markus Reinert, June 2019
#
# Extensions for the EMShell, the command-line interface of the fluid2d Experiment Management System (EMS)

import numpy as np
import netCDF4 as nc
from scipy.fftpack import fft, fftshift, fftfreq


def get_strongest_wavenumber(his_filename):
    """Calculate the most intense wavenumber in y-direction.

    This function opens the given history-file, performs a Fourier
    transform in y on the masked streamfunction psi and returns the
    highest wavenumber which has at some point in time the highest
    intensity apart from the wavenumber zero."""
    # Open history file and load the data
    dataset_his = nc.Dataset(his_filename)
    ny = dataset_his.ny
    dy = dataset_his.Ly / ny
    # Save the data as a masked numpy array
    psi = dataset_his["psi"][:]
    psi.mask = 1 - dataset_his["msk"][:]
    # Set length of zero-padded signal (use ny for no zero-padding)
    fft_ny = ny
    # Caculate zero-padded Fourier-transform in y
    fft_psi = fftshift(fft(psi, n=fft_ny, axis=1), axes=1)
    # Its sampling frequency is dky = 1/dy/fft_ny
    # Calculate the corresponding axis in Fourier-space
    ky = fftshift(fftfreq(fft_ny, dy))
    # Remove the zero-frequency because it is always very large
    fft_psi[:, ky==0] = 0
    # Calculate the frequency of maximal intensity (apart from the zero frequency)
    ky_max = np.abs(ky)[np.argmax(np.max(np.abs(fft_psi), axis=2), axis=1)]
    return np.max(ky_max)
