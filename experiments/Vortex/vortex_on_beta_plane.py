"""QG-simulation of a vortex on the beta-plane.

This experiment provides an example for the usage of custom diagnostics.
Custom diagnostics are used to get the position of the vortex with a
high temporal resolution (here: daily).  This can then be used to
calculate accurately the velocity of the vortex due to the beta-effect.

To make additional information appear in the diagnostics file, we use
the parameter custom_diag which we set to a dictionary as in the example
below.  The keys of the dictionary are strings which represent the names
under which the information is stored in the NetCDF file.  The values
are functions -- typical Python style.  They take as argument the
current state of the model and return a single value, usually a float.

While the same information could be extracted later from the history
file, the advantage of using custom diagnostics is that we get a good
resolution in time with a small size of the output files.  For this, we
choose a high frequency output to write diagnostics, so we get the
information we want with very good temporal resolution, and we choose a
low frequency output for the history file, so that its file size does
not become too big.

This is extremely useful in bulk studies with the Experiment Management
System of Fluid2D, which makes it easy to run hundreds of experiments
automatically with small changes in parameters.

"""

from fluid2d import Fluid2d
from param import Param
from grid import Grid

import numpy as np


# Define custom diagnostic functions
def get_vortex_x(state):
    """Calculate the x-position of the center of a vortex.

    The implementation is for a vortex with positive relative vorticity
    and uses the global variables ‘param’ and ‘grid’.
    """
    pvanom = state[param.varname_list.index("pvanom")]
    # Take the maximum along y, then get the argmax along x and return
    # the corresponding coordinate
    return grid.x1d[np.argmax(np.max(pvanom, axis=0))]

def get_vortex_y(state):
    """Calculate the y-position of the center of a vortex.

    The implementation is for a vortex with positive relative vorticity
    and uses the global variables ‘param’ and ‘grid’.
    """
    pvanom = state[param.varname_list.index("pvanom")]
    return grid.y1d[np.argmax(np.max(pvanom, axis=1))]


param = Param("default.xml")

# Choose a name for the output of the experiment
param.expname = "QG_vortex"

# Set the type of model and domain, its size and resolution
param.modelname = "quasigeostrophic"
param.geometry = "xchannel"
param.nx = 64 * 2
param.ny = 64 * 2
# We can think of lengths as being in kilometer
param.Lx = 1000
param.Ly = 1000
# Set the number of CPU cores used
param.npx = 1
param.npy = 1

# Configure physics
param.forcing = False
param.noslip = False
param.diffusion = False
param.Rd = 100
param.beta = 0.001

# Set time settings; we can think of one time unit as one day
param.tend = 365
param.adaptable_dt = True
# for adaptable timestepping:
param.cfl = 0.8
param.dtmax = 1
# for fixed timestepping:
param.dt = 0.5

# Set the discretization
param.order = 5
param.timestepping = "RK3_SSP"
param.exacthistime = True

# Configure the output with our own diagnostics
param.var_to_save = ["psi", "pv", "pvanom"]
param.list_diag = "all"
param.custom_diag = {"x_pos": get_vortex_x, "y_pos": get_vortex_y}
# Set a high frequency for ‘diagnostics’ (every day) to have a high
# temporal resolution; set a low frequency for ‘history’ (every 10 days)
# to keep the filesize small
param.freq_his = 10
param.freq_diag = 1

# Set plot settings
param.plot_interactive = True
param.generate_mp4 = False  # requires plot_interactive to be True
param.freq_plot = 20
param.plot_psi = False
param.plot_var = "pv"
param.cmap = "plasma"

# Create model
grid = Grid(param)
f2d = Fluid2d(param, grid)
model = f2d.model

# Set the initial vorticity field
pv = model.var.get("pv")
# Get centered coordinates
x = grid.xr - param.Lx / 2
y = grid.yr - param.Ly / 2
# Create a vortex that is n-times stronger than the vorticity background
n = 2
amplitude = n * param.beta * param.Ly
width = 50
pv[:] = model.pvback + amplitude * np.exp(-(x**2 + y**2) / width**2)
model.ope.fill_halo(pv[:])
model.set_psi_from_pv()

# Start simulation
f2d.loop()
