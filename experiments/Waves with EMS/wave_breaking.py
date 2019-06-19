"""Fluid2d-experiment with EMS: wave breaking at the coast (Python file)

This file, together with its corresponding experiment file, aims to be
both, an interesting experiment and a helpful introduction to the
Experiment Management System (EMS) of fluid2d.

Compared to a simple fluid2d experiment without EMS, three changes are
necessary to profit from the full functionality of the EMS:

 1. Create Param by specifying the path to the experiment file.
 2. Fetch and use the experiment parameters from Param.
 3. Call the finalize-method of Param at the end.

Furthermore, the variable `param.expname` should not be set in the
Python file, since the EMS takes care of setting it.  In this way, the
loss of old output files by accidentally replacing them is avoided.
Nevertheless, it is possible to explicitly ask the EMS to replace old
experiments.  See the experiment file for more details.

When the Experiment Management System is activated like this, an entry
in the experiment-database is created.  The database is stored in the
same directory as the output of the experiments, which is defined in
`param.datadir`.  Every entry contains the unique ID of the experiment,
the date and time of the run, the last point in time of the integration,
the sizes of the output files in MB, a comment or description and
-- most importantly -- the parameters that are chosen by the user to
keep track off.  These parameters are defined in the experiment file.
After the EMS is set-up, only the experiment file needs to be changed by
the user.  Thanks to the EMS, it is not necessary to modify the Python
file more than once at the beginning to set it up.

Author: Markus Reinert, April/May/June 2019
"""

from fluid2d import Fluid2d
from param import Param
from grid import Grid

import numpy as np


### STEP 1
# Load default parameter values, load specific parameters from the given
# experiment file and set-up the EMS.
param = Param(ems_file="wave_breaking.exp")

# Do not set `param.expname` because it will be ignored.

### STEP 2
# There are two ways to get the dictionary with the experiment parameters.
# (a) The quick and easy way (using "get"):
# To run exactly one experiment with only the values specified at first, use
#   EP = param.get_experiment_parameters()
# In this case, it is not necessary to indent the whole file as below.
# (b) The recommended way (using "loop"):
# To run one experiment for every possible combination of all the values
# specified in the experiment file, use
#   for EP in param.loop_experiment_parameters():
# and put the rest of the Python file within this loop by indenting every line.
# The behaviour of both ways is the same if only one value is given for every
# parameter in the experiment file.  Therefore it is recommended to use the
# multi-run implementation, since it is more versatile.
#
# In both cases, the dictionary EP is now available.  Its keys are the names
# of the parameters in the experiment file.  Every key refers to the
# corresponding value.  If multiple values are given in the experiment file,
# in every iteration of the loop, another combination of these values is in EP.
# Within the following lines of Python code, the values of the EP are used to
# implement the desired behaviour of the experiment.
for EP in param.loop_experiment_parameters():
    # Set model type, domain type, size and resolution
    param.modelname = "boussinesq"
    param.geometry = "closed"
    param.ny = 2 * 64
    param.nx = 3 * param.ny
    param.Ly = 2
    param.Lx = 3 * param.Ly
    # Set number of CPU cores used
    # Remember to set a fixed (initial) ID in the experiment file, if multiple
    # cores are used.
    param.npx = 1
    param.npy = 1

    # Use a fixed time stepping for a constant framerate in the mp4-file
    param.tend = 20.0
    param.adaptable_dt = False
    param.dt = 0.02

    # Choose discretization
    param.order = 5

    # Set output settings
    param.var_to_save = ["vorticity", "buoyancy", "psi"]
    param.list_diag = "all"
    param.freq_his = 0.2
    param.freq_diag = 0.1

    # Set plot settings
    param.plot_var = "buoyancy"
    param.plot_interactive = True
    param.generate_mp4 = True
    param.freq_plot = 10
    param.colorscheme = "imposed"
    param.cax = [0, 1]
    param.cmap = "Blues_r"  # reversed blue colour axis

    # Configure physics
    param.gravity = 1.0
    param.forcing = False
    param.noslip = False
    param.diffusion = EP["diffusion"]
    param.Kdiff = EP["Kdiff"] * param.Lx / param.nx

    # Initialize geometry
    grid = Grid(param)
    xr, yr = grid.xr, grid.yr

    LAND = 0
    AQUA = 1

    # Add linear sloping shore
    m = EP["slope"]
    t = EP["height"] - EP["x_start"] * m
    grid.msk[(yr <= m*xr + t) & (yr < EP["height"])] = LAND
    grid.finalize_msk()

    # Create model
    f2d = Fluid2d(param, grid)
    model = f2d.model

    # Set initial perturbation of the interface
    buoy = model.var.get("buoyancy")
    buoy[:] = 1
    if EP["perturbation"].lower() == "sin":
        # Use a sinusoidal perturbation
        buoy[
            (yr < EP["intensity"] * np.sin(2 * np.pi * xr/EP["sigma"]) + 1.0)
            & (xr < EP["x_start"])
        ] = 0
    elif EP["perturbation"].lower() == "gauss":
        # Use a Gaussian perturbation
        buoy[
            (yr < EP["intensity"] * np.exp(-(xr/EP["sigma"])**2) + 1.0)
            & (xr < EP["x_start"])
        ] = 0
    else:
        raise ValueError("unknown type of perturbation: {}.".format(EP["perturbation"]))
    buoy *= grid.msk

    # Start simulation
    f2d.loop()

    ### STEP 3
    # Finish off the EMS database entry.
    # Without this call to `finalize`, the size of the output files cannot be
    # saved in the database.  It is important to have this line within the
    # for-loop if the multi-run strategy (b) is used.
    param.finalize(f2d)
