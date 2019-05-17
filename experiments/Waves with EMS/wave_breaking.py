# Markus Reinert, April/May 2019
#
# Fluid2d-experiment with EMS: wave breaking at the coast (Python file)
#
# This file aims to be both, an interesting experiment and a helpful
# introduction to the Experiment Management System (EMS) of fluid2d.
#
# To use the EMS, we import the module `ems` into the Python script of
# our experiment and create an instance of the class `EMS`, giving the
# usual Param-object and an experiment file.  This creates an entry in
# the experiment-database.  The database is stored as an SQLite file
# with the name `experiments.db` in the data directory of fluid2d, i.e.,
# the path defined in `param.datadir`.  Furthermore, it sets the value
# of `param.expname` to refer to this entry.  Therefore, when the EMS is
# used, the value of `param.expname` must not be set manually.
# Exceptions to this rule are explained below.
# At the beginning, some information is missing in the new entry of the
# database.  To complete it, we call the method `finalize` of the EMS at
# the end of the simulation.  Through this, the actual integration time
# and the size of the output files (in MB) is stored in the database.
# The field datetime is also updated and set to the current time.
# The parameters defined in the experiment file are accessible in the
# Python script through the dictionary of experiment parameters (EP).
#
# If the experiment parameters are also needed somewhere else in the
# code where the EMS-object is not accessible, for example in a forcing
# module, then a new EMS-object can be created.  This new EMS-object
# must be initialized with the same Param-object and without an
# experiment file.  In this way, both EMS-objects refer to the same
# entry of the database and share the same experiment parameters.
#
# It is possible to re-run experiments with overwriting the files, for
# example because the integration was stopped too early.
# To do this, we specify as `param.expname` the name of the experiment
# we want to run again, consisting of the name of the experiment class
# and a 3-digit number.  Then we call the constructor of EMS without an
# experiment file.  Example:
#   param.expname = "Breaking_Waves_001"
#   ems = EMS(param)
# This replaces the files in the folder `Breaking_Waves_001` and updates
# the entry 1 of the table `Breaking_Waves` at the end of the run.
#
# In general, whenever an EMS-object is created with an experiment file
# given, it creates a new entry in the database with the information
# from the experiment file and it sets `param.expname` to refer to this
# new entry.  In contrast, if an EMS-object is created without an
# experiment file, it loads the data from the database entry which
# corresponds to the value of `param.expname`.


from fluid2d import Fluid2d
from param import Param
from grid import Grid
from ems import EMS

import numpy as np


# Load default values and set type of model
param = Param("default.xml")
param.modelname = "boussinesq"

# Activate the Experiment Management System (EMS) with an experiment file
ems = EMS(param, "wave_breaking.exp")
# Fetch experiment parameters (EP)
EP = ems.get_parameters()

# Set domain type, size and resolution
param.geometry = "closed"
param.ny = 2 * 64
param.nx = 3 * param.ny
param.Ly = 2
param.Lx = 3 * param.Ly
# Set number of CPU cores used
param.npx = 1
param.npy = 1

# Set time settings
param.tend = 100.0
param.cfl = 1.2
param.adaptable_dt = False
param.dt = 0.01
param.dtmax = 0.1

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

# Set initial perturbation of the surface (or interface)
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

# Make it a complete database entry
ems.finalize(f2d)
