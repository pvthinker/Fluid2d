
=====
Howto
=====

.. toctree::
   :maxdepth: 2

.. _install:

Download and install
--------------------

#. Download the code from this `beta official release
   <http://pagesperso.univ-brest.fr/~roullet/fluid2d.1.2.zip>`_ (zip
   file)

   You may prefer to download it from Github and thus benefit from the
   latest update (by simply doing a git pull once you have cloned the
   package)

   https://github.com/pvthinker/Fluid2d

#. Unzip the file at a suitable place, e.g. $HOME/models/. The zip
   will create a folder `fluid2d` and place everything under it.

#. Go into fluid2d and type make. This will compile the Fortran
   routines and transform them into python modules.

   Alternatively, you may replace this step by the instructions given
   in the INSTALL file.

#. In the folder, type `source activate.sh`. This makes python aware
   of the model's path. This stage should be done everytime you open a
   new terminal. To avoid this stage, add the model's path to your
   PYTHONPATH environment variable. I recommend this latter
   approach. If you use the bash shell, this can be done by adding the
   following line in your $HOME/.bashrc

.. code::

   export PYTHONPATH=$PYTHONPATH:$HOME/models/fluid2d/core:$HOME/models/core/gmg/
..

   assuming you store fluid2d in $HOME/models

You are ready to do your first experiment.  To summarize

.. code::

   cd $HOME/mycode
   # you need to put fluid2d.1.0.zip here
   unzip fluid2d.1.0.zip
   cd fluid2d
   make
   source activate.sh

Possible issues: (tell me, I'll start a list)

Requierements
-------------

To have the code run you need to have the following python packages installed

- numpy
- matplotlib
- netCDF4
- mpi4py >= 2.0.0


Do my first experiment
----------------------

Do

.. code::

   cd fluid2d/experiments/DipoleWall
   python dipole_wall.py

and the model should run smoothly. The model will generate two NetCDF
files: `DW_0_his.nc` and `DW_0_diag.nc`. The first contains a
snapshots collection of selected models variables. The second contains
timeseries of domain integrated quantities (energy, enstrophy etc.). The prefix of the NetCDF is controlled by

.. code::

   param.expname = 'DW_0'

The selected variables are defined by

.. code::

   param.var_to_save=['vorticity','psi','v','tauw','wshear']
   param.list_diag=['ke','vorticity','enstrophy']

They list of available variables depend on the model equations.

Disable the interactive plotting
--------------------------------

The code can be run with interactive plotting turned on or off. This
is controlled by

.. code::

   param.plot_interactive = True #  or False

If you run with more than one core, each core displays its
subdomains. Windows are overlapping, it is hard to comprehend. It is
smarter to deactivate and look at the results in the NetCDF file.

Use several cpus
----------------

The code is parallelized with a domain partition. The control
parameters are `param.npx` and `param.npy`. They have to be a power of
2, not nessarily equal though. For instance to run the `myexp.py`
script with four subdomains, set

.. code::

   param.npx = 2
   param.npy = 2

and call the code with

.. code::

   mpirun -np 4 python myexp.py

Four history NetCDF are generated, one per subdomain. If you run with eight cores or less, the files are joined together automatically at the end of the experiment.

Add a passive tracer
--------------------

It is possible to enrich your experiment with one or several
tracers. It is also possible to put some gentle dynamics on it, apart
from the transport. To do it, add

.. code::

    param.additional_tracer = ['mytracer']

in the param definition section before the line

.. code::

    f2d = Fluid2d(param,grid)

Add 'mytracer' to the list of param.var_to_save if you want to have it saved in the netcdf.

Initialize the tracer distribution. For instance

.. code::

    trac = model.var.get('mytracer')
    trac[:,:] = grid.xr

to have a tracer distribution that is linear in x.

You may use a more meaningful tracer name. You may add as many tracers as you want.

You may want to have an *age tracer*
or a *dye source*. In that case look at how it is done in the :ref:`von-karman` experiment.

Restart from a past experiment
------------------------------

It is possible to split a long calculus in several jobs. For that the model need to write a `restart` file at the end of each job. And then to restart from this state for the next job. This is done by commenting the last line of the experiment file `f2d.loop()` and by adding a new line, specifically

.. code::

   #f2d.loop()
   res=Restart(param,grid,f2d)

The history and diagnostic files will have a job index in their name, starting with index 0.

Know all the possible parameters
--------------------------------

The script :download:`experiment_template.py <experiment_template.py>`
contains in front of each parameter the list of possible parameters.

For a more thorough documentation you may read the pdf file
:download:`parameters.pdf <parameters.pdf>` that is generated
automatically from the list of default parameters.

Know the default parameters
---------------------------

If a parameter is not defined in the experiment script, the model
assigns it a default value that is read from the xml file, located in
`core/defaults.xml`



Know how the variables are stored
---------------------------------

The model has a very non-standard way to store the 2D variables. It is
handled by the module `core/variables.py`. The variables are stored in
an object `var`. The 2D arrays all have the same size `[ny,nx]`. They
are packed together in a single 3D array. The first dimension is the
variable `index`. The object contains a lookup table for the
correspondance between the `index` and the variable name.  There is
only the `get` method available that returns a pointer to the 2D array
given the variable name. Example

.. code::

   vor = model.var.get('vorticity')

returns `vor` pointing on the 2D slab corresponding to the
vorticity. Since `vor` is a pointer, it is important to not overwrite
it, hence the notation

.. code::

   vor[:]-= vortex( param,grid,0.3,0.52,sigma)

where the function `vortex` returns a 2D that is assigns to the `vor` container.

The list of variables model is defined in the model's module in
the `core` directory, e.g. `core/euler.py` for the Euler model.


Know how grid coordinates are defined
-------------------------------------

All the grid information is stored in the object `grid`. The x and y
coordinates are `grid.xr` and `grid.yr`. There are cell centers
coordinates. Cells are uniformaly distributed in the rectangle whose
southwest corner is (0,0) and northeast corner is (Lx,Ly). The model
imposes the strict condition to have square cells, i.e. equal grid
size in x and y.

There is alternatively `grid.xr0` and `grid.yr0` that are shifted
coordinates with the origin at the barycenter of the domain (that
depends on the mask).

.. _howto-mask:

How to customize the mask
-------------------------

The mask is handled via the `grid.msk` array. A one indicates a fluid
cell, a zero a solid cell. Yo u can find examples of customized mask
in the :ref:`2D-turb` experiment folder.

It is important to finalize the setup by a call to

.. code::

   grid.finalize_msk()

This operation flags the boundary cells (if no-slip is activated) and computes the domain area.

Howto do a periodic channel
---------------------------

This is done by setting

.. code::

   param.geometry = 'xchannel'



.. _howto-viscosity:

How to have an explicit diffusion
---------------------------------

Simply turn on the `param.diffusion` and set the value for the
diffusion parameter `param.Kdiff`. Note that this coefficient is
likely set *after* the grid is known because you may want to scale the
diffusion coefficient as :math:`\Delta x^2`, the grid step size.

How to run an accelerated flow
------------------------------

With the Boussinesq model you may start with a rest state. The code has then no way to guess which time step to use. This is the role of

.. code::

   param.dtmax = 1.

This value will be used by default if the velocity is zero or very
small. It is your responsibility to control this value. For a
stratified flow, it is related to the Brunt-Vaisala frequency of your system.

Do a forced-dissipated experiment
---------------------------------

Look at the Rayleigh Benard experiment. Basically you need to define
the forcing in dedicated forcing file and set the
`param.forcing_module` to the name of this file.

Why is the model unstable despite a valid courant number?
---------------------------------------------------------

That's a very good question. This might be because your system
supports waves that travel faster than your current speed. In that
case, you need to figure out what are the waves of your system,
e.g. Rossby waves or internal waves, and to adjust `param.dtmax`
accordingly. This question is particularly relevant if you have a
background state (e.g. stratification or PV) and your flow has a small
amplitude.


How is the code structured?
---------------------------

Well, it tries to use object-oriented programming. The :download:`flow
chart <fluid2d_schematic.pdf>` is therefore a little bit unsual
compared to models coded in Fortran for instance. The code is stored
in the folder `core`. The subfolder `core/gmg` contains a multigrid
solver.

How to generate a movie?
------------------------

There is nothing easier. Typically you want to have something like

.. code::

   param.generate_mp4 = True
   param.freq_plot = 10
   param.adaptable_dt = False
   param.dt = 0.1 # <- experiment and resolution dependent, pick the right one
   param.colorscheme='imposed' # <- you want to have an imposed colorscale
   param.cax=[-1,1]


This will generate a `expname.mp4` file on the fly, as it is plotted
in the interactive window. By default `generate_mp4` is False. Images
update rate is controlled by `freq_plot` which means here, one frame
every 10 time iterations. To have a uniform time flow, you better
deactivate the `adaptable_dt` This means that you need to know which
time step `dt` is stable. For the colorscale `cax`, you also want to
optimize it. This step requieres some trials and errors. A caveat
though, you need to run on a single core.
