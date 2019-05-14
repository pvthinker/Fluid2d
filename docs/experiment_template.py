#
# TEMPLATE for a fluid2d experiment
#
# with explanations of parameters and list of choices, if should be the case
#
from param import Param
from grid import Grid
from fluid2d import Fluid2d
from numpy import exp,sqrt,pi,cos,random,shape



param = Param('default.xml') # the xml file contains all the defaults values + a doc on variable
param.modelname = 'euler' # advection, euler, boussinesq, quasigeostrophic
param.expname = 'myfirstexp' # string used to name the netcdf files:
                             # expname_his.nc and expname_diag.nc



# *** domain and resolution ***
param.nx = 32*2 # nb of grid points, should of the form 2**p, 2**p*3 or 2**p*5
param.ny = 32   # idem
param.npy=4     # nb of processors in y: should be a power of 2
param.Lx = 2.   # domain size in x [Lx can be dimensional if the user
                # wishes]
param.Ly = param.Lx/2 # domain size in y. The model imposes that dx=dy
                      # => Lx/nx = Ly/ny

param.geometry = 'closed'  # closed, disc, xchannel, ychannel, perio.
# For a more complex domain, adjust the grid.msk once grid is available



# *** time ***
param.tend=4e4 # time of integration [could be dimensional of not]
param.adaptable_dt=True
param.cfl=1.   # desired cfl if adaptable_dt=True
param.dt = 1. # desired dt if adaptable_dt=False
param.dtmax=100. # max dt, used for accelerated flows or flows
                 # dominated by waves, in these case, dtmax requires
                 # a little bit of tuning



# *** discretization ***
param.order=5 # 1,2,3,4,5 : order of the spatial discretization. Even
              # cases are centered fluxes, odd cases are upwinded
              # fluxes
param.timestepping='RK3_SSP'# LFAM3, Heun, EF, LF, RK3_SSP, AB2,
                            # AB3. See core/timestepping.py to see the
                            # list, you may add your own in there



# *** output ***
param.var_to_save=['pv','psi'] # variables to save in the history. The
                               # name of variables depends on the
                               # model
param.list_diag=['ke','pv','pv2'] # diagnostics to save in the
                                  # diag. They are global
                                  # quantities. The name of the
                                  # diagnostics depends on the model

param.freq_his=50 # time interval between two histories (same unit as tend and dt)
param.freq_diag=10 # time interval between two diagnostics (same unit as tend and dt)



# *** plot ***
param.plot_var='pv' # variable to plot (depends on the model)
param.freq_plot=10 # number of time step between two plots
param.cax=[-1,1] # colorscale if colorscheme='imposed'
param.plot_interactive=False # activate the interactive plot
param.colorscheme='imposed' # 'imposed', 'minmax', 'symmetric'
param.plotting_module='plotting_adv' # name of the plotting
                                     # script. The default one is
                                     # core/plotting.py



# *** physics ***
param.beta=1.  # beta parameter (for QG model only)
param.Rd=1.    # Rossby deformation radius (for QG only)
param.gravity=1. # for Boussinesq model
param.forcing = True  # activate the forcing
param.forcing_module='forcing_dbl_gyre' #python module name where the forcing is defined
param.noslip = False # active the no-slip boundary condition (causes energy dissipation)
param.diffusion=False # activate a diffusion on tracer.
param.additional_tracer=['tracer'] # you may add an additional passive
                                   # tracer, by simply adding this
                                   # line. Example in inverse_cascade.py


grid  = Grid(param) # define the grid coordinates grid.xr, grid.yr and
                    # the mask, grid.msk. You may modify the mask just below

param.Kdiff=0.5e-4*grid.dx # diffusion coefficient (same for all
                           # tracers), there is a possibility to
                           # assign a different value for each tracer,
                           # see the experiment rayleigh_benard.py
                           # watch out, Kdiff should be adjusted when
                           # the resolution is changed





f2d = Fluid2d(param,grid) # define everything
model = f2d.model

xr,yr = grid.xr,grid.yr 
vor = model.var.get('pv') # way to access the 2D array of a variable,
                          # vor is 2 array


# set an initial tracer field
vor[:] = 1e-2*random.normal(size=shape(vor))*grid.msk
y=vor[:]*1.
model.ope.fill_halo(y) # each variable is surrounded with a halo. Halo
                       # width is param.nh=3 grid points. The halo is
                       # used in the periodic case (xchannel,
                       # ychannel, perio) and when the model is used with MPI
vor[:]=y

model.add_backgroundpv() # for the QG model

model.set_psi_from_pv()

f2d.loop() # launch the model
