from __future__ import print_function
from netCDF4 import Dataset
import glob


class Restart(object):
    """Class to manage the execution of Fluid2d with restarts

    In the user script, instead of calling f2d.loop() use

    Restart(param, grid, f2d)

    => add 'from restart import Restart' in the header to load the class

    If no restart exists, the code will start from scratch, otherwise
    the code will start from the last restart found. Note that
    'param.tend' is interpreted as the 'time of integration' of the
    run, not the 'final time' of the run.

    restart files are overwritten during a job

    Works with mpi and several subdomains.

    Each 'his', 'diag', and 'restart' file have the restart index
    added in their name, e.g.

    vortex_00_his_000.nc : is the 'his' of rank 000 during job 0 (first)
    vortex_01_his_000.nc : is the 'his' of rank 000 during job 1 (second)
    ...

    vortex_00_restart_000.nc : is the restart created at the end of
                               job 0, to start job 1.
    ...

    output.txt is over-written by each job
    """

    def __init__(self, param, grid, f2d, launch=True):
        self.list_param = ['expname', 'expdir', 'nbproc', 'myrank',
                           'tend', 'varname_list', 'ninterrestart', 'diag_fluxes']
        param.copy(self, self.list_param)

        self.list_grid = ['nh', 'nxl', 'nyl']
        grid.copy(self, self.list_grid)

        self.template = self.expdir + '/%s_%02i_restart_%03i.nc'

        self.get_lastrestart()

        self.timelength = f2d.tend

        if self.myrank == 0:
            print('-'*50)
        if self.lastrestart is not None:
            self.restart_file = self.template % (
                self.expname, self.lastrestart, self.myrank)
            if self.myrank == 0:
                print(' Restart found')
                print(' Restarting from %s' % self.restart_file)
            tend, t, dt, kt, tnextdiag, tnexthis = self.read(f2d.model.var)
            f2d.tend += round(t)
            f2d.t = t
            f2d.dt = dt
            f2d.kt = kt
            f2d.output.tnextdiag = tnextdiag
            f2d.output.tnexthis = tnexthis

            self.nextrestart = self.lastrestart+1

        else:
            if self.myrank == 0:
                print(' No restart')
                print(' Starting from scratch')
            self.nextrestart = 0

        f2d.output.diagfile = self.expdir + '/%s_%02i_diag.nc' % (
            self.expname, self.nextrestart)
        f2d.output.template = self.expdir +'/%s_%02i_his' % (
            self.expname, self.nextrestart)+'_%03i.nc'
        f2d.output.hisfile = f2d.output.template % (self.myrank)
        f2d.output.hisfile_joined = self.expdir + '/%s_%02i_his.nc' % (
            self.expname, self.nextrestart)
        if self.diag_fluxes:
            f2d.output.template = self.expdir +'/%s_%02i_flx' % (
                self.expname, self.nextrestart)+'_%03i.nc'
            f2d.output.flxfile = f2d.output.template % (self.myrank)
            f2d.output.flxfile_joined = self.expdir + '/%s_%02i_flx.nc' % (
                self.expname, self.nextrestart)
            
            

        # split the integration in 'ninterrestart' intervals and
        # save a restart at the end of each

        self.lengthsubint = self.timelength/self.ninterrestart

        if launch:
            self.launchf2d(f2d)

    def launchf2d(self, f2d):
        """launch fluid2d and write a restart at the end

        if ninterrestart>1 do it several times. This is safer if you
        launch a long job on a cluster. If you submit a ten hours job,
        you may want to have intermediary restart points in case
        things go wrong.

        """
        for kres in range(self.ninterrestart):
            f2d.tend = f2d.t + self.lengthsubint
            f2d.loop(joinhis=(kres == self.ninterrestart-1),
                     keepplotalive=(kres < self.ninterrestart-1))
            self.restart_file = self.template % (
                self.expname, self.nextrestart, self.myrank)
            if self.myrank == 0:
                print('writing restart %i in %s' % (kres, self.restart_file))
            tend = f2d.tend
            t = f2d.t
            dt = f2d.dt
            kt = f2d.kt
            tnextdiag = f2d.output.tnextdiag
            tnexthis = f2d.output.tnexthis
            self.write(tend, t, dt, kt, tnextdiag, tnexthis, f2d.model.var)

    def get_lastrestart(self):
        """ determine the index of the last restart"""
        files = glob.glob(self.expdir + '/%s_*_restart_000.nc' % self.expname)
        nbrestart = len(files)
        if nbrestart == 0:
            self.lastrestart = None
        else:
            k = 0
            for f in files:
                pos = f.find('restart')
                idx = int(f[pos-3:pos-1])
                if idx > k:
                    k = idx
            self.lastrestart = k

    def write(self, tend, t, dt, kt, tnextdiag, tnexthis, var):
        """ write the restart file"""
        nc = Dataset(self.restart_file, 'w')

        nc.setncattr('tend', tend)
        nc.setncattr('t', t)
        nc.setncattr('dt', dt)
        nc.setncattr('kt', kt)
        nc.setncattr('tnextdiag', tnextdiag)
        nc.setncattr('tnexthis', tnexthis)

        nc.createDimension('x', self.nxl)
        nc.createDimension('y', self.nyl)

        for v in self.varname_list:
            nc.createVariable(v, 'd', ('y', 'x'))  # save in double precision
            z2d = var.get(v)
            nc.variables[v][:, :] = z2d[:, :]

        nc.close()

    def read(self, var):
        """ read the restart file"""
        nc = Dataset(self.restart_file, 'r')

        tend = nc.getncattr('tend')
        t = nc.getncattr('t')
        dt = nc.getncattr('dt')
        kt = nc.getncattr('kt')
        tnextdiag = nc.getncattr('tnextdiag')
        tnexthis = nc.getncattr('tnexthis')

        for v in self.varname_list:
            z2d = var.get(v)
            z2d[:, :] = nc.variables[v][:, :]

        nc.close()
        return tend, t, dt, kt, tnextdiag, tnexthis
