from __future__ import print_function
import numpy as np
from param import Param
from netCDF4 import Dataset
import os


class Output(Param):
    def __init__(self, param, grid, diag, flxlist=None):
        self.list_param = ['expname', 'myrank', 'nbproc',
                           'nx', 'ny', 'npx', 'npy', 'nprint',
                           'var_to_save', 'varname_list',
                           'expdir', 'tracer_list',
                           'diag_fluxes',
                           'freq_his', 'freq_diag', 'list_diag']
        param.copy(self, self.list_param)

        self.list_grid = ['nxl', 'nyl', 'nh']
        grid.copy(self, self.list_grid)

        self.param = param  # to store all the parameters in the netcdf file
        self.grid = grid
        self.flxlist = flxlist

        # prepare the 'history' file
        self.template = self.expdir+'/%s_his_%03i.nc'
        if self.nbproc > 1:
            self.hisfile = self.template % (self.expname, self.myrank)
            self.hisfile_joined = '%s/%s_his.nc' % (self.expdir, self.expname)
        else:
            self.hisfile = '%s/%s_his.nc' % (self.expdir, self.expname)

        if self.var_to_save == 'all':
            self.var_to_save = [v for v in self.varname_list]
            
        for v in self.var_to_save:
            if not(v in self.varname_list):
                raise ValueError('%s is not a model variable' % v\
                + ' => modify param.var_to_save')

        # prepare the 'fluxes' file        
        if self.diag_fluxes:
            template = self.expdir+'/%s_flx_%03i.nc'
            if self.nbproc > 1:
                self.flxfile = template % (self.expname, self.myrank)
                self.flxfile_joined = '%s/%s_flx.nc' % (self.expdir, self.expname)
            else:
                self.flxfile = '%s/%s_flx.nc' % (self.expdir, self.expname)

        # prepare the 'diagnostics' file
        if self.list_diag == 'all':
            self.list_diag = diag.keys()

        self.diagfile = '%s/%s_diag.nc' % (self.expdir, self.expname)

        self.tnextdiag = 0.
        self.tnexthis = 0.
        self.first = True
        self.grid = grid
        self.diag = diag

    def do(self, data, t, kt):
        """ write history and diag if it's timely 

        the model state 'x' and the scalar 'diag' are passed via the
        dictionnary 'data'. This allows to add new variables, e.g. 'fluxes'

        x = data['his']
        diag = data['diag']
        """
        
        if self.first:
            self.first = False
            # self.create_his(self.grid)
            self.nchis = NcfileIO(self.param, self.grid, self.hisfile,
                                  self.param.var_to_save, self.param.varname_list)
            
            self.nchis.create()
            if self.diag_fluxes:
                self.ncflx = NcfileIO(self.param, self.grid, self.flxfile, self.flxlist, self.flxlist)
                self.ncflx.create()
            if self.myrank == 0:
                self.create_diag(self.diag)

        if (t >= self.tnextdiag):
            self.tnextdiag += self.freq_diag
            if self.myrank == 0:
                self.write_diag(data['diag'], t, kt)

        if (t >= self.tnexthis):
            self.tnexthis += self.freq_his
            # self.write_his(x, t, kt)
            self.nchis.write(data['his'], t, kt)
            if self.diag_fluxes:
                self.ncflx.write(data['flx'], t, kt)


        if (self.myrank == 0) and (kt % self.nprint == 0):
            print(' / diag: %-4i - his: %-4i' %
                  (self.kdiag, self.nchis.khis), end='')

    def create_diag(self, diag):
        with Dataset(self.diagfile, 'w', format='NETCDF4') as nc:
            nc.createDimension('t', None)

            d = nc.createVariable('t', 'f', ('t',))
            d.long_name = 'model time'

            d = nc.createVariable('kt', 'i', ('t',))
            d.long_name = 'model iteration'

            for v in self.list_diag:
                d = nc.createVariable(v, 'f', ('t',))
                d.long_name = v

        self.kdiag = 0

        # set up internal buffer to avoid too frequent disk access
        self.ndiags = len(self.list_diag)+2
        self.buffersize = 10
        self.buffer = np.zeros((self.buffersize, self.ndiags))


    def write_diag(self, diag, t, kt):

        # store diag into the buffer
        k = self.kdiag % self.buffersize
        self.buffer[k, 0] = t
        self.buffer[k, 1] = kt
        j = 2
        for v in self.list_diag:
            self.buffer[k, j] = diag[v]  # getattr(diag,v)
            j += 1
        self.kdiag += 1

        # write buffer into netcdf if full
        if self.kdiag % self.buffersize == 0:
            self.dump_diag()

    def dump_diag(self):
        start = self.buffersize*((self.kdiag-1)//self.buffersize)
        last = self.kdiag-start
        k = range(start, self.kdiag)

        nc = Dataset(self.diagfile, 'r+')
        nc.variables['t'][k] = self.buffer[:last, 0]
        nc.variables['kt'][k] = self.buffer[:last, 1]
        j = 2
        for v in self.list_diag:
            nc.variables[v][k] = self.buffer[:last, j]
            j += 1
        nc.close()

    def join(self):
        if self.nbproc > 1:
            filename = self.hisfile.split('his')[0]+'his'
            join(filename)
            if self.diag_fluxes:
                filename = self.flxfile.split('flx')[0]+'flx'
                join(filename)
                
            
            
class NcfileIO(object):
    """Allow to create() and write() a Netcdf file of 'history' type,
    i.e. a set of model 2D snapshots. Variables were originally the
    model state: vorticity, psi, buoyancy etc.

    The tools are now generic and also handle the 'flux' diagnostics

    """
    def __init__(self, param, grid, ncfile, var_to_save, varname_list):
        self.ncfile = ncfile
        self.var_to_save = var_to_save
        self.varname_list = varname_list

        self.list_grid = ['nxl', 'nyl', 'nh']
        grid.copy(self, self.list_grid)

        self.param = param  # to store all the parameters in the netcdf file
        self.grid = grid

    def create(self):
        with Dataset(self.ncfile, 'w', format='NETCDF4') as nc:

            # store all the parameters as NetCDF global attributes
            dparam = self.param.__dict__
            for k in dparam.keys():
                if type(dparam[k]) == type([0, 1]):
                    # it is not straightforward to store a list in a netcdf file
                    pass

                elif type(dparam[k]) == type({}):
                    pass

                elif type(dparam[k]) == type(True):
                    if dparam[k]:
                        value = 1
                    else:
                        value = 0
                    nc.setncattr(k, value)

                else:
                    nc.setncattr(k, dparam[k])

            nc.createDimension('t', None)
            nc.createDimension('x', self.nxl-2*self.nh)
            nc.createDimension('y', self.nyl-2*self.nh)

            d = nc.createVariable('x', 'f', ('x'))
            d.long_name = 'coordinate in the x-direction (at cell center)'

            d = nc.createVariable('y', 'f', ('y'))
            d.long_name = 'coordinate in the y-direction (at cell center)'

            d = nc.createVariable('msk', 'i', ('y', 'x'))
            d.long_name = 'mask at cell center'

            d = nc.createVariable('t', 'f', ('t',))
            d.long_name = 'model time'

            d = nc.createVariable('kt', 'i', ('t',))
            d.long_name = 'model iteration'

            for v in self.var_to_save:
                d = nc.createVariable(v, 'f', ('t', 'y', 'x'))
                d.long_name = v

            # now write the mask and the coordinates
            nh = self.nh
            nc.variables['msk'][:, :] = self.grid.msk[nh:-nh, nh:-nh]
            nc.variables['x'][:] = self.grid.xr[0, nh:-nh]
            nc.variables['y'][:] = self.grid.yr[nh:-nh, 0]

        self.khis = 0

    def write(self, x, t, kt):
        nh = self.nh
        with Dataset(self.ncfile, 'r+') as nc:
            nc.variables['t'][self.khis] = t
            nc.variables['kt'][self.khis] = kt
            for v in self.var_to_save:
                iv = self.varname_list.index(v)
                nc.variables[v][self.khis, :, :] = x[iv][nh:-nh, nh:-nh]

        self.khis += 1
    
    
def join(filename):
    ''' Join history files without having to mpirun

    Useful when the run has been broken and one wants to join
    things from an interactive session '''

    template = filename+'_%03i.nc'
    hisfile_joined = filename+'.nc'

    nc0 = Dataset(template % 0, 'r')
    npx = nc0.getncattr('npx')
    npy = nc0.getncattr('npy')
    nx = nc0.getncattr('nx')
    ny = nc0.getncattr('ny')
    nxproc = nx//npx
    nyproc = ny//npy
    nbproc = npx*npy

    print(' Joining history files')
    print('-'*50)
    print('  - join his_*.nc files into a single %s' % hisfile_joined)
    print('  - found %i cores with %i x %i subdomains' % (nbproc, npx, npy))

    varlist1D = []
    varlist2D = []
    varlist3D = []
    for var in list(nc0.variables):
        v = nc0.variables[var]
        if len(v.dimensions) == 1:
            varlist1D.append(var)
        elif len(v.dimensions) == 2:
            varlist2D.append(var)
        elif len(v.dimensions) == 3:
            varlist3D.append(var)

    nt = len(nc0.dimensions['t'])  # and extract time dimension

    print('  - found %i snapshots' % nt)
    print('  - now joining the history files')

    # creates global nc file
    ncglo = Dataset(hisfile_joined, 'w', format='NETCDF4')

    ncglo.createDimension('t', None)
    ncglo.createDimension('x', nx)
    ncglo.createDimension('y', ny)

    # 1/ create the global netcdf file with all its structure

    # copy the attributes
    for att in nc0.ncattrs():
        ncglo.setncattr(att, nc0.getncattr(att))
    
    for var in varlist1D:
        vloc = nc0.variables[var]
        dim = vloc.dimensions[0]
        vglo = ncglo.createVariable(var, vloc.dtype, (dim))
        for i in range(len(vloc.ncattrs())):  # copy-paste variables attributes
            vglo.setncattr(vloc.ncattrs()[i],
                           vloc.getncattr(vloc.ncattrs()[i]))

    for var in varlist2D:
        vloc = nc0.variables[var]
        print('  - creating variable ', var)
        vglo = ncglo.createVariable(var, vloc.dtype, ('y', 'x'))
        vloc = nc0.variables[var]
        for i in range(len(vloc.ncattrs())):
            vglo.setncattr(vloc.ncattrs()[i],
                           vloc.getncattr(vloc.ncattrs()[i]))

    for var in varlist3D:
        vloc = nc0.variables[var]
        print('  - creating variable ', var)
        vglo = ncglo.createVariable(var, vloc.dtype, ('t', 'y', 'x'))
        for i in range(len(vloc.ncattrs())):
            vglo.setncattr(vloc.ncattrs()[i],
                           vloc.getncattr(vloc.ncattrs()[i]))

    ncglo.close()

    # 2/ write the data
    ncglo = Dataset(hisfile_joined, 'r+')

    for var in varlist1D:
        v = ncglo.variables[var]
        v2 = nc0.variables[var]
        d = v2.dimensions[0]
        if d == 't':
            ncglo.variables[var][:] = nc0.variables[var][:]
        if d == 'x':
            for proc in range(npx):
                ncloc = Dataset(template % proc, 'r')
                ip = proc
                ncglo.variables[var][ip*nxproc:(ip+1)*nxproc]\
                    = ncloc.variables[var][:]
                ncloc.close()
        if d == 'y':
            for proc in range(npy):
                ncloc = Dataset(template % (proc*npx), 'r')
                jp = proc
                ncglo.variables[var][jp*nyproc:(jp+1)*nyproc]\
                    = ncloc.variables[var][:]
                ncloc.close()

    nc0.close()

    for var in varlist2D:
        for proc in range(nbproc):
            ncloc = Dataset(template % proc, 'r')
            ip, jp = proc % npx, proc//npx
            ncglo.variables[var][jp*nyproc:(jp+1)*nyproc,
                                 ip*nxproc:(ip+1)*nxproc]\
                = ncloc.variables[var][:]
            ncloc.close()

    print('  - processing 3D variables')
    ncloc = []
    for proc in range(nbproc):
        ncloc.append(Dataset(template % proc, 'r'))

    for var in varlist3D:
        for kt in range(nt):
            print('\r     - %9s / kt = %i' % (var, kt), end='')
            for proc in range(nbproc):

                varname = var  # varloc.name
#                if (nbproc>16):
#                    print('\r - %s / proc = %i'%(varname,proc))
                ip, jp = proc % npx, proc//npx
                ncglo.variables[varname][kt, jp*nyproc:(jp+1)*nyproc,
                                         ip*nxproc:(ip+1)*nxproc]\
                    = ncloc[proc].variables[varname][kt]
        print()
    for proc in range(nbproc):
        ncloc[proc].close()

    ncglo.close()

    # remove subdomains files
    for k in range(nbproc):
        ncfile = template % k
        os.remove(ncfile)

    print('-'*50)
