from __future__ import print_function
from numpy import max,zeros,array,round
from netCDF4  import Dataset
import os
import glob

class Restart(object):
    def __init__(self,param,grid,f2d,launch=True):
        self.list_param=['expname','nbproc','myrank','tend','varname_list','ninterrestart']
        param.copy(self,self.list_param)

        self.list_grid=['nh','nxl','nyl']
        grid.copy(self,self.list_grid)

        self.template = '%s_%02i_restart_%03i.nc'

        self.get_lastrestart()
                
        self.timelength = f2d.tend

        if self.lastrestart!=None:
            self.restart_file = self.template%(self.expname,self.lastrestart,self.myrank)
            if self.myrank==0:
                print('restarting from %s'%self.restart_file)
            tend,t,dt,kt,tnextdiag,tnexthis = self.read(f2d.model.var)
            f2d.tend += round(t)
            f2d.t  = t
            f2d.dt = dt
            f2d.kt = kt
            f2d.output.tnextdiag = tnextdiag
            f2d.output.tnexthis  = tnexthis

            self.nextrestart = self.lastrestart+1

        else:
            if self.myrank==0:
                print('starting from scratch')
            self.nextrestart = 0

        f2d.output.diagfile = '%s_diag_%02i.nc'%(self.expname,self.nextrestart)
        f2d.output.template = '%s_%02i_his'%(self.expname,self.nextrestart)+'_%03i.nc'
        f2d.output.hisfile = f2d.output.template%(self.myrank)
        f2d.output.hisfile_joined = '%s_%02i_his.nc'%(self.expname,self.nextrestart)

        # split the integration in 'ninterrestart' intervals and
        # save a restart at the end of each
        #self.ninterrestart = 10
        self.lengthsubint = self.timelength/self.ninterrestart

        if launch:            
            self.launchf2d(f2d)



    def launchf2d(self,f2d):

        for kres in range(self.ninterrestart):
            f2d.tend = f2d.t + self.lengthsubint
            f2d.loop(joinhis=(kres==self.ninterrestart-1))
            self.restart_file = self.template%(self.expname,self.nextrestart,self.myrank)
            if self.myrank==0:
                print('writing restart %i in %s'%(kres,self.restart_file))
            tend  = f2d.tend  
            t  = f2d.t  
            dt = f2d.dt 
            kt = f2d.kt 
            tnextdiag  = f2d.output.tnextdiag
            tnexthis   = f2d.output.tnexthis
            self.write(tend,t,dt,kt,tnextdiag,tnexthis,f2d.model.var)


    def get_lastrestart(self):
        files = glob.glob('%s_*_restart_000.nc'%self.expname)
        nbrestart = len(files)
        if nbrestart == 0:
            self.lastrestart=None
        else:
            k=0
            for f in files:
                pos = f.find('restart')
                idx = int(f[pos-3:pos-1])
                if idx>k: k = idx
            self.lastrestart = k

    # def write_param(self,param):
    #     nc = Dataset(self.restart_file,'r+')
    #     for atr in dir(param):
    #         if atr[0]!='_':
    #             nc.setncattr(atr,getattr(param,atr))
    #     nc.close()

    # def read_param(self,param):
    #     nc = Dataset(self.restart_file,'r')
    #     for atr in nc.ncattrs():
    #         print(atr)
    #         setattr(param,atr,nc.getncattr(atr))
    #     nc.close()

    def write(self,tend,t,dt,kt,tnextdiag,tnexthis,var):
        nh = self.nh
        nc = Dataset(self.restart_file,'w')

        nc.setncattr('tend',tend)
        nc.setncattr('t',t)
        nc.setncattr('dt',dt)
        nc.setncattr('kt',kt)
        nc.setncattr('tnextdiag',tnextdiag)
        nc.setncattr('tnexthis',tnexthis)

        nc.createDimension('x',self.nxl)
        nc.createDimension('y',self.nyl)

        for v in self.varname_list:
            nc.createVariable(v,'d',('y','x')) # save in double precision
            z2d = var.get(v)
            nc.variables[v][:,:]=z2d[:,:]

        nc.close()             

    def read(self,var):
        nh = self.nh
        nc = Dataset(self.restart_file,'r')

        tend  = nc.getncattr('tend')
        t  = nc.getncattr('t')
        dt = nc.getncattr('dt')
        kt = nc.getncattr('kt')
        tnextdiag = nc.getncattr('tnextdiag')
        tnexthis  = nc.getncattr('tnexthis')

        y = zeros((self.nyl,self.nxl))
        for v in self.varname_list:
            z2d = var.get(v)
            z2d[:,:] = nc.variables[v][:,:]

        nc.close()  
        return tend,t,dt,kt,tnextdiag,tnexthis
        
        
# def create_restart_from_his(hisfile):
#     # get information from the core 0 history file

#     # loop over history files and create the associated restart
#     for proc in range(nbproc):
        


#     def create_from_his(self):

# #        template = self.get_template_from_netcdf(self.history)
# #        self.create_netcdf_from_template(self.restart,template)    

#         snapshot = self.read_one_snapshot_from_netcdf(self.history)
#         self.write_one_snaposhot(self.restart,snapshot)

# # snapshot contains both the data and the metadata

#     def read(self,x):
#         nc = Dataset(self.restartfile,'r')
#         t  = nc.variables['t'][-1]
#         kt = nc.variables['kt'][-1]
#         nh = self.nh
#         y = zeros((nyl,nxl))

#         for v in self.var_to_save:
#             iv = self.varname_list.index(v)
#             y[nh:-nh,nh:-nh] = nc.variables[v][-1,:,:]
#             self.fill_halo(y)
#             x[iv] = y
            
#         nc.close()
#         return t,kt
