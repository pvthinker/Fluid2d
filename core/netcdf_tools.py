from netCDF4  import Dataset
from numpy import size,array

def read_one_full_snapshot_from_netcdf(ncfile,idx):

    snapshot={}
    nc = Dataset(ncfile,'r')

    dim = list(nc.dimensions)
    
    
    timedim=''
    dims=[]
    for d in dim:
        if nc.dimensions[d].isunlimited():
            timedim = d
            length = 1
        else:
            length = nc.dimensions[d].size
        dims.append( (d,length) )

    print(dims)
            
    if not(timedim ==''):
        print('time dimension is %s',timedim)
        time = nc.variables[timedim][idx]
        print('idx %i is time = %f'%(idx,time))
    else:
        print('found no time dimension')
        print('I''m forced to abort')
        exit(0)

    variables=[]
    for var in list(nc.variables):
        if timedim in nc.variables[var].dimensions:
            value = nc.variables[var][idx]
            if size(value)==1:
                print('var %s = %f'%(var,value))
        else:
            value = nc.variables[var][idx]
            print(type(value))
            
        d=list(nc.variables[var].dimensions)
        variables.append( (var,d,value) )

    snapshot={'dimensions':dims,'variables':variables}

    return snapshot

def write_one_snaposhot(ncfile,snapshot):
    nc = Dataset(ncfile,'w')
    
    dims = snapshot['dimensions']
    variables = snapshot['variables']

    for d in dims:
        nc.createDimension(d[0],d[1])


    type_dict={"<class 'numpy.float32'>":'f', "<class 'numpy.ndarray'>":'f'}
    for v in variables:
        print('size=',size(v[2]))
        if size(v[2])>1:
            s=str(type(v[2][0]))
        else:
            s=str(type(v[2]))
        print(v[0],s)
        typ=type_dict[s]
        nc.createVariable(v[0],typ,v[1])
        if len(v[1])==1:
            nc.variables[v[0]]=v[2]
        if len(v[1])==2:
            nc.variables[v[0]][:,:]=v[2][:,:]
        if len(v[1])==3:
            nc.variables[v[0]][:,:,:]=v[2][:,:,:]
    nc.close()

if __name__ == "__main__":
    
    ncfile = '/home/roullet/project/fluid2d/practicals/turbulence/rayleighbenard_horconvhr_his.nc'

    snap = read_one_full_snapshot_from_netcdf(ncfile,-10)

#    print(snap)

    write_one_snaposhot('toto.nc',snap)
