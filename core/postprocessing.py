#!/bin/python
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt
import sys
import glob
import scipy.io as io

def plot_numvisc(diagfile):
    plt.figure()
    nc = Dataset(diagfile)
    t=nc.variables['t'][:]
    ke=nc.variables['ke'][:]
    dkdt=np.diff(ke)/np.diff(t)
    ens=nc.variables['enstrophy'][:]
    ensm=0.5*(ens[1:]+ens[:-1])
#    deltake[visc,res]=-(ke[-1]-ke[0])
        
#    deltaens[visc,res]=max(medfilt(ens,21))-ens[5]
    
    visc_tseries = -dkdt/ensm*4.4*np.pi
    visc_num = max(visc_tseries[t[1:]>0.02])
    #print('N=%4i / visc = %4.1e / num = %4.2e'%(N[res],Kdiff[visc],visc_num[res]))
    plt.semilogy(t[1:],visc_tseries)
    plt.xlabel('time')
    plt.ylabel('viscosity (-(1/2V)dE/dt)')
    plt.grid('on')
    plt.show()

def read_diags(tempdiag):
    if 1==1:
        name_dim='time'
        name_t='t_enstrophy'
        name_ke='energy'
    else:
        name_dim='t'
        name_t='t'
        name_ke='ke'
        

    diagfiles = glob.glob('%s*diag*.nc'%tempdiag)
    nd = len(diagfiles)
    print('found %i diag files',nd)

    if nd>1:
        rescaled_time = True
        print('time is considered as rescaled')
    else:
        rescaled_time = False
        print('time will be rescaled')
            
    nc = []
    diagfile=diagfiles[0]#'%s_%02i.nc'%(tempdiag,k)
    nc.append(Dataset(diagfile))
    try:
        name_dim='time'
        name_t='t_enstrophy'
        name_ke='energy'
        nt = len(nc[0].dimensions[name_dim])
        print('time dimension = %s'%name_dim)
    except:
        name_dim='t'
        name_t='t'
        name_ke='ke'
        print('time dimension = %s'%name_dim)
    nc[0].close()
    
    nc = []
    nt = 0
    for k in range(nd):
        diagfile=diagfiles[k]#'%s_%02i.nc'%(tempdiag,k)
        nc.append(Dataset(diagfile))
        nt += len(nc[k].dimensions[name_dim])

    tv=np.zeros((nt,))
    ke=np.zeros((nt,))
    ens=np.zeros((nt,))
    k0=0
    for k in range(nd):
        n = len(nc[k].dimensions[name_dim])
        tv[k0:k0+n]=nc[k].variables[name_t][:]
        ke[k0:k0+n]=nc[k].variables[name_ke][:]
        ens[k0:k0+n]=nc[k].variables['enstrophy'][:]
        nc[k].close()
        print('diag file %i : t0=%.2f - t1=%.2f'%(k,tv[k0],tv[k0+n-1]))
        k0+=n
        
    if not(rescaled_time):
        t0 = tv
        tv[1:]=np.cumsum(np.diff(tv)*np.sqrt(ens[1:]))
        #tv[1:]=np.cumsum(np.diff(tv)*ens[1:])

    L = 2*np.pi*np.sqrt( ke/ens)
    data = {'tv':tv,'ke':ke,'ens':ens,'L':L}

    return data

def diags2mat(tempdiag,data):
    filemat = '%s.mat'%tempdiag
    print('save data into %s'%filemat)
    io.savemat(filemat,data)

def plot_De(data):
    ens=data['ens']
    De = -np.diff(np.log(data['ke']))/np.diff(data['tv'])
    nu = De/np.sqrt(ens[1:])*4*np.pi * data['ke'][1:]
    Re = np.sqrt(data['ke'][1:]) / nu # assuming L=1
    #Re = data['ke'][1:] / ( np.sqrt(data['ens'][1:])* nu) # assuming L=1
    plt.figure()
    plt.semilogy(data['tv'][1:],De,label=r'$D_E$')
    plt.semilogy(data['tv'][1:],nu,label=r'$\nu$')
    plt.semilogy(data['tv'][1:],1/Re,label=r'$Re^{-1}$')
    #plt.semilogy(data['tv'][:],data['ke'][:],label=r'$E$')
    plt.xlabel(r'$t_V$')
    plt.ylabel(r'$D_E,\ \nu$')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    diagfile = sys.argv[1]
    #plot_numvisc(diagfile)
    data = read_diags(diagfile)
    plot_De(data)
    #diags2mat(diagfile,data)
