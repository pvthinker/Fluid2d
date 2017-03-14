from __future__ import print_function
from numpy import max,zeros,array,arange,int8,intersect1d,linspace,where,median,sign,log,nan,min,max,log10,round,shape,flipud,isnan
import matplotlib.pyplot as plt
import matplotlib.text as text
from netCDF4  import Dataset
import matplotlib.animation as animation
from redblue import *

class Exp(object):
    def __init__(self,repos,expname):
        self.repos = repos
        self.expname = expname
        
        self.histemplate = repos+'/'+expname+'_%03i.nc'

        nc0		= Dataset(self.histemplate%(0),'r') # open proc 0 nc file   
        self.npx = nc0.getncattr('npx')
        self.npy = nc0.getncattr('npy')
        self.nx = nc0.getncattr('nx')
        self.ny = nc0.getncattr('ny')
        self.nt = len(nc0.dimensions['t'])
        self.nxproc          = self.nx//self.npx
        self.nyproc          = self.ny//self.npy
        self.nbproc = self.npx*self.npy
        self.x = nc0.variables['x'][:]
        self.y = nc0.variables['y'][:]
        #print('found %i cores with geometry %i x %i'%(nbproc,npx,npy))
        
        nc0.close()
        
        self.nc=[]
        for proc in range(self.nbproc):
            self.nc.append( Dataset(self.histemplate%(proc),'r'))

        self.cropped = False
        self.msk = zeros( (self.ny,self.nx), dtype=int8)
        self.x = zeros( (self.nx,) )
        self.y = zeros( (self.ny,) )
        #varname = 'msk'
        for proc in range(self.nbproc):
            #varloc  = self.nc[proc].variables[varname]
            ip,jp   = proc%self.npx,proc//self.npx
            #self.msk[jp*self.nyproc:(jp+1)*self.nyproc, \
            #        ip*self.nxproc:(ip+1)*self.nxproc] \
            #        = varloc[:,:]
            if jp==0:
                self.x[ip*self.nxproc:(ip+1)*self.nxproc] = self.nc[proc].variables['x'][:]
            if ip==0:
                self.y[jp*self.nyproc:(jp+1)*self.nyproc] = self.nc[proc].variables['y'][:]
        
        self.xidx=linspace(0,self.nx)
        self.yidx=linspace(0,self.ny)
        self.domain = [min(self.x),max(self.x),min(self.y),max(self.y)]
        print('Load experiment : %s'%(self.histemplate))
        print('(nx ,ny)  = (%i,%i)'%(self.nx,self.ny))
        print('(npx,npy) = (%i,%i)'%(self.npx,self.npy))
        print('nt = %i'%self.nt)               

    def crop(self,domain):
        self.cropped = True
        self.domain = domain
        x0,x1,y0,y1 = domain[0],domain[1],domain[2],domain[3]
        self.xidx = where( (self.x>=x0) & (self.x<=x1) )[0]
        self.yidx = where( (self.y>=y0) & (self.y<=y1) )[0]

    def read(self,varname,kt):
        if self.cropped:
            return self.read_crop(self,varname,kt)
        else:
            return self.read_full(self,varname,kt)

    def read_crop(self,varname,kt):
        z2d = zeros((len(self.yidx),len(self.xidx)))
        time = self.nc[0].variables['t'][kt]
        #print('nt=%i / kt=%i'%(self.nt,kt))
        if (kt>=self.nt):
            print('kt is out of range max(kt)=%i'%(self.nt-1))
        else:
            for proc in range(self.nbproc):
                varloc  = self.nc[proc].variables[varname]
                print('\r %12s - kt =%i - proc = %i'%(varname,kt,proc),end='')
                ip,jp   = proc%self.npx,proc//self.npx

                iiglo = arange( ip*self.nxproc,(ip+1)*self.nxproc )
                jjglo = arange( jp*self.nyproc,(jp+1)*self.nyproc )
                ii = intersect1d( self.xidx, iiglo)
                jj = intersect1d( self.yidx, jjglo)
                
                if( (len(ii)>0) & (len(jj)>0)):

                    i0,i1=ii[0]-ip*self.nxproc,ii[-1]-ip*self.nxproc+1
                    j0,j1=jj[0]-jp*self.nyproc,jj[-1]-jp*self.nyproc+1
                                            
                    zz = varloc[kt,j0:j1,i0:i1]

                    i0,i1=ii[0]-self.xidx[0], ii[-1]+1-self.xidx[0]
                    j0,j1=jj[0]-self.yidx[0], jj[-1]+1-self.yidx[0]


                    z2d[j0:j1,i0:i1] = zz

        return time,z2d
        
    def gen_anim(self,varname,cax):
        plt.ion()
        plt.figure()
        plt.set_cmap('Spectral')
        t,z2d = self.read_crop(varname,-1)
        im=plt.imshow(z2d,vmin=cax[0],vmax=cax[1],extent=self.domain)
        plt.colorbar()
        for kt in range(self.nt):
            t,z2d=self.read_crop(varname,kt)
            im.set_data(z2d)
            plt.title('t = %4.2f'%t)
            #plt.draw()
            plt.savefig('%s_%02i.png'%(varname,kt))


    def read_full(self,varname,kt):
        z2d = zeros((self.ny,self.nx))
        time = self.nc[0].variables['t'][kt]
        print('nt=%i / kt=%i'%(self.nt,kt))
        if (kt>=self.nt):
            print('kt is out of range max(kt)=%i'%(self.nt-1))
        else:
            for proc in range(self.nbproc):
                varloc  = self.nc[proc].variables[varname]
                print('\r - %s proc = %i'%(self.histemplate%(proc),proc),end='')
                ip,jp   = proc%self.npx,proc//self.npx
                z2d[jp*self.nyproc:(jp+1)*self.nyproc, \
                    ip*self.nxproc:(ip+1)*self.nxproc] \
                    = varloc[kt,:,:]
                
        return time,z2d

    def imvort(self,kt,cax=[-1,1],scaled=None,maxi=None):
        t,z2d = self.read_crop('vorticity',kt)
        t,psi = self.read_crop('psi',kt)
        z2d[z2d==0]=nan
        cm=redblue(nstripes=10)
        cm.set_bad((0.3,0.3,0.3,1))
        if scaled==None:
            im=plt.imshow(flipud(z2d),vmin=cax[0],vmax=cax[1],cmap=cm,extent=self.domain)
        else:
            idx=where(~isnan(z2d))
            print('\n',shape(idx))
            #wc = median(abs(z2d[idx[0],idx[1]]))
            wc = median(abs(z2d[idx]))
            print('wc=%g'%wc)
            z2d = flipud(z2d)
            im=plt.imshow(sign(z2d)*log(1+(z2d/wc)**2),vmin=-12,vmax=12,cmap=cm,extent=self.domain)


        plt.colorbar(im)
        if maxi == None:
            maxi = roundlog( max(abs(psi)))
        #print(linspace(-maxi,maxi,21))
        ci = maxi/10
        if scaled==None:
            plt.annotate('CI=%2g'%ci,(0.95,0.05),xycoords='axes fraction',color='g',fontsize=16,horizontalalignment='right')            
        plt.contour(self.x[self.xidx],self.y[self.yidx],psi,linspace(-maxi,maxi,21),colors='g',linewidths=2)

        #'ci=%2g'%ci
        #plt.contour(psi,[0],colors='k',linewidths=2)
        if scaled==None:
            plt.title('N = %i / t = %4.2f'%(self.nx,t),fontsize=14)
            plt.xlabel('x')
            plt.ylabel('y')
        else:
            plt.title(r'$t_v=%4.0f$'%t,fontsize=16)
            plt.xlabel(r'$x$',fontsize=16)
            plt.ylabel(r'$y$',fontsize=16)


def roundlog(x):
    k = 10**round(log10(x)) 
    return round(x/k,1)*k

def read(filename,varname,kt):
    template = filename+'_%03i.nc'

    nc0		= Dataset(template%(0),'r') # open proc 0 nc file   
    npx = nc0.getncattr('npx')
    npy = nc0.getncattr('npy')
    nx = nc0.getncattr('nx')
    ny = nc0.getncattr('ny')
    nt = len(nc0.dimensions['t'])
    nxproc          = nx//npx
    nyproc          = ny//npy
    nbproc = npx*npy
    #print('found %i cores with geometry %i x %i'%(nbproc,npx,npy))

    time = nc0.variables['t'][kt]
    nc0.close()
    z2d = zeros((ny,nx))

    print('nt=%i / kt=%i / t=%f'%(nt,kt,time))
    if (kt>=nt):
        print('kt is out of range max(kt)=%i'%(nt-1))
    else:

        for proc in range(nbproc):
            ncloc   = Dataset(template%(proc),'r')
            varloc  = ncloc.variables[varname]
            #print('\r - %s / proc = %i'%(template%(proc),proc),end='')
            ip,jp   = proc%npx,proc//npx
            z2d[jp*nyproc:(jp+1)*nyproc,ip*nxproc:(ip+1)*nxproc]\
                = ncloc.variables[varname][kt,:,:]
            ncloc.close()
    return time,z2d


def anim(filename,varname,cax):
    template = filename+'_%03i.nc'
    nc0		= Dataset(template%(0),'r') # open proc 0 nc file   
    nt = len(nc0.dimensions['t'])
    nc0.close()

    def animate(i):
        global kt
        time,z2d=read(filename,varname,kt)
        im.set_data(z2d)
        ti.set_text('%.0f     '%time)
        kt +=5
        return im,ti
    fig=plt.figure()
    ax = fig.add_subplot(111)
    time,z2d=read(filename,varname,0)
    im=ax.imshow(z2d,vmin=cax[0],vmax=cax[1])
    ti=ax.text(100,-50,'%.0f        '%time)
    print('launch the animation')
    global kt
    kt = 0
    ani = animation.FuncAnimation(fig, animate,arange(nt),interval=5, blit=True)
    plt.show()


if __name__ == "__main__":

    #plt.ion()
    fig=plt.figure()
    plt.show(False)
    fig.show()
    e=Exp('/home/roullet/data/DWall/DW1024_1e-05','dw_his')
    e.crop([0,.5,.5,1]) 
    for kt in range(e.nt):
        plt.clf()
        e.imvort(kt,cax=[-1500,1500.],maxi=0.5)     
        plt.savefig('%s/vort1_%02i.png'%(e.repos,kt))
