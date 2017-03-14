############################################################
#
# functions relative to subdomains
#
############################################################
try:
    from mpi4py import MPI
    #print('- mpi4py : found')
except:
    print('- mpi4py : not found, please install it')
    exit(0)

from numpy import zeros,arange,ones,cumsum,sqrt,array,meshgrid
from gmg.fortran_multigrid import buffertodomain
from time import time
#from plotutils import plot2d

def set_family(myrank,np,mp,ix,iy):
    procs=arange(np*mp)
    
    i = procs%np
    j = procs//np

    col=((i//ix)%2)+2*((j//iy)%2)

    if col[myrank]==0:
        rank0=myrank
    if col[myrank]==1:
        rank0=myrank-ix
    if col[myrank]==2:
        rank0=myrank-iy*np
    if col[myrank]==3:
        rank0=myrank-ix-iy*np

#    print(col[myrank])

    
    if (ix<np) and (iy<mp):
        family=array([rank0,rank0+ix,rank0+iy*np,rank0+ix+iy*np])
        family=family.reshape((2,2))
    elif (ix<np):
        if col[myrank] in (0,1):
            family=array([rank0,rank0+ix])
        else:
            family=array([rank0+iy*np,rank0+ix+iy*np])
        family=family.reshape((1,2))
    elif (iy<mp):
        if col[myrank] in (0,2):
            family=array([rank0,rank0+iy*np])
        else:
            family=array([rank0+ix,rank0+iy*np+ix])
        family=family.reshape((2,1))
    else:
        if myrank==0:            
            print('pb with family')
            print(ix,iy,np,mp)
            print(col)
        exit(0)

#    if myrank==0:
#        print('defined a new family shape %s / ix,iy=%i,%i / np,mp=%i,%i'%(family.shape,ix,iy,np,mp))
    return family

class Subdomains(object):
    def __init__(self,nh,n,m,comm,family,method=2):
        """ (n,m) is the shape of the small subdomain before gathering """
#        print('family shape=',family)
        np = family.shape[1]
        mp = family.shape[0]        

        sizes = ones(np*mp)*(n*m)
        offsets = zeros(np*mp)
        offsets[1:] = cumsum(sizes)[:-1]

        self.nh = nh
        self.n  = n
        self.m  = m
        self.family = family        
        self.np = np
        self.mp = mp
        self.n1 = 2*nh+(n-2*nh)*np
        self.m1 = 2*nh+(m-2*nh)*mp        

        self.method = method
        self.nbtimes = 1 # redundancy factor for timing purpose (should be 1 except for timing)

        myrank = comm.Get_rank()
        self.myrank = myrank

        j1,i1=(family==myrank).nonzero()
        self.i1=i1[0]*(n-2*nh)
        self.j1=j1[0]*(m-2*nh)

#        if self.myrank==0:
#            print("define buffers for family",family,"np*mp=",np*mp,"n,m=",n,m)

        self.localcomm =  MPI.COMM_WORLD.Split(family[0,0],0)       
        self.sbuff = zeros(n*m)
        self.rbuff = zeros(n*m*np*mp)

        self.sizes=sizes
        self.offsets=offsets

    def gather(self,x,y):
#        if self.myrank==0:
#            print("gather",x.shape,self.sbuff.shape,self.rbuff.shape,self.np,self.mp,self.n1,self.m1)

        for k in range(self.nbtimes):
            self.localcomm.Allgatherv(x.ravel(),
                                  [self.rbuff,self.sizes,self.offsets,MPI.DOUBLE])


            b = self.rbuff.reshape( (self.mp,self.np,self.m,self.n))
            buffertodomain(b,y,self.nh,self.m1,self.n1)

        
    def split(self,x,y):
        
        # it's important to *copy* the subarray to x to ensure
        # contiguity of the data in x
        # this instruction might work best in Fortran
        #x=zeros((self.m,self.n))
        y[:,:]=x[self.j1:self.j1+self.m,self.i1:self.i1+self.n]
        
        #print(self.m1,self.n1)
#        return x

    def timeit(self):
        nt = 100
        x=zeros( (self.m,self.n))+self.myrank

        x=zeros( (self.m,self.n))+self.myrank
        y=zeros( (self.m1,self.n1))
        z=zeros( (self.m,self.n))
        t0 = time()
        for kt in xrange(nt):
            self.gather(x,y)
        t1 = time()
        dta=t1-t0
        t0 = time()
        for kt in xrange(nt):
            subdomains.split(y,z)
        t1 = time()
        dtb=t1-t0
        return dta,dtb

    
if  __name__ =='__main__':
    comm=MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nbofproc = MPI.COMM_WORLD.Get_size()

    nh=3
    n=40
    m=80
    method = 2
    ix=1
    iy=1
    np = 2#int(sqrt(nbofproc))
    mp = nbofproc//np


    n1 = 2*nh+(n-2*nh)*np
    m1 = 2*nh+(m-2*nh)*mp        


    # if np*mp!=nbofproc:
    #     if myrank==0:
    #         print('Need exactly 4 processes to run the test, use')
    #         print('mpirun -np 4 python subdomains.py')
    #     exit()

#    family = arange(4).reshape((2,2),order='F')
    family = set_family(myrank,np,mp,ix,iy)

    if myrank==0:
        print('----------------------------------------')
        print('Family:')
        print(family)
        


    subdomains = Subdomains(nh,n,m,comm,family,method)

    if myrank==0:
        print('----------------------------------------')
        print('Time:')


    # dta,dtb=subdomains.timeit()
    # if myrank==0:
    #     print('gather : %f'%(dta))
    #     print('split  : %f'%(dtb))
    # comm.Barrier()
    

    x=zeros( (m,n))+myrank
    xp,yp=meshgrid(arange(n)*1.,arange(m)*1.)
    if myrank==0:
        x=xp.astype(float)#xp.swapaxes(0,1)
    else:
        x=x*0.
    x=yp.astype(float)

    z=zeros( (m,n))
    y=zeros((m1,n1))
    subdomains.gather(x,y)
#    subdomains.split(y,z)

    if myrank==0:
        print('----------------------------------------')
        print('y:')
#        print(y)
        print('----------------------------------------')
        print('z:')
#        print(z)
        plot2d(y,'y')
