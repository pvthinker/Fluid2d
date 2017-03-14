############################################################
#
# functions relative to halo
#
############################################################
try:
    from mpi4py import MPI
    #print('- mpi4py : found')
except:
    print('- mpi4py : not found, please install it')
    exit(0)
       
from numpy import zeros,arange,sqrt,int32,meshgrid,where
from time import time
from gmg.fortran_multigrid import *
import numpy 
#from fillhalo import fillhalo_fortran


def set_neighbours(myrank,np,mp,ix,iy):
    i=myrank%np
    j=myrank//np
    im=(i-ix)%np
    jm=(j-iy)%mp
    ip=(i+ix)%np
    jp=(j+iy)%mp

    sw = im+jm*np
    s  = i +jm*np
    se = ip+jm*np

    w  = im+j*np
    e  = ip+j*np

    nw = im+jp*np
    n  = i +jp*np
    ne = ip+jp*np


    return [sw,s,se,w,e,nw,n,ne]
#    return [sw,w,nw,s,n,se,e,ne]



class Halo(object):
    def __init__(self,nh,n,m,comm,neighbours,method=0):
        self.n = n # number of in x, including halo
        self.m = m # number of in y, including halo
        self.nh = nh # halo width
        self.method = method
        self.comm = comm
        self.myrank = self.comm.Get_rank()
        self.nbhalo = 1 # redundancy factor for timing purpose (should be 1 except for timing)

        self.neighbours = neighbours

        n2 = n-2*nh # inner point = without halo
        m2 = m-2*nh
        self.shapes=[(nh,nh),(nh,n2),(nh,nh),(m2,nh),(m2,nh),(nh,nh),(nh,n2),(nh,nh)]
#        self.shapes=[(nh,nh),(m2,nh),(nh,nh),(nh,n2),(nh,n2),(nh,nh),(m2,nh),(nh,nh)]

        # preallocated buffers 
        self.sbuff=[]
        self.rbuff=[]
        for k in range(8):
            self.sbuff.append( zeros(self.shapes[k]))
            self.rbuff.append( zeros(self.shapes[k]))

        self.reqs = []
        self.reqr = []
        for k in range(8):
            self.reqs.append( comm.Send_init(self.sbuff[k],neighbours[k], k ) )
            self.reqr.append( comm.Recv_init( self.rbuff[k],neighbours[7-k], k ) )        

        self.requests={}

    def init_fill_optimized(self,x,vname):
        """ initialize the optimized send/recv messages = bypassing
        the intermediate buffers """
        n = self.nh
        comm = self.comm
        neighbours = self.neighbours
        reqs = []
        reqs = []
        if vname=='x':
            reqs.append( comm.Ssend_init(x[n:n+n,n:n+n],neighbours[0], 0 ) )
            reqs.append( comm.Ssend_init(x[n:n+n,n:-n],neighbours[1], 1 ) )
            reqs.append( comm.Ssend_init(x[n:n+n,-n-n:-n],neighbours[2], 2 ) )
            reqs.append( comm.Ssend_init(x[n:-n,n:n+n],neighbours[3], 3 ) )
            reqs.append( comm.Ssend_init(x[n:-n,-n-n:-n],neighbours[4], 4 ) )
            reqs.append( comm.Ssend_init(x[-n-n:-n,n:n+n],neighbours[5], 5 ) )
            reqs.append( comm.Ssend_init(x[-n-n:-n,n:-n],neighbours[6], 6 ) )
            reqs.append( comm.Ssend_init(x[-n-n:-n,-n-n:-n],neighbours[7], 7 ) )

            reqr.append( comm.Recv_init(x[-n:,-n:] ,neighbours[7-0], 0 ) ) 
            reqr.append( comm.Recv_init(x[-n:,n:-n] ,neighbours[7-1], 1 ) ) 
            reqr.append( comm.Recv_init(x[-n:,:n] ,neighbours[7-2], 2 ) ) 
            reqr.append( comm.Recv_init(x[n:-n,-n:] ,neighbours[7-3], 3 ) ) 
            reqr.append( comm.Recv_init(x[n:-n,:n] ,neighbours[7-4], 4 ) ) 
            reqr.append( comm.Recv_init(x[:n,-n:] ,neighbours[7-5], 5 ) ) 
            reqr.append( comm.Recv_init(x[:n,n:-n] ,neighbours[7-6], 6 ) ) 
            reqr.append( comm.Recv_init(x[:n,:n] ,neighbours[7-7], 7 ) ) 
        

            self.requests['xr']=reqr
            self.requests['xs']=reqs


    def fill_optimized(self,x,vname):
        """ an optimized halo filling that directly maps variables
        without relying on a third party buffer """
        
        ok = True

        if vname=='x':
            reqr = self.requests['xr']
            reqs = self.requests['xs']
        # elif vname=='r':
        #     reqr = self.requests['rr']
        #     reqs = self.requests['rs']
        # elif vname=='b':
        #     reqr = self.requests['br']
        #     reqs = self.requests['bs']
        else:
            print('unrecognized vname in fill_optimized')
            ok = False

        if ok:
            MPI.Prequest.Startall(reqr)            
            MPI.Prequest.Startall(reqs) 
            MPI.Prequest.Waitall(reqr)        
            MPI.Prequest.Waitall(reqs)
        else:
            self.fill(x)

            
    def fill(self,x):
        """ fill the halo of x"""
        # fill the halo if proc only knows a subdomain
        if self.method==0:
            for k in range(self.nbhalo):
                fillhalo(x,self.nh)
#            x=zeros((self.m,self.n))
        elif self.method==1:

            n = self.nh

            reqs = self.reqs
            reqr = self.reqr

#            x=zeros((self.m,self.n))
#            MPI.Prequest.Startall(reqr)            

            # dump target into preallocated buff0
            # the [:] is mandatory to stream the data 
            # in that particular memory place
            self.sbuff[0][:,:]=x[   n:n+n,   n:n+n]
            MPI.Prequest.Start(reqs[0])

            self.sbuff[1][:,:]=x[   n:n+n,   n:-n] 
            MPI.Prequest.Start(reqs[1])

            self.sbuff[2][:,:]=x[   n:n+n,   -n-n:-n] 
            MPI.Prequest.Start(reqs[2])
                                                                           
            self.sbuff[3][:,:]=x[    n:-n,  n:n+n]
            MPI.Prequest.Start(reqs[3])

            self.sbuff[4][:,:]=x[     n:-n,-n-n:-n] 
            MPI.Prequest.Start(reqs[4])

            self.sbuff[5][:,:]=x[    -n-n:-n,n:n+n]
            MPI.Prequest.Start(reqs[5])

            self.sbuff[6][:,:]=x[    -n-n:-n,n:-n] 
            MPI.Prequest.Start(reqs[6])

            self.sbuff[7][:,:]=x[ -n-n:-n,-n-n:-n]  
            MPI.Prequest.Start(reqs[7])

            # if self.myrank==1:
            #     for k in range(8):
            #         print(k,self.sbuff[k])

            MPI.Prequest.Wait(reqr[0])
            x[-n-n:-n, -n-n:-n] = self.rbuff[0]

            MPI.Prequest.Wait(reqr[1])
            x[    -n-n:-n,n:-n]=self.rbuff[1]

            MPI.Prequest.Wait(reqr[2])
            x[    -n-n:-n,n:n+n]=self.rbuff[2]

            MPI.Prequest.Wait(reqr[3])
            x[   n:-n,-n-n:-n] =self.rbuff[3]
                                 
            MPI.Prequest.Wait(reqr[4])
            x[       n:-n,n:n+n]=self.rbuff[4]

            MPI.Prequest.Wait(reqr[5])
            x[   :n,-n-n:-n]=self.rbuff[5]
                                 
            MPI.Prequest.Wait(reqr[6])
            x[     :n,n:-n]=self.rbuff[6]

            MPI.Prequest.Wait(reqr[7])
            x[      :n,:n]=self.rbuff[7]
            
            # if self.myrank==1:
            #     for k in range(8):
            #         print(k,self.rbuff[k])

        # by-pass MPI exchange if level is coarse enough so that 
        # the whole domain fits in
        elif self.method==2:

            reqs = self.reqs
            reqr = self.reqr
            
            for k in range(self.nbhalo):

                MPI.Prequest.Startall(reqr)            
#                for l in range(8):
#                    MPI.Prequest.Start(reqr[l])


    #            Print(x.shape,b0.shape,b1.shape,b2.shape,b3.shape,b4.shape,b5.shape,b6.shape,b7.shape)
    #            print(b0.flags,x.flags)

                # n  = self.n
                # m  = self.m
                # nh = self.nh
                # n2 = n-2*nh
                # m2 = m-2*nh

                # b0=self.sbuff[0]
                # b1=self.sbuff[1]
                # b2=self.sbuff[2]
                # b3=self.sbuff[3]
                # b4=self.sbuff[4]
                # b5=self.sbuff[5]
                # b6=self.sbuff[6]
                # b7=self.sbuff[7]

                halotobuffer(x,
                             self.sbuff[0],self.sbuff[1],self.sbuff[2],
                             self.sbuff[3],              self.sbuff[4],
                             self.sbuff[5],self.sbuff[6],self.sbuff[7])


                # for k in range(8):
                #     self.sbuff[k][:,:]=result[k]        

                # if self.myrank==1:
                #     for k in range(8):
                #         print(k,self.sbuff[k])

                # if self.myrank==1:
                #     print(b0,b1,b2,b3,b4,b5,b6,b7)

                # for l in [0,2,5,7,1,3,4,6]:
                #     MPI.Prequest.Start(reqs[l])
                #     MPI.Prequest.Wait(reqr[l])


                MPI.Prequest.Startall(reqs) 
                MPI.Prequest.Waitall(reqr)
#                MPI.COMM_WORLD.Barrier()


                # b0=self.rbuff[0]
                # b1=self.rbuff[1]
                # b2=self.rbuff[2]
                # b3=self.rbuff[3]
                # b4=self.rbuff[4]
                # b5=self.rbuff[5]
                # b6=self.rbuff[6]
                # b7=self.rbuff[7]

                # if self.myrank==1:
                #     for k in range(8):
                #         print(k,self.rbuff[k])
                # if self.myrank==1:
                #     print(b0,b1,b2,b3,b4,b5,b6,b7)

    #            x=buffertohalo(x,b0,b1,b2,b3,b4,b5,b6,b7)

                buffertohalo(x,
                             self.rbuff[0],self.rbuff[1],self.rbuff[2],
                             self.rbuff[3],self.rbuff[4],
                             self.rbuff[5],self.rbuff[6],self.rbuff[7])

                MPI.Prequest.Waitall(reqs)
#                MPI.COMM_WORLD.Barrier()

        elif self.method==3:
            comm = MPI.COMM_WORLD
            for l in range(self.nbhalo):
                halotobuffer(x,
                             self.sbuff[0],self.sbuff[1],self.sbuff[2],
                             self.sbuff[3],              self.sbuff[4],
                             self.sbuff[5],self.sbuff[6],self.sbuff[7])
            
                req=[]
                for k in range(8):
                    # exchange data
                    req.append(comm.issend(self.sbuff[k],dest=self.neighbours[k],tag=k))

                for k in range(8):
                    self.rbuff[k]=comm.recv(source=self.neighbours[7-k],tag=k)

                for k in range(8):
                    req[k].Wait(MPI.Status())

                buffertohalo(x,
                             self.rbuff[0],self.rbuff[1],self.rbuff[2],
                             self.rbuff[3],self.rbuff[4],
                             self.rbuff[5],self.rbuff[6],self.rbuff[7])



        elif self.method==5:
            comm = MPI.COMM_WORLD
            myrank=self.myrank
            status = MPI.Status()
            neighb = self.neighbours
            n = self.nh
            reqs = self.reqs
            reqr = self.reqr

            #y = zeros((self.m,self.n))
            #y[:,:] = x[:,:].copy()
            
            reqr = []
            for k in range(8):
                reqr.append(comm.Irecv(self.rbuff[k],neighb[7-k], k ))

            #MPI.Prequest.Startall(reqr) 

            self.sbuff[7][:,:]=x[ -n-n:-n,-n-n:-n]  
            comm.Isend(self.sbuff[7],neighb[7],7)

            self.sbuff[6][:,:]=x[    -n-n:-n,n:-n] 
            comm.Isend(self.sbuff[6],neighb[6],6)

            self.sbuff[5][:,:]=x[    -n-n:-n,n:n+n]
            comm.Isend(self.sbuff[5],neighb[5],5)

            self.sbuff[4][:,:]=x[     n:-n,-n-n:-n] 
            comm.Isend(self.sbuff[4],neighb[4],4)

            self.sbuff[3][:,:]=x[    n:-n,  n:n+n]
            comm.Isend(self.sbuff[3],neighb[3],3)

            self.sbuff[2][:,:]=x[   n:n+n,   -n-n:-n] 
            comm.Isend(self.sbuff[2],neighb[2],2)

            self.sbuff[1][:,:]=x[   n:n+n,   n:-n] 
            comm.Isend(self.sbuff[1],neighb[1],1)

            self.sbuff[0][:,:]=x[   n:n+n,   n:n+n]
            comm.Isend(self.sbuff[0],neighb[0],0)

#            MPI.Prequest.Start(reqs[0])

#            MPI.Prequest.Start(reqs[1])

#            MPI.Prequest.Start(reqs[2])
                                                                           
#            MPI.Prequest.Start(reqs[3])

#            MPI.Prequest.Start(reqs[4])

#            MPI.Prequest.Start(reqs[5])

#            MPI.Prequest.Start(reqs[6])

#            MPI.Prequest.Start(reqs[7])


            # halotobuffer(x,
            #              self.sbuff[0],self.sbuff[1],self.sbuff[2],
            #              self.sbuff[3],              self.sbuff[4],
            #              self.sbuff[5],self.sbuff[6],self.sbuff[7])
            # self.sbuff[0][:,:]=x[   n:n+n,   n:n+n]
            # self.sbuff[1][:,:]=x[   n:n+n,   n:-n] 
            # self.sbuff[2][:,:]=x[   n:n+n,   -n-n:-n] 
            # self.sbuff[3][:,:]=x[    n:-n,  n:n+n]
            # self.sbuff[4][:,:]=x[     n:-n,-n-n:-n] 
            # self.sbuff[5][:,:]=x[    -n-n:-n,n:n+n]
            # self.sbuff[6][:,:]=x[    -n-n:-n,n:-n] 
            # self.sbuff[7][:,:]=x[ -n-n:-n,-n-n:-n]  

#            MPI.Prequest.Startall(reqs) 
#            MPI.Prequest.Waitall(reqr)

            ok = False
            while not(ok):
#            for j in range(8):
                out = MPI.Prequest.Waitany(reqr,status=status)        
                ok = (out<0)
                j = status.tag
                # k = status.source
                # if self.myrank==0:
                #     print('tag,source,out=>',j,k,out)
                #     print('status.source=',status.source,'tag=',status.tag,'error=',status.error)

                if j==7:
                    #if myrank==0:
                    #    print(x[  :n ,  :n],self.rbuff[7])
                    x[  :n ,  :n] = self.rbuff[7][:,:]
                    #if myrank==0:
                    #    print('*',x[  :n ,  :n])
                if j==6:
                    #print(j,self.rbuff[6][:,:])
                    x[  :n , n:-n]= self.rbuff[6][:,:]
                if j==5:
                     x[  :n ,-n:]= self.rbuff[5][:,:]
                if j==4:
                     x[ n:-n,  :n]= self.rbuff[4][:,:]
                if j==3:
                     x[ n:-n,-n:]= self.rbuff[3][:,:]
                if j==2:
                     x[-n:  ,  :n]= self.rbuff[2][:,:]
                if j==1:
                     x[-n:  , n:-n]= self.rbuff[1][:,:]
                if j==0:
                     x[-n:  ,-n:]= self.rbuff[0][:,:]
#                if (self.myrank==0) and not(ok):
#                    print(j)
                #     #print(j,self.rbuff[j])
                #     print(j,self.n,self.m)
                #     #print(x)

            # numpy.set_printoptions(threshold=numpy.nan)

            # buffertohalo(y,
            #              self.rbuff[0],self.rbuff[1],self.rbuff[2],
            #              self.rbuff[3],self.rbuff[4],
            #              self.rbuff[5],self.rbuff[6],self.rbuff[7])

            # diff = where( x!=y )  

            # if len( diff[0] )>0 and self.myrank==0:
            #     print('***diff',len(diff),self.myrank,self.n)
            #     print(diff)
            #     msk=zeros((self.m,self.n),int)
            #     msk[diff]=1
            #     print(msk)
            #     #print(x[diff])
            #     #print(y[diff])
            # elif self.myrank==0:
            #     print('ok!')
            # x = y.copy()
            #MPI.Prequest.Waitall(reqs)


        elif self.method==4:
            # don't use preallocated buffers
            # send directly blocks of x and 
            # receive directly msg into x
            n = self.nh

            req = []
            
            sbuff[0]=x[   n:n+n,   n:n+n]
            rbuff[0]=x[:]
            req[0]=comm.issend(sbuff[0],dest=self.neighbours[0],tag=0)
            rbuff[0][:,:]=comm.recv(source=self.neighbours[7],tag=0)
            

            self.sbuff[1][:,:]=x[   n:n+n,   n:-n] 
            MPI.Prequest.Start(reqs[1])

            self.sbuff[2][:,:]=x[   n:n+n,   -n-n:-n] 
            MPI.Prequest.Start(reqs[2])
                                                                           
            self.sbuff[3][:,:]=x[    n:-n,  n:n+n]
            MPI.Prequest.Start(reqs[3])

            self.sbuff[4][:,:]=x[     n:-n,-n-n:-n] 
            MPI.Prequest.Start(reqs[4])

            self.sbuff[5][:,:]=x[    -n-n:-n,n:n+n]
            MPI.Prequest.Start(reqs[5])

            self.sbuff[6][:,:]=x[    -n-n:-n,n:-n] 
            MPI.Prequest.Start(reqs[6])

            self.sbuff[7][:,:]=x[ -n-n:-n,-n-n:-n]  
            MPI.Prequest.Start(reqs[7])

            # if self.myrank==1:
            #     for k in range(8):
            #         print(k,self.sbuff[k])

            MPI.Prequest.Wait(reqr[0])
            x[  :n ,  :n]  = self.rbuff[0]

            MPI.Prequest.Wait(reqr[1])
            x[  :n , n:-n]=self.rbuff[1]

            MPI.Prequest.Wait(reqr[2])
            x[  :n ,-n:]=self.rbuff[2]

            MPI.Prequest.Wait(reqr[3])
            x[ n:-n,  :n] =self.rbuff[3]
                                 
            MPI.Prequest.Wait(reqr[4])
            x[ n:-n,-n:]=self.rbuff[4]

            MPI.Prequest.Wait(reqr[5])
            x[-n:  ,  :n]=self.rbuff[5]
                                 
            MPI.Prequest.Wait(reqr[6])
            x[-n:  , n:-n]=self.rbuff[6]

            MPI.Prequest.Wait(reqr[7])
            x[-n:  ,-n:]=self.rbuff[7]

        elif self.method==6:
            comm=MPI.COMM_WORLD.py2f()
            m = self.m
            n = self.n
            nh=self.nh
            neighb=self.neighbours
#            y=numpy.array(x,order='F')
            fillhalo_fortran(x,comm,neighb,nh)
#            if self.myrank==1:
#                print(x)

        else:
            n = self.nh

            x[  :n ,  :n] =0.
            x[  :n , n:-n]=0.
            x[  :n ,-n:]  =0.
                                
            x[ n:-n,  :n] =0.
            x[ n:-n,-n:]  =0.
                                
            x[-n:  ,  :n] =0.
            x[-n:  , n:-n]=0.
            x[-n:  ,-n:]  =0.
            
#        return x

    def timeit(self):
        nt = 1000
        x=zeros( (self.m,self.n),order='C')+self.myrank
        t0 = time()
        for kt in range(nt):
            self.fill(x)
        t1 = time()
        return (t1-t0)/nt
        

if  __name__ =='__main__':
    comm=MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nbofproc = MPI.COMM_WORLD.Get_size()

    nh=3
    n=256
    m=256

    ix=1
    iy=1
    np = int(sqrt(nbofproc))
    mp = nbofproc/np


    neighbours = set_neighbours(myrank,np,mp,ix,iy)
    if myrank==0:
        print('----------------------------------------')
        print('neighbours:')
    comm.Barrier()    
    print(myrank,neighbours)
    comm.Barrier()
    halo = Halo(nh,n,m,comm,neighbours,method=5)
    
    if myrank==0:
        print( halo.sbuff[0].flags)

    if myrank==0:
        print('----------------------------------------')
        print('Time:')

    for k in [2,6]:#range(-1,3):
        halo.method=k
        dt=halo.timeit()
        if myrank==0:
            print('method %i : %e'%(k,dt))
    comm.Barrier()

    halo.method=6
    x=zeros( (m,n))+myrank*1.
    #x,y=meshgrid(arange(n)*1.,arange(m)*1.)

    x=zeros( (m,n))+myrank*1.
    #x=numpy.random.randn(m,n)
    for k in range(10):
        halo.fill(x)

    comm.Barrier()

    for k in [0]:#range(nbofproc):
        if myrank==k:
            print('----------------------------------------')
            print('x: ',x.shape)
            print(x)
        comm.Barrier()
