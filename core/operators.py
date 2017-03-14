from param import Param
from numpy import zeros, shape,zeros_like,sum,roll,exp,where
from gmg.hierarchy import Gmg
from fortran_advection import *
from fortran_operators import *
from fortran_diag import *


class Operators(Param):

    def __init__(self,param,grid):
        self.list_param=['varname_list','tracer_list','whosetspsi','mpi','npx','npy',
                         'nh','gravity','f0','beta','Rd','qgoperator','order','Kdiff',
                         'enforce_momentum','isisland']
        param.copy(self,self.list_param)

        self.list_grid=['msk','nxl','nyl','dx','dy','bcarea','mpitools','msknoslip',
                        'mskbc','domain_integration','nh','xr0','yr0']

        grid.copy(self,self.list_grid)


        # internal work array for the inversion
        self.work = zeros((self.nyl,self.nxl))
        self.work2 = zeros((self.nyl,self.nxl))

        pp = {'np':param.npx,'mp':param.npy,'nh':param.nh,
              'n':param.nx//param.npx,'m':param.ny//param.npy,
              'omega':8./9.,'npmpmax':1,'verbose':False,
              'dx':grid.dx,'n1':32,'n0':8,'method':'deep','nagglo':2}


        # load the multigrid solver
        #
        #WARNING: the multigrid needs the mask at cell corners!!!
        #         not at cell centers
        mskr = self.msk*1.

        # this piece is a bit awkward: to initialize gmg, we need
        # a mask with a halo properly filled but the fill_halo method
        # belongs to gmg. We have a circular definition.
        # the trick: define a dummy gmg first a msk=1 everywhere
        # then grab the fill_halo method and redefine once again the
        # multigrid, this time with the proper mask
        #self.gmg = Gmg(pp,mskr) 
        # borrow the fill_halo from the multigrid
        #self.fill_halo = self.gmg.grid[0].halo.fill

        celltocorner(mskr,self.work)
        #self.fill_halo(self.work)

        #del self.gmg
        #del self.fill_halo

        self.work[self.work<1.]=0.        
        self.mskp=self.msk*0
        self.mskp[self.work==1.]=1
        pp['verbose']=True
        if hasattr(self,'qgoperator'):
            pp['qgoperator']=True
            pp['Rd']=self.Rd
            self.gmg = Gmg(pp,self.work)
        else:
            self.gmg = Gmg(pp,self.work) 

        # borrow the fill_halo from the multigrid
        self.fill_halo = self.gmg.grid[0].halo.fill
        grid.fill_halo = self.gmg.grid[0].halo.fill
        
        self.blwidth = param.Lx*0.05

        self.set_boundary_msk()


        # select the proper flux discretization
        if self.order%2==0:
            self.fortran_adv = adv_centered
        else:
            self.fortran_adv = adv_upwind


        self.cst=zeros(3,)
        self.cst[0]=grid.dx
        self.cst[1]=0.05
        # should be updated at each timestep
        #self.cst[2]=param.umax

        # controls the flux splitting method
        # 0 = min/max
        # 1 = parabolic
        self.method=1

        # these coefficients below are used for the thermalwind model
        coef = zeros_like(self.msk)
        coef[1:-1,1:-1] = self.msk[1:-1,2:]+self.msk[1:-1,0:-2]
        coef[coef<2]=0.
        coef[coef==2]=0.5
        self.coefb=coef*1.
        
        coef = zeros_like(self.msk)
        coef[1:-1,1:-1] = self.msk[2:,1:-1]+self.msk[0:-2,1:-1]
        coef[coef<2]=0.
        coef[coef==2]=0.5
        self.coefV=coef*1.

        if type(self.Kdiff) != dict:
            K = self.Kdiff
            self.Kdiff={}
            for trac in self.tracer_list:
                self.Kdiff[trac]=K


    def set_boundary_msk(self):
        """ for the no slip boundary source term """
        nh = self.nh
        msk=self.msknoslip
        z = roll(msk,-1,axis=1)+roll(msk,-1,axis=0)+roll(msk,+1,axis=1)+roll(msk,+1,axis=0)-4*msk
        z = z*msk
        self.mskbc = self.msk*0
        self.mskbc[z<0]=1            
        # the halo will be filled later in operator.py when fill_halo will become available
        # we can now fix the boundary mask
        self.mskbc *= self.msknoslip
        self.fill_halo(self.mskbc)

        # idx in the 2D array where boundary terms are computed
        # used for storage and i/o
        #self.idxbc = where(self.mskbc==1)

        self.bcarea = self.domain_integration(self.mskbc)
        self.x2bc = self.domain_integration((self.xr0)**2*self.mskbc*self.msknoslip) 
        self.y2bc = self.domain_integration((self.yr0)**2*self.mskbc*self.msknoslip) 

        return

        def smooth(msk,msk0,k):
            y = roll(msk,-1,axis=1)+roll(msk,-1,axis=0)+roll(msk,+1,axis=1)+roll(msk,+1,axis=0)
            z = msk*1.
            z[(y>0) ]=k+1
            z[msk>0]=msk[msk>0]
            z[msk0==0]=0
            self.fill_halo(z)
            return z


        z0=1-msk*1.        
        nk = int(round( (self.blwidth/self.dx)))
        nk=1
        for k in range(nk):
            z = z0*1.
            z0 = smooth(z,msk,k)
        z0[z0==0]=nk
        z0=z0/nk
        z0*=msk
        z0=(1.-z0)**(nk/2.)
        z0[msk==0]=1
        self.cv = (roll(z0,-1,axis=1)+z0)*.5
        self.cu = (roll(z0,-1,axis=0)+z0)*.5
#        self.mskbc = z0
#        self.bcarea = self.domain_integration(z0)


    def rhs_adv(self,x,t,dxdt):
        """ compute -div(u*tracer) using finite volume flux discretization 
        the flux is computed at edge cells using p-th order interpolation
        for p even, the flux is centered
        for p odd, the flux is upwinded (more points on the upwind side) """
        iu = self.varname_list.index('u')
        iv = self.varname_list.index('v')
        u   = x[iu]
        v   = x[iv]

        for trac in self.tracer_list:
            ik = self.varname_list.index(trac)
            y = dxdt[ik]
            self.fortran_adv(self.msk,x[ik],y,u,v,
                       self.cst,self.nh,self.method,self.order)
            self.fill_halo(y)
            # for an unknown reason dxdt[ik] is
            # not updated by the Fortran routine
            # it should be done manually
            # (this yields an excessive data movement)
            dxdt[ik] = y
            #self.fill_halo(dxdt[ik])

    def wallshear(self,x,shear):
        ip = self.varname_list.index('psi')

        meansource = computewallshear(self.msk,x[ip],shear,self.dx,self.nh)
        
    def rhs_noslip(self,x,source):
        """ add the vorticity source term along the boundary to enforce
        zero tangential velocity (=no-slip) """

        ip = self.varname_list.index('psi')
        iu = self.varname_list.index('u')
        iv = self.varname_list.index('v')
        iw = self.varname_list.index(self.whosetspsi)

        cornertocell(x[ip],self.work)

#        y = source*0
 #       meansource = 0.
        meansource = computenoslipsourceterm(self.msknoslip,x[ip],self.work,self.dx,self.nh)

#        meansource = computenoslipsourceterm_regularized(self.msk,self.work,y,self.cu,self.cv,self.dx,self.nh)
        source[:,:] = self.work



        # this step is SUPER important to ensure GLOBAL vorticity conservation
        meansource=self.domain_integration(source) / self.bcarea

        source -= meansource*self.mskbc

        if self.enforce_momentum:
            xr = self.xr0
            yr = self.yr0
            # this step ensures the zero momentum
            px = computedotprod(self.msk,source,xr,self.nh)
            py = computedotprod(self.msk,source,yr,self.nh)
            cst = self.mpitools.local_to_global( [(px,'sum'),(py,'sum')])

            px,py = cst[0]/self.x2bc, cst[1]/self.y2bc
            source -= (px*xr+py*yr)*self.mskbc

        self.fill_halo(source)
        x[iw] -= source 


    def rhs_diffusion(self,x,t,dxdt):
        """ add a diffusion term on the tracer variables """

        for trac in self.tracer_list:
            ik = self.varname_list.index(trac)
            y = dxdt[ik]
            add_diffusion(self.msk,x[ik],self.dx,self.nh,self.Kdiff[trac],y)
            self.fill_halo(y)
            dxdt[ik] = y

    def rhs_torque(self,x,t,dxdt):
        """ compute g*db/dx for the Boussinesq model """
        ib = self.varname_list.index('buoyancy')
        iw = self.varname_list.index('vorticity')
        
        y = dxdt[iw]

        # nh=self.nh
        # z = y*self.msk
        # z=z[nh:-nh,nh:-nh]
        # s1 = sum( z.ravel())
        b = x[ib]
        #self.fill_halo(b)

        add_torque(self.msk,b,self.dx,self.nh,self.gravity,y)
        self.fill_halo(y)
        dxdt[iw] = y
        # z = y*self.msk
        # z=z[nh:-nh,nh:-nh]
        # s2 = sum(z.ravel() )
        # print('ds=',s2-s1,s2,s1)

    def rhs_thermalwind(self,x,t,dxdt):

        iu = self.varname_list.index('u')
        ib = self.varname_list.index('buoyancy')
        iw = self.varname_list.index('vorticity')
        iV = self.varname_list.index('V')
        
        # add the themal wind balance
        # g*db/dx + f0*dV/dz
        # to domega/dt
        b = x[ib]
        V = x[iV]
        dw= dxdt[iw]

        dw[1:-1,1:-1] += self.coefb[1:-1,1:-1]*( b[1:-1,2:]-b[1:-1,0:-2] )*(self.gravity/self.dx)
#        dw[1:-1,1:-1] += self.coefV[1:-1,1:-1]*( V[2:,1:-1]-V[0:-2,1:-1] )*(self.f0/self.dx)
        dw[1:-1,1:-1] += ( V[2:,1:-1]-V[0:-2,1:-1] )*(self.f0/self.dx)
        #y = dxdt[iw]
        #add_bxfvz(self.msk,x[ib],x[iV],self.dx,self.nh,self.gravity,self.f0,y)
        #dxdt[iw] = y

        # add f0*u to dV/dt
        u=x[iu]
        dxdt[iV][:,1:] += 0.5*self.f0*( u[:,:-1]+u[:,1:])
        #y = dxdt[iV]
        #add_corio(self.msk,x[iu],self.dx,self.nh,self.f0,y)
        #dxdt[iV] = y


    def invert_vorticity(self,x,flag='full',island=False):
        """ compute psi from vorticity (or 'whosetspsi' in general)

        this routine interpolates the vorticity from cell centers to
        cell corners (where psi is defined)

        it then solves div*grad psi = omega with psi=0 along the boundary
        (Dirichlet condition) using a multigrid

        the non-divergent velocity is computed from psi"""
        iu = self.varname_list.index('u')
        iv = self.varname_list.index('v')
        ip = self.varname_list.index('psi')
        iw = self.varname_list.index(self.whosetspsi)

        u   = x[iu]
        v   = x[iv]
        psi = x[ip]
        
        celltocorner(x[iw],self.work)
        if island:
#            print('adding rhsp')
            self.work[:,:] -= self.rhsp            
            

        if flag=='fast':
            ite,res=self.gmg.twoVcycle(psi,self.work,
                                       {'maxite':1,'tol':1e-6,'verbose':True})    
        else:
            # compute to machine accuracy
            ite,res=self.gmg.solve(psi,self.work,
                                   {'maxite':8,'tol':1e-11,'verbose':True})

        # don't apply the fill_halo on it [because fill_halo as it is is applying periodic BC]
        psi = psi*self.mskp         
        if island:
#            print('adding psi')
            psi += self.psi

        # compute (u,v) @ U,V points from psi @ cell corner
        computeorthogradient(self.msk,psi,self.dx,self.nh,u,v)
        #self.fill_halo(u)
        #self.fill_halo(v)
        x[iu]=u
        x[iv]=v
        x[ip]=psi       


    def invert_twolayers(self,x,flag='full'):
        """ compute psi from vorticity (or 'whosetspsi' in general)

        this routine interpolates the vorticity from cell centers to
        cell corners (where psi is defined)

        it then solves div*grad psi = omega with psi=0 along the boundary
        (Dirichlet condition) using a multigrid

        the non-divergent velocity is computed from psi"""
        iu = self.varname_list.index('u')
        iv = self.varname_list.index('v')
        ip = self.varname_list.index('psi')
        iw = self.varname_list.index(self.whosetspsi)

        # psi_t = H1*psi_1 + H_2*psi_2
        # psi_c = 

        u   = x[iu]
        v   = x[iv]
        psi = x[ip]
        
        # interpolate vorticity from cell centers
        # to cell corners (where psi is defined)
        celltocorner(x[iw],self.work)

        if flag=='fast':
            ite,res=self.gmg.twoVcycle(psi,self.work,
                            {'maxite':1,'tol':1e-6,'verbose':False})    
        else:
            # compute to machine accuracy
            ite,res=self.gmg.solve(psi,self.work,
                            {'maxite':8,'tol':1e-11,'verbose':False})

        # don't apply the fill_halo on it [because fill_halo as it is is applying periodic BC]
        psi = psi*self.mskp         
        # compute (u,v) @ U,V points from psi @ cell corner
        computeorthogradient(self.msk,psi,self.dx,self.nh,u,v)
        #self.fill_halo(u)
        #self.fill_halo(v)
        x[iu]=u
        x[iv]=v
        x[ip]=psi       
