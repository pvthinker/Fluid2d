!----------------------------------------
      subroutine computeorthogradient(msk,psi,dx,dy,nh,n,m,u,v)

! compute u,v= (-dpsi/dy,dpsi/dx)
!
 
      implicit none

      integer,intent(in):: n,m,nh
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n) :: psi,u,v
      real*8::dx,dy

!f2py intent(inplace)::psi,msk,u,v

      real*8::zdx,zdy
      integer:: i,j
      integer*1:: mm

      zdx = 1./dx
      zdy = 1./dy
      do j=2,m-1!nh-1,m-nh+1
         do i=2,n-1!nh-1,n-nh+1
            mm=msk(j,i)+msk(j,i+1)
            if(mm.eq.2)then
               u(j,i) = zdy*(psi(j-1,i)-psi(j,i))
            else
               u(j,i) = 0.
            endif
            mm=msk(j,i)+msk(j+1,i)
            if(mm.eq.2)then
               v(j,i) = zdx*(psi(j,i)-psi(j,i-1))
            else
               v(j,i) = 0.
            endif
         enddo
      enddo

      end subroutine



!----------------------------------------
      subroutine celltocorner(xr,xp,n,m)

! average 4 cell centers at upper right cell corner
!
 
      implicit none

      integer,intent(in):: n,m
      real*8,dimension(m,n) :: xr,xp

!f2py intent(inplace)::xr,xp

      integer:: i,j

      do j=1,m-1
         do i=1,n-1
            xp(j,i)=0.25*(xr(j,i)+xr(j,i+1)+xr(j+1,i)+xr(j+1,i+1))
         enddo
      enddo

      end subroutine

!----------------------------------------
      subroutine celltocornerbicubic(xr,xp,n,m)

! average 4x4 cell centers at upper right cell corner
! using the finite difference interpolation [-1,9,9,1]/16 in each direction
 
      implicit none

      integer,intent(in):: n,m
      real*8,dimension(m,n) :: xr,xp

!f2py intent(inplace)::xr,xp
      real:: cff
      integer:: i,j,jj,k1,k2,k3,k4
      real*8,dimension(n,4)::z

      cff = 1./12.
      do j=1,m
         k1 = mod(j-1,4)+1
         do i=2,n-2
            z(i,k1)=cff*(-xr(j,i-1)+7*(xr(j,i)+xr(j,i+1))-xr(j,i+2))
         enddo
         if (j.ge.4)then
            k2 = mod(j-2,4)+1
            k3 = mod(j-3,4)+1
            k4 = mod(j-4,4)+1
            jj = j-2
            do i=2,n-2
               xp(jj,i)=cff*(-z(i,k1)+7*(z(i,k2)+z(i,k3))-z(i,k4))
            enddo
         endif
      enddo

      end subroutine

!----------------------------------------
      subroutine cornertocell(xp,xr,n,m)

! average 4 cell centers at upper right cell corner
!
 
      implicit none

      integer,intent(in):: n,m
      real*8,dimension(m,n) :: xr,xp

!f2py intent(inplace)::xr,xp

      integer:: i,j

      do j=2,m
         do i=2,n
            xr(j,i)=0.25*(xp(j,i)+xp(j,i-1)+xp(j-1,i)+xp(j-1,i-1))
         enddo
      enddo

      end subroutine

!----------------------------------------
      subroutine add_diffusion(msk,trac,dx,nh,Kdiff,dtrac,n,m)

! add a diffusion term (with homogenous Neumann BC) on trac
!
 
      implicit none

      integer,intent(in):: n,m,nh
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n) :: trac,dtrac
      real*8,intent(in)::dx,Kdiff

!f2py intent(inplace)::msk,trac,dtrac

      real*8::coef,dbl,dbr
      integer:: i,j
      integer*1:: mm,ml,mr

      coef = Kdiff/(dx*dx)
      do j=2,m-1
         do i=2,n-1
            if (msk(j,i).eq.1) then
               dtrac(j,i) = dtrac(j,i) + coef*( &
                    + msk(j,i-1) * ( trac(j,i-1)-trac(j,i) ) &
                    + msk(j,i+1) * ( trac(j,i+1)-trac(j,i) ) &
                    + msk(j-1,i) * ( trac(j-1,i)-trac(j,i) ) &
                    + msk(j+1,i) * ( trac(j+1,i)-trac(j,i) ) )
            endif
         enddo
      enddo

      end subroutine


!----------------------------------------
      subroutine computewallshear(msk,x,y,total,dx,nh,m,n)

! compute y=b-A*x

      implicit none

      integer,intent(in):: n,m,nh
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n) :: x
      real*8,dimension(m,n) :: y
      real*8::dx,total

!f2py intent(inplace)::msk,x
!f2py intent(inplace)::y
!f2py intent(out)::total

!      real*8,dimension(n,3) :: xo
      integer:: i,j
      integer*1::msku,mskv
      real*8::cff,u,v

      cff=1./(2*dx*dx) 

!      y(:,:)=0.
      total=0.
      y(1:nh,:)=0.
      do j=nh+1,m-nh+1
         y(j,1:nh)=0.
         do i=nh+1,n-nh+1
            y(j,i)=0.
            msku=msk(j,i-1)+msk(j,i)
            if (msku.eq.1)then
               if(msk(j,i).ne.0)then
                  v=(x(j,i+1)+x(j-1,i+1))*cff
                  y(j,i)=y(j,i)-v
                  total=total-v
               else
                  v=(x(j,i-3)+x(j-1,i-3))*cff
                  y(j,i-1)=y(j,i-1)+v
                  total=total+v
               endif
            endif
            mskv=msk(j-1,i)+msk(j,i)
            if (mskv.eq.1)then
               if(msk(j,i).ne.0)then
                  u=(x(j+1,i)+x(j+1,i-1))*cff
                  y(j,i)=y(j,i)+u
                  total=total+u
               else
                  u=(x(j-3,i)+x(j-3,i-1))*cff
                  y(j-1,i)=y(j-1,i)-u
                  total=total-u
               endif
            endif
         enddo
      enddo


      end subroutine 

!----------------------------------------
      subroutine computenoslipsourceterm(msk,x,y,total,dx,dy,nh,m,n)

! compute y=b-A*x

      implicit none

      integer,intent(in):: n,m,nh
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n) :: x
      real*8,dimension(m,n) :: y
      real*8::dx,dy,total

!f2py intent(inplace)::msk,x
!f2py intent(inplace)::y
!f2py intent(out)::total

!      real*8,dimension(n,3) :: xo
      integer:: i,j
      integer*1::msku,mskv
      real*8::cff,u,v

      cff=1./(2*dx*dy)

!      y(:,:)=0.
      total=0.
      y(1:nh,:)=0.
      do j=nh+1,m-nh+1
         y(j,1:nh)=0.
         do i=nh+1,n-nh+1
            y(j,i)=0.
            msku=msk(j,i-1)+msk(j,i)
            if (msku.eq.1)then
               v=(x(j,i)+x(j-1,i)-x(j,i-2)-x(j-1,i-2))*cff
               if(msk(j,i).ne.0)then
                  y(j,i)=y(j,i)-v
                  total=total-v
               else
                  y(j,i-1)=y(j,i-1)+v
                  total=total+v
               endif
            endif
            mskv=msk(j-1,i)+msk(j,i)
            if (mskv.eq.1)then
               u=-(x(j,i)+x(j,i-1)-x(j-2,i)-x(j-2,i-1))*cff
               if(msk(j,i).ne.0)then
                  y(j,i)=y(j,i)+u
                  total=total+u
               else
                  y(j-1,i)=y(j-1,i)-u
                  total=total-u
               endif
            endif
         enddo
      enddo


      end subroutine 


!----------------------------------------
      subroutine computenoslipsourceterm_regularized(msk,x,y,cu,cv,total,dx,nh,m,n)

! compute y=b-A*x

      implicit none

      integer,intent(in):: n,m,nh
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n) :: x,cu,cv
      real*8,dimension(m,n) :: y
      real*8::dx,total

!f2py intent(inplace)::msk,x,cu,cv
!f2py intent(inplace)::y
!f2py intent(out)::total

!      real*8,dimension(n,3) :: xo
      integer:: i,j
      integer*1::msku,mskv
      real*8::cff,ul,ur,vl,vr,z

      cff=2./(dx*dx)

!      y(:,:)=0.
      total=0.
      y(1:nh,:)=0.
      do j=nh+1,m-nh+1
         y(j,1:nh)=0.
         do i=nh+1,n-nh+1
            if(msk(j,i).ne.0)then               
               vl = (x(j,i)-x(j,i-1))*cv(j,i-1)/(msk(j,i)+msk(j,i-1))
               vr = (x(j,i)-x(j,i+1))*cv(j,i)/(msk(j,i)+msk(j,i+1))

               ul = (x(j,i)-x(j-1,i))*cu(j-1,i)/(msk(j,i)+msk(j-1,i))
               ur = (x(j,i)-x(j+1,i))*cu(j,i)/(msk(j,i)+msk(j+1,i))
               
               z = -(vl+ur+vr+ul)*cff
               y(j,i) = z
               total=total+z
            else
               y(j,i)=0.
            endif
         enddo
      enddo


      end subroutine 

!----------------------------------------
      subroutine add_torque(msk,buoy,dx,nh,gravity,domega,n,m)

! add torque db/dx to domega
!
 
      implicit none

      integer,intent(in):: n,m,nh
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n) :: buoy,domega
      real*8,intent(in)::dx,gravity

!f2py intent(inplace)::msk,buoy,domega

      real*8::coef,dbl,dbr,sum,z
      integer:: i,j,ii
      integer*1:: ml,mr,mm,one

      coef = 0.5*gravity/(dx)

      one = 1
      do j=1+nh,m-nh
!         sum=0.
         i=1+nh
         ml  =  max(one,msk(j,i)+ msk(j,i-1))
         dbl = (buoy(j,i)*msk(j,i)+buoy(j,i-1)*msk(j,i-1))/ml
         do i=1+nh,n-nh

            mr = max(one,msk(j,i+1)+ msk(j,i))
            dbr=(buoy(j,i+1)*msk(j,i+1)+buoy(j,i)*msk(j,i))/mr

            mm=(ml+mr)
!            if(mm.gt.0) then
!               z = (dbr-dbl)*(coef)
!            else
!               z=0.
!            endif
            if (mm.eq.4) then
!                 domega(j,i) = domega(j,i) + (dbr-dbl)*(coef)*msk(j,i)
                 domega(j,i) = domega(j,i) + (buoy(j,i+1)-buoy(j,i-1))*(coef)*msk(j,i)
            endif
!            sum=sum+z
!            if(mm.gt.0) domega(j,i) =  domega(j,i)+(dbr-dbl)*(coef/mm)
            ml=mr
            dbl=dbr
!            ii=ii+1
         enddo
!         write(*,*)"sum=",j,ii,sum
      enddo


      end subroutine 


