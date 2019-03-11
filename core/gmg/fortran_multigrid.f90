!----------------------------------------
      subroutine smoothtwicewithA(msk,A,x,b,coef,m,n,yo)
        
! compute 
!        z=x+(A*x-b)*coef
!        y=z+(A*z-b)*coef        
! using a buffer technique that
! minimizes memory access
!
! result is only valid in (i,j)=(3..n-2)x(3..n-2)
! need to have a halo width of 2 to be working

      implicit none

      integer,intent(in):: n,m
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n,5) :: A
!      real*8,dimension(m,n) :: D
      real*8,dimension(m,n) :: b
      real*8,dimension(m,n) :: x
      real*8,intent(in):: coef

!f2py intent(inplace)::msk,R,D,b,yo
!f2py intent(inplace)::x

      real*8,dimension(3,n) :: yo
!      real*8,dimension(n,2) :: bo
      integer:: i,j,k,jm,km,ko,kp,lm,lo,lp,jo,jmo,jmm,im,ip,jp
!      real*8::w
      integer*1::s
      
      real*8::c1,c2,c3

!      xo(:,1)=x(1,:)
!      xo(:,2)=x(2,:)
!      bo(:,1)=b(1,:)

!      km=1
!      ko=2
!      kp=3

      lm=1
      lo=2
      lp=3

      c1=coef
      c2=1.-c1

      do i=1,n
	yo(lo,i)=x(1,i)
      enddo	

      do j=2,m-1
         ! read x from the memory and store it into xo
!         xo(:,kp)=x(j+1,:)
!         bo(:,jo)=b(j,:)
         jm=j-1
         jp=j+1
         yo(lo,1)=x(j-1,1)
	 yo(lp,1)=x(j,1)
!	 yo(lp,n)=x(j,n)
         do i=2,n-1
            ip=i+1
            im=i-1
            ! compute the inner yo=A*x at row j
            if (msk(j,i).ne.0)then
               c3=c1/abs(A(j,i,5))
            yo(lp,i)= x(j,i)*c2 + c3*(
     $           +A(j,i,1)*x(j-1,i-1)
     $           +A(j,i,2)*x(j-1,i  )
     $           +A(j,i,3)*x(j-1,i+1)
     $           +A(j,i,4)*x(j,i-1)
!     $           +A(j,i,5)*x(j,i  )
     $           +A(j,ip,4)*x(j,ip)
     $           +A(jp,im,3)*x(jp,im)
     $           +A(jp,i,2)*x(jp,i  )
     $           +A(jp,ip,1)*x(jp,ip)
!!$     $           +A(j,i,6)*x(j,i+1)
!!$     $           +A(j,i,7)*x(j+1,i-1)
!!$     $           +A(j,i,8)*x(j+1,i  )
!!$     $           +A(j,i,9)*x(j+1,i+1) 
     $           -b(j,i) )
             else
                yo(lp,i)=0.
             endif
         enddo
!	 yo(lp,n)=x(j,n)
         jmm=jm-1
!	 x(jmm,2)=yo(lm,2)
!	 x(jmm,n-1)=yo(lm,n-1)
!	 x(jm,2)=yo(lo,2)
         if (jm.ge.3) then
             x(jmm,2)=yo(lm,2)
             x(jmm,n-1)=yo(lm,n-1)
             x(jm,2)=yo(lo,2)
            do i=3,n-2
               ip=i+1
               im=i-1
               ! yo is known at row j-2,j-1,j
               ! so we can compute y=A*yo at row jm=j-1
               if (msk(jm,i).ne.0)then

               c3=c1/abs(A(jm,i,5))
               x(jm,i)= yo(lo,i  )*c2 + c3*(
     $           +A(jm,i,1)*yo(lm,i-1)
     $           +A(jm,i,2)*yo(lm,i  )
     $           +A(jm,i,3)*yo(lm,i+1)
     $           +A(jm,i,4)*yo(lo,i-1)
     $           +A(jm,ip,4)*yo(lo,ip )
     $           +A(j,im,3) *yo(lp,im )
     $           +A(j,i,2)  *yo(lp,i   )
     $           +A(j,ip,1) *yo(lp,ip )
     $           -b(jm,i) )
               else
                  x(jm,i)=0.
               endif
              enddo
          endif
         ! swap the buffers by permuting the indices

         lm=mod(lm,3)+1
         lo=mod(lo,3)+1
         lp=mod(lp,3)+1

      enddo

      end subroutine

!----------------------------------------
      subroutine smoothoncewithA(msk,A,x,b,coef,m,n,yo)
        
! compute 
!        z=x+(A*x-b)*coef
!        y=z+(A*z-b)*coef        
! using a buffer technique that
! minimizes memory access
!
! result is only valid in (i,j)=(3..n-2)x(3..n-2)
! need to have a halo width of 2 to be working

      implicit none

      integer,intent(in):: n,m
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n,5) :: A
!      real*8,dimension(m,n) :: D
      real*8,dimension(m,n) :: b
      real*8,dimension(m,n) :: x
      real*8,intent(in):: coef

!f2py intent(inplace)::msk,R,D,b,yo
!f2py intent(inplace)::x

      real*8,dimension(3,n) :: yo
!      real*8,dimension(n,2) :: bo
      integer:: i,j,k,jm,km,ko,kp,lm,lo,lp,jo,jmo,jmm,im,ip,jp
!      real*8::w
      integer*1::s
      
      real*8::c1,c2,c3

!      xo(:,1)=x(1,:)
!      xo(:,2)=x(2,:)
!      bo(:,1)=b(1,:)

!      km=1
!      ko=2
!      kp=3

      lm=1
      lo=2
      lp=3

      c1=coef
      c2=1.-c1

!      do i=1,n
!	yo(lo,i)=x(1,i)
!      enddo	

      do j=2,m-1
         ! read x from the memory and store it into xo
!         xo(:,kp)=x(j+1,:)
!         bo(:,jo)=b(j,:)
         jm=j-1
         jp=j+1
!         yo(lo,1)=x(j-1,1)
!	 yo(lp,1)=x(j,1)
!	 yo(lp,n)=x(j,n)
         do i=2,n-1
            ip=i+1
            im=i-1
            ! compute the inner yo=A*x at row j
            if (msk(j,i).ne.0)then
               c3=c1/abs(A(j,i,5))
            x(j,i)= x(j,i)*c2 + c3*(
     $           +A(j,i,1)*x(jm,im)!x(j-1,i-1)
     $           +A(j,i,2)*x(jm,i)!x(j-1,i  )
     $           +A(j,i,3)*x(jm,ip)!x(j-1,i+1)
     $           +A(j,i,4)*x(j,im)!x(j,i-1)
     $           +A(j,ip,4)*x(j,ip)
     $           +A(jp,im,3)*x(jp,im)
     $           +A(jp,i,2)*x(jp,i  )
     $           +A(jp,ip,1)*x(jp,ip)
     $           -b(j,i) )
             else
                x(j,i)=0.
             endif
         enddo
      enddo

      end subroutine

!----------------------------------------
      subroutine smoothtridiag(msk,A,x,b,m,n)
        
! compute 
!        y=V^{-1}*(b-H*x)
! where      
! Laplacian = V+H
! V is the vertical part (including the diagonal term)
!

      implicit none

      integer,intent(in):: n,m
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n,5) :: A
!      real*8,dimension(m,n) :: D
      real*8,dimension(m,n) :: b
      real*8,dimension(m,n) :: x

!f2py intent(inplace)::msk,b
!f2py intent(inplace)::x

      real*8,dimension(m) :: y,d,ud,rhs

      integer:: i,j,im,ip,jm,jp

      do i=2,n-1
         ip=i+1
         im=i-1
         ! here we assume a halo width = 3
         ! since psi(j,i) is on the (top right) corner
         ! of cell (j,i)
         ! the bottom bry is x(4,:)=0
         ! the top    bry is x(m-3,:)=0
         ! the tridiagonal inversion is done
         ! on interior points 4 <= j <= m-4
         ! watch out Fortran starts with index 1...
         ! whereas Python starts with index 0
         if(msk(4,i).ne.0) then
         y(:)=0.
         do j=4,m-4
            jm=j-1
            jp=j+1

! y = b - A*x reads
!!$            y(j,i)=b(j,i)
!!$     $           -A(j,i,1)*x(j-1,i-1)
!!$     $           -A(j,i,2)*x(j-1,i  )
!!$     $           -A(j,i,3)*x(j-1,i+1)
!!$     $           -A(j,i,4)*x(j,i-1)
!!$     $           -A(j,i,5)*x(j,i  )
!!$     $           -A(j,i+1,4)*x(j,i+1)
!!$     $           -A(j+1,i-1,3)*x(j+1,i-1)
!!$     $           -A(j+1,i,2)*x(j+1,i  )
!!$     $           -A(j+1,i+1,1)*x(j+1,i+1)
            
            rhs(j)= b(j,i)
     $       -A(j,i,1)*x(jm,im)
     $       -A(j,i,3)*x(jm,ip)
     $       -A(j,i,4)*x(j,im)
     $       -A(j,ip,4)*x(j,ip)
     $       -A(jp,im,3)*x(jp,im)
     $       -A(jp,ip,1)*x(jp,ip)
     
            d(j) = A(j,i,5)
            ud(j) = A(jp,i,2)
         enddo
         call tridiag(d,ud,rhs,y,m)
         do j=4,m-4
            x(j,i)=y(j)
         enddo
         endif
      enddo

      end subroutine
!----------------------------------------
      subroutine tridiag(d,dd,b,xc,l)
        !     Axc = b
        !     Solve tridiagonal system
        implicit none
        integer,intent(in)  :: l
        real*8,dimension(l) :: d,b
        real*8,dimension(l) :: dd
        real*8,dimension(l) :: xc
        
!f2py intent(inplace)::d,dd,b,xc

        integer:: k,k0
        real*8,dimension(l):: gam
        real*8             :: bet

        k0 = 4
        bet   = 1. / d(k0)
        xc(k0) = b(k0)*bet
        do k=k0+1,l-4
           gam(k)= dd(k-1)*bet
           bet     = 1. /(d(k)-dd(k-1)*gam(k))
           xc(k) = (b(k)-dd(k-1)*xc(k-1))*bet
        enddo
        do k=l-5,k0,-1
           xc(k) = xc(k)-gam(k+1)*xc(k+1)
        enddo

      end subroutine tridiag

!----------------------------------------
      subroutine computeresidualwithA(msk,A,x,b,y,m,n)

! compute y=b-A*x

      implicit none

      integer,intent(in):: n,m
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n,5) :: A
      real*8,dimension(m,n) :: x,b
      real*8,dimension(m,n) :: y

!f2py intent(inplace)::msk,A,x,b
!f2py intent(inplace)::y

!      real*8,dimension(n,3) :: xo
      integer:: i,j,km,ko,kp

      do j=2,m-1
         do i=2,n-1
            if(msk(j,i).ne.0)then
            y(j,i)=b(j,i)
     $           -A(j,i,1)*x(j-1,i-1)
     $           -A(j,i,2)*x(j-1,i  )
     $           -A(j,i,3)*x(j-1,i+1)
     $           -A(j,i,4)*x(j,i-1)
     $           -A(j,i,5)*x(j,i  )
     $           -A(j,i+1,4)*x(j,i+1)
     $           -A(j+1,i-1,3)*x(j+1,i-1)
     $           -A(j+1,i,2)*x(j+1,i  )
     $           -A(j+1,i+1,1)*x(j+1,i+1)
!!$     $           -A(j,i,6)*x(j,i+1)
!!$     $           -A(j,i,7)*x(j+1,i-1)
!!$     $           -A(j,i,8)*x(j+1,i  )
!!$     $           -A(j,i,9)*x(j+1,i+1)
             else
                y(j,i)=0.
             endif
         enddo
      enddo


      end subroutine

!----------------------------------------
      subroutine fillhalo(x,nh,m,n)

! compute y=sum(x*x)

      implicit none

      integer,intent(in):: n,m,nh
      real*8,dimension(m,n) :: x

!f2py intent(inplace)::x

      integer:: i,j,jj,n2,m2

      n2=n-2*nh
      m2=m-2*nh
      do j=1,nh
         do i=1,nh
            x(j,i)=x(m2+j,n2+i)
         enddo
         do i=1,n2
            x(j,i+nh)=x(m2+j,i+nh)
         enddo
         do i=1,nh
            x(j,i+n-nh)=x(m2+j,nh+i)
         enddo
      enddo
      do j=1,m2
         jj=j+nh
         do i=1,nh
            x(jj,i)=x(jj,n2+i)
            x(jj,i+n-nh)=x(jj,i+nh)
         enddo
      enddo
      do j=1,nh
         jj=j+m2+nh
         do i=1,nh
            x(jj,i)=x(j+nh,n2+i)
         enddo
         do i=1,n-2*nh
            x(jj,i+nh)=x(j+nh,i+nh)
         enddo
         do i=1,nh
            x(jj,i+n-nh)=x(j+nh,i+nh)
         enddo
      enddo


      end subroutine

!----------------------------------------
      subroutine interpolate(msk1,msk2,x2,nh,n2,m2,n1,m1,x1)

! compute x1=I*x2 where x2 is on the coarse grid, x1 is on the fine grid

      implicit none

      integer,intent(in):: n1,m1,n2,m2,nh
      integer*1,dimension(m1,n1) :: msk1
      integer*1,dimension(m2,n2) :: msk2
      real*8,dimension(m2,n2) :: x2
      real*8,dimension(m1,n1) :: x1

!f2py intent(inplace)::x2,msk1,msk2
!f2py intent(inplace)::x1

!      real*8,dimension(n2,2) :: xo
      integer:: i1,j1,i2,j2,ko,kp
      real*8,dimension(0:1)::c1
      real*8,dimension(0:2)::c2
      real*8,dimension(0:4)::c3
      integer*1::s

      data c1/0.,1./
      data c2/0.,1.,0.5/
      data c3/0.,1.,0.5,0.3333333333333333333333333333,0.25/

!!$      data c1/0.,1./
!!$      data c2/0.,0.5,0.5/
!!$      data c3/0.,0.25,0.25,0.25,0.25/

!      xo(:,1)=x2(2,:)
      ko=1
      kp=2
!      j1=nh 
      if (nh.eq.2)j1=2
      if (nh.eq.3)j1=1
      ! mapping is (i1,j1)=(2,2) coincides with (i2,j2)=(2,2)

!      do j2=nh,m2-nh
      do j2=2,m2-2
!         i1=nh
         if (nh.eq.2)i1=2
         if (nh.eq.3)i1=1
!         xo(:,kp)=x2(j2+1,:)
!         x1(j1,:)=0.
!         x1(j1+1,:)=0.
!         do i2=nh,n2-nh
         do i2=2,n2-2
            !s=msk2(j2,i2)
            if (msk1(j1,i1).gt.0)then
               x1(j1,i1) = x2(j2,i2)
            else
               x1(j1,i1) = 0.
            endif

            if (msk1(j1,i1+1).gt.0) then    
               s=msk2(j2,i2)+msk2(j2,i2+1)
               x1(j1,i1+1) = (x2(j2,i2)+x2(j2,i2+1))*c2(s)!*c2(s)!*0.5!*msk1(j1,i1+1)*c2(s)
            else
               x1(j1,i1+1) = 0.
            endif

            if (msk1(j1+1,i1).gt.0) then
               s=msk2(j2,i2)+msk2(j2+1,i2)
               x1(j1+1,i1) = (x2(j2,i2)+x2(j2+1,i2))*c2(s)!*0.5!msk1(j1+1,i1)*c2(s)
            else
               x1(j1+1,i1) = 0.
            endif

            if (msk1(j1+1,i1+1).gt.0) then    
               s=msk2(j2,i2)+msk2(j2,i2+1)+msk2(j2+1,i2)+msk2(j2+1,i2+1)
               x1(j1+1,i1+1) = c3(s)*(x2(j2,i2)+x2(j2,i2+1)
     $                     +x2(j2+1,i2)+x2(j2+1,i2+1))!*msk1(j1+1,i1+1)
             else
                x1(j1+1,i1+1) = 0.
             endif
               i1=i1+2
         enddo
         j1=j1+2
!         ko=mod(ko,2)+1
!         kp=mod(kp,2)+1
      enddo

      end subroutine

!----------------------------------------
      subroutine restrict(msk2,x1,nh,n2,m2,n1,m1,x2)

! compute x1=I*x2 where x2 is on the coarse grid, x1 is on the fine grid

      implicit none

      integer,intent(in):: n1,m1,n2,m2,nh
      integer*1,dimension(m2,n2) :: msk2
      real*8,dimension(m2,n2) :: x2
      real*8,dimension(m1,n1) :: x1

!f2py intent(inplace)::x1,msk2
!f2py intent(inplace)::x2

!      real*8,dimension(n1,3) :: xo
      integer:: i1,j1,i2,j2,km,ko,kp
      real*8:: w

      km=1
      ko=2
      kp=3
      j1=nh
!      xo(:,1)=x1(j1-1,:)
      ! mapping is (i1,j1)=(2,2) coincides with (i2,j2)=(2,2)
      do j2=nh,m2-nh
         i1=nh
!         xo(:,km)=x1(j1-1,:)
!         xo(:,ko)=x1(j1,:)
!         xo(:,kp)=x1(j1+1,:)
         do i2=nh,n2-nh
            if(msk2(j2,i2).ne.0)then
               x2(j2,i2) = 0.25*x1(j1,i1)
     $ +0.125 *(x1(j1  ,i1-1)+x1(j1  ,i1+1)+x1(j1-1,i1)  +x1(j1+1,i1)  )
     $ +0.0625*(x1(j1-1,i1-1)+x1(j1-1,i1+1)+x1(j1+1,i1-1)+x1(j1+1,i1+1))
            else
              x2(j2,i2) = 0.
           endif
            i1=i1+2
         enddo
         j1=j1+2
         ! km=mod(km+1,3)+1
         ! ko=mod(ko+1,3)+1
         ! kp=mod(kp+1,3)+1
      enddo

      end subroutine

!----------------------------------------
      subroutine halotobuffer(x,b0,b1,b2,b3,b4,b5,b6,b7,nh,n,m,n2,m2)

! copy halo of array x into buffers
! x is assumed in 2D shape, Fortran order

      implicit none

      integer,intent(in):: n,m,nh,n2,m2
      real*8,dimension(m,n) :: x
      ! inout is important because buffers are preallocated
      real*8,dimension(nh,nh) :: b0,b2,b5,b7
      real*8,dimension(nh,n2) :: b1,b6
      real*8,dimension(m2,nh) :: b3,b4

!f2py intent(inplace,c)::b0,b2,b5,b7
!f2py intent(inplace,c)::b1,b6
!f2py intent(inplace,c)::b3,b4
!f2py intent(inplace)::x

      integer:: i,j,jj
      
      do j=1,nh
         do i=1,nh
            b0(j,i)=x(j+nh,i+nh)
         enddo
      enddo
      do j=1,nh
         do i=1,n2
            b1(j,i)=x(j+nh,i+nh)
         enddo
      enddo
      do j=1,nh
         do i=1,nh
            b2(j,i)=x(j+nh,i+n2)
         enddo
      enddo
      do j=1,m2
         jj=j+nh
         do i=1,nh
            b3(j,i)=x(jj,i+nh)
         enddo
      enddo
      do j=1,m2
         jj=j+nh
         do i=1,nh
            b4(j,i)=x(jj,i+n2)
         enddo
      enddo
      do j=1,nh
         jj=j+m2
         do i=1,nh
            b5(j,i)=x(jj,i+nh)
         enddo
      enddo
      do j=1,nh
         jj=j+m2
         do i=1,n-2*nh
            b6(j,i)=x(jj,i+nh)
         enddo
      enddo
      do j=1,nh
         jj=j+m2
         do i=1,nh
            b7(j,i)=x(jj,i+n-2*nh)
         enddo
      enddo

      end subroutine

!----------------------------------------
      subroutine buffertohalo(x,b0,b1,b2,b3,b4,b5,b6,b7,nh,n,m,n2,m2)

! copy halo of array x into buffers
! x is assumed in 2D shape, Fortran order

      implicit none

      integer,intent(in):: n,m,nh,n2,m2
      real*8,dimension(m,n) :: x
      ! inout is important because buffers are preallocated
      real*8,dimension(nh,nh) :: b0,b2,b5,b7
      real*8,dimension(nh,n2) :: b1,b6
      real*8,dimension(m2,nh) :: b3,b4

!f2py intent(inplace,c)::b0,b1,b2,b3,b4,b5,b6,b7
!f2py intent(inplace)::x

      integer:: i,j,jj
      
      do j=1,nh
         do i=1,nh
            x(j,i)=b7(j,i)
         enddo
         do i=1,n-2*nh
            x(j,i+nh)=b6(j,i)
         enddo
         do i=1,nh
            x(j,i+n-nh)=b5(j,i)
         enddo
      enddo

      do j=1,m2
         jj=j+nh
         do i=1,nh
            x(jj,i)=b4(j,i)
            x(jj,i+n-nh)=b3(j,i)
         enddo
      enddo
      do j=1,nh
         jj=j+m-nh
         do i=1,nh
            x(jj,i)=b2(j,i)
         enddo
         do i=1,n2
            x(jj,i+nh)=b1(j,i)
         enddo
         do i=1,nh
            x(jj,i+n-nh)=b0(j,i)
         enddo
      enddo

      end subroutine

!----------------------------------------
      subroutine buffertodomain(b,x,nh,m1,n1,m,n,mp,np)

! copy buffer

      implicit none

      integer,intent(in):: n1,m1,n,m,nh,np,mp
!      real*8,dimension(2*nh+np*(n-2*nh),2*nh+mp*(m-2*nh)),intent(out) :: x
      real*8,dimension(m1,n1) :: x
      real*8,dimension(mp,np,m,n) :: b

!f2py intent(inplace)::b
!f2py intent(inplace)::x

      integer:: i,j,k,l,ii,jj

      do l=1,mp
         jj=1+(l-1)*(m-2*nh)
         do j=1,m
            do k=1,np
               ii=1+(k-1)*(n-2*nh)
               do i=1,n                  
                  x(jj,ii) = b(l,k,j,i)
                  ii=ii+1
               enddo
            enddo
            jj=jj+1
         enddo
      enddo

      end subroutine

!----------------------------------------
      subroutine coarsenmatrix(Afine,Acoarse,msk1,msk2,nh,n1,m1,n2,m2)

! compute the matrix A on the coarse grid

      implicit none

      integer,intent(in):: n1,m1,n2,m2,nh
!      real*8,dimension(2*nh+np*(n-2*nh),2*nh+mp*(m-2*nh)),intent(out) :: x
      real*8,dimension(m2,n2,9) :: Acoarse
      real*8,dimension(m1,n1,9) :: Afine
      integer*1,dimension(m1,n1):: msk1
      integer*1,dimension(m2,n2):: msk2

!f2py intent(out)::Acoarse
!f2py intent(in)::Afine,msk1,msk2

      real*8::w
      real*8,dimension(-1:1,-1:1)::coef,z3
      real*8,dimension(-2:2,-2:2)::z5
      real*8,dimension(0:3):: w2
      integer:: i,j,k,l
      integer::i1,j1,i2,j2,ii,jj,di2,dj2
      integer::ki,kj

      integer::loc,m

      data coef/0.125,0.25,0.125,0.25,0.5,0.25,0.125,0.25,0.125/
      data w2/1,0.5,0.5,0.25/

      do j2=nh,m2-nh
         do i2=nh,n2-nh
            if(msk2(j2,i2).eq.1)then
            do l=1,9 ! loop on each 9 points to be connected with
               ! (i2,j2) location on the coarse grid
               ! (i1,j1) location on the fine grid
               i1= 2*(i2-nh)+nh
               j1= 2*(j2-nh)+nh
               !
               di2=mod(l-1,3)-1
               dj2=   (l-1)/3-1
               ! fill the 5x5 matrix (interpolation phase)
               z5(:,:)=0.
               do kj=-1,1
                  do ki=-1,1
                     i=2*di2+ki
                     j=2*dj2+kj
                     if ((abs(i).le.2).and.(abs(j).le.2))then
                     loc=mod(i+2,2)+mod(j+2,2)*2
!                     write(*,*)i,j,loc
                       if (msk1(j1+j,i1+i).eq.1)
     $                   z5(i,j)=2.*coef(ki,kj)
!                        z5(i,j)=2.*coef(ki,kj)
                        if(msk1(j1+j,i1+i).eq.1)then
                           if(loc.eq.0)
     $       m=msk2(j2+dj2,i2+di2)
                           if(loc.eq.1)
     $       m=(msk2(j2+dj2,i2+di2)+msk2(j2+dj2,i2+di2+1))
                           if(loc.eq.2)
     $       m=(msk2(j2+dj2,i2+di2)+msk2(j2+dj2+1,i2+di2))
                           if(loc.eq.3)
     $       m=(msk2(j2+dj2  ,i2+di2)+msk2(j2+dj2  ,i2+di2+1)
     $         +msk2(j2+dj2+1,i2+di2)+msk2(j2+dj2+1,i2+di2+1))
!                          if(m.gt.0) z5(i,j)=1./m
!                          z5(i,j)=w2(loc)!*msk2(j2+dj2,i2+di2)
!                          m=m-1./w2(loc)
!                          if(m.ne.0)write(*,*)i1,j1,m
                        endif
                     endif
                  enddo
               enddo
               ! compute the laplacian on the 3x3 matrix
               z3(:,:)=0.
               do jj=-1,1
                  do ii=-1,1
                     do kj=-1,1
                        do ki=-1,1
                           k=1+(ki+1)+(kj+1)*3               
!                           if((i1+ii+ki.gt.0).and.(j1+jj+kj.gt.0))
                           if(msk1(j1+jj,i1+ii).eq.1)
     $                     z3(ii,jj) = z3(ii,jj) 
     $                    + Afine(j1+jj,i1+ii,k)*z5(ii+ki,jj+kj)
                        enddo
                     enddo
                  enddo
               enddo
               ! now do the restriction using the 3x3 matrix
               w=0.
               if (msk2(j2,i2).ne.0) then
                  do jj=-1,1
                     do ii=-1,1
                        w = w + 0.5*coef(ii,jj)*z3(ii,jj)
                     enddo
                  enddo
               endif
               Acoarse(j2,i2,l) = w
            enddo
            else
               do l=1,9
                  Acoarse(j2,i2,l) = 0.
               enddo
            endif
         enddo
      enddo
      
      
      end subroutine
!----------------------------------------
      subroutine computenorm(msk,x,nh,n,m,y)

! compute y=sum(x*x)
!
 
      implicit none

      integer,intent(in):: n,m,nh
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n) :: x
      real*8,intent(out) :: y

!f2py intent(inplace)::x,msk
!f2py intent(out)::y

      integer:: i,j

      y=0.
      do j=nh+1,m-nh
         do i=nh+1,n-nh
            if(msk(j,i).ne.0)then
               y = y + x(j,i)**2
            endif
         enddo
      enddo

      end subroutine
!----------------------------------------
      subroutine computeinner(msk,x,y,nh,n,m,z)

! compute z=sum(x*y)
!
 
      implicit none

      integer,intent(in):: n,m,nh
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n) :: x,y
      real*8,intent(out) :: z

!f2py intent(inplace)::x,y,msk
!f2py intent(out)::z

      integer:: i,j

      z=0.
      do j=nh+1,m-nh
         do i=nh+1,n-nh
            if(msk(j,i).ne.0)then
               z = z + x(j,i)*y(j,i)
            endif
         enddo
      enddo

      end subroutine

