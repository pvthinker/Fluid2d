
!----------------------------------------
      subroutine computedotprod(msk,x,y,nh,n,m,z)

! compute y=sum(x*x)
!
! TODO: add the mask = the norm is computed only on the fluid points
 
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

!----------------------------------------
      subroutine computemax(msk,x,nh,n,m,y)

! compute y=sum(x*x)
!
! TODO: add the mask = the norm is computed only on the fluid points
 
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
               y = max(y,abs(x(j,i)))
            endif
         enddo
      enddo

      end subroutine

!----------------------------------------
      subroutine computesum(msk,x,nh,n,m,y)

! compute y=sum(x)
!
! TODO: add the mask = the norm is computed only on the fluid points
 
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
               y = y + x(j,i)
            endif
         enddo
      enddo

      end subroutine

!----------------------------------------
      subroutine computesumandnorm(msk,x,nh,n,m,y,y2)

! compute y=sum(x) and y2=sum(x*x)
!
 
      implicit none

      integer,intent(in):: n,m,nh
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n) :: x
      real*8,intent(out) :: y,y2

!f2py intent(inplace)::x,msk
!f2py intent(out)::y,y2

      integer:: i,j

      y=0.
      y2=0.
      do j=nh+1,m-nh
         do i=nh+1,n-nh
            if(msk(j,i).eq.1)then
               y = y + x(j,i)
               y2 = y2 + x(j,i)**2
            endif
         enddo
      enddo
      end subroutine

!----------------------------------------
      subroutine computenormmaxu(msk,x,nh,n,m,y,ymax)

! compute y=sum(x*x)
!
! TODO: add the mask = the norm is computed only on the fluid points
 
      implicit none

      integer,intent(in):: n,m,nh
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n) :: x
      real*8,intent(out) :: y,ymax

!f2py intent(inplace)::x,msk
!f2py intent(out)::y,ymax

      integer:: i,j

      y=0.
      ymax=0.
      do j=nh+1,m-nh
         do i=nh+1,n-nh
            if( (msk(j,i)+msk(j,i+1)).eq.2)then
               y = y + x(j,i)**2
               ymax = max(ymax,abs(x(j,i)))
            endif
         enddo
      enddo

      end subroutine

!----------------------------------------
      subroutine computekemaxu(msk,u,v,nh,n,m,ke,maxu)

! compute y=sum(x*x)
!
! TODO: add the mask = the norm is computed only on the fluid points
 
      implicit none

      integer,intent(in):: n,m,nh
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n) :: u,v
      real*8,intent(out) :: ke,maxu
      real*8:: z,zu,zv,maxv,um,vm

!f2py intent(inplace)::u,v,msk
!f2py intent(out)::ke,maxu

      integer:: i,j

      ke=0.
      maxu=0.
      maxv=0.
      do j=nh+1,m-nh
         do i=nh+1,n-nh
            if (msk(j,i).eq.1)then
               !z = (u(j,i)+u(j,i-1))**2+(v(j,i)+v(j-1,i))**2
               zu = u(j,i)**2+u(j,i-1)**2
               zv = v(j,i)**2+v(j-1,i)**2
               um = abs(u(j,i)+u(j,i-1))
               vm = abs(v(j,i)+v(j-1,i))
               ke = ke + zu+zv
               !maxu = max(maxu,zu)!zu+zv)
               !maxv = max(maxv,zv)!zu+zv)
               maxu = max(maxu,um+vm)
            endif
         enddo
      enddo
      ke = ke * 0.25
      maxu = maxu*0.5

      end subroutine

!----------------------------------------
      subroutine computekemaxuv(msk,u,v,nh,n,m,ke,maxu,maxv)

! compute y=sum(x*x)
!
! TODO: add the mask = the norm is computed only on the fluid points
 
      implicit none

      integer,intent(in):: n,m,nh
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n) :: u,v
      real*8,intent(out) :: ke,maxu,maxv
      real*8:: z,zu,zv

!f2py intent(inplace)::u,v,msk
!f2py intent(out)::ke,maxu,maxv

      integer:: i,j

      ke=0.
      maxu=0.
      maxv=0.
      do j=nh+1,m-nh
         do i=nh+1,n-nh
            if (msk(j,i).eq.1)then
               !z = (u(j,i)+u(j,i-1))**2+(v(j,i)+v(j-1,i))**2
               zu = u(j,i)**2+u(j,i-1)**2
               zv = v(j,i)**2+v(j-1,i)**2
               ke = ke + zu+zv
               maxu = max(maxu,zu)!zu+zv)
               maxv = max(maxv,zv)!zu+zv)
            endif
         enddo
      enddo
      ke = ke * 0.25
      maxu=sqrt(maxu/2.)
      maxv=sqrt(maxv/2.)

      end subroutine

!----------------------------------------
      subroutine computekewithpsi(msk,omega,psi,nh,n,m,ke)

! compute y=sum(x*x)
!
! TODO: add the mask = the norm is computed only on the fluid points
 
      implicit none

      integer,intent(in):: n,m,nh
      integer*1,dimension(m,n) :: msk
      real*8,dimension(m,n) :: omega,psi
      real*8,intent(out) :: ke

!f2py intent(inplace):: omega,psi,msk
!f2py intent(out)::ke

      integer:: i,j

      ke=0.
      do j=nh+1,m-nh
         do i=nh+1,n-nh
            if (msk(j,i).eq.1) then
               ke = ke - 0.125*(psi(j,i)+psi(j-1,i)
     & +psi(j-1,i-1)+psi(j,i-1))*omega(j,i)
            endif
         enddo
      enddo

      end subroutine

