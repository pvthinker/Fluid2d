!----------------------------------------
  subroutine adv_upwind(msk,x,y,u,v,cst,nh,method,order,m,n)

    ! method 0: min/max, 1: parabolic, 2: logcosh

    implicit none

    integer,intent(in):: n,m,nh,method,order
    integer*1,dimension(m,n) :: msk
    real*8,dimension(m,n) :: x,y,u,v
    real*8,dimension(5):: cst

!f2py intent(inplace)::msk,x,y,u,v

    ! fluxes (fym is the flux through the bottom edge)
    real*8,dimension(n) :: fx,fy,fym
    
    integer*1:: mx5,mx3,mx2,my2
    integer*1,dimension(n)::my5,my3

    real*8::zdx,dx,zdy,dy,a,logcosh,umax,Uparab,u1,aa,bb,aparab

    real*8::KK,UU,up,um,qp,qm

    real*8::d1,d2,d3,d4,d5
    real*8::c1,c2,c3
    
    integer::i,j

    if (nh.ne.3)then
       write(*,*) "NHALO = 3 is compulsory with UP5"
       write(*,*) "catastrophic ABORT!!!!"
       stop
    endif

    d1 =   1./30.
    d2 = -13./60.
    d3 =  47./60.
    d4 =   9./20.
    d5 =  -1./20.

    c1 = -1./6.
    c2 =  5./6.
    c3 =  2./6.

    dx      = cst(1)
    dy      = cst(2)
    logcosh = cst(3)
    umax    = cst(4)
    aparab  = cst(5)

    zdx = 1./dx
    zdy = 1./dy
    !
    a   = logcosh/umax
    KK   = ( umax - log(cosh(umax * a))/a) * 0.5
    !
    u1 = aparab * umax
    aa = 1./(2.*u1)
    bb = u1*0.5
    !
    do i=nh,n-nh
       my5(i)=msk(1,i)+msk(2,i)+msk(3,i)+msk(4,i)+msk(5,i)
       my3(i)=         msk(2,i)+msk(3,i)+msk(4,i)
    enddo

    ! main and single double do-loop
    do j=nh,m-nh
       mx5=msk(j,1)+msk(j,2)+msk(j,3)+msk(j,4)+msk(j,5)
       mx3=         msk(j,2)+msk(j,3)+msk(j,4)
       ! first compute the fluxes across cells
       do i=nh,n-nh
          mx2=msk(j,i)+msk(j,i+1)
          if (mx2.eq.2)then
             
             ! flux splitting method
             if (method.eq.0) UU = abs(u(j,i))
             if (method.eq.1) then
                UU = abs(u(j,i))
                if(UU.lt.u1) UU = aa*u(j,i)**2+bb
             endif
             ! the overhead of logcosh is GIGANTIC compared to parabolic
             if (method.eq.2) UU = 2*KK+log(cosh(u(j,i) * a)) / a

             up=0.5*(u(j,i)+UU)
             um=0.5*(u(j,i)-UU)
             ! right-going x flux
             if ((mx5.eq.5).and.(order.eq.5)) then
                qp=d1*x(j,i-2)+d2*x(j,i-1)+d3*x(j,i)+d4*x(j,i+1)+d5*x(j,i+2)
             else if ((mx3.eq.3).and.(order.ge.3)) then
                qp=c1*x(j,i-1)+c2*x(j,i)+c3*x(j,i+1)
             else
                qp=x(j,i)
             endif
             ! left-going x flux
             mx5=mx5-msk(j,i-2)+msk(j,i+3)
             mx3=mx3-msk(j,i-1)+msk(j,i+2)
             if ((mx5.eq.5).and.(order.eq.5))  then
                qm=d5*x(j,i-1)+d4*x(j,i)+d3*x(j,i+1)+d2*x(j,i+2)+d1*x(j,i+3)
             else if ((mx3.eq.3).and.(order.ge.3)) then
                qm=c3*x(j,i)+c2*x(j,i+1)+c1*x(j,i+2)
             else
                qm=x(j,i+1)
             endif                
             fx(i)=up*qp+um*qm
             ! flux is the sum of right and left-going
          else
             mx5=mx5-msk(j,i-2)+msk(j,i+3)
             mx3=mx3-msk(j,i-1)+msk(j,i+2)
             fx(i)=0.
          endif
          ! same but for the along y flux
          my2=msk(j,i)+msk(j+1,i)
          if (my2.eq.2)then
             !
             if (method.eq.0) UU = abs(v(j,i))
             if (method.eq.1) then
                UU = abs(v(j,i))
                if(UU.lt.u1) UU  = aa*v(j,i)**2+bb
             endif
             if (method.eq.2) UU = 2*KK+log(cosh(v(j,i) * a)) / a
             !
             up=0.5*(v(j,i)+UU)
             um=0.5*(v(j,i)-UU)
             if ((my5(i).eq.5).and.(order.eq.5)) then
                qp=d1*x(j-2,i)+d2*x(j-1,i)+d3*x(j,i)+d4*x(j+1,i)+d5*x(j+2,i)
             else if ((my3(i).eq.3).and.(order.ge.3)) then
                qp=c1*x(j-1,i)+c2*x(j,i)+c3*x(j+1,i)
             else
                qp=x(j,i)
             endif
             my5(i)=my5(i)-msk(j-2,i)+msk(j+3,i)
             my3(i)=my3(i)-msk(j-1,i)+msk(j+2,i)
             if ((my5(i).eq.5).and.(order.eq.5)) then
                qm=d5*x(j-1,i)+d4*x(j,i)+d3*x(j+1,i)+d2*x(j+2,i)+d1*x(j+3,i)
             else if ((my3(i).eq.3).and.(order.ge.3)) then
                qm=c3*x(j,i)+c2*x(j+1,i)+c1*x(j+2,i)
             else
                qm=x(j+1,i)
             endif
             fy(i)=up*qp+um*qm
          else
             my5(i)=my5(i)-msk(j-2,i)+msk(j+3,i)
             my3(i)=my3(i)-msk(j-1,i)+msk(j+2,i)
             fy(i)=0.
          endif
       enddo
       ! fluxes are ok we can move to the divergence
       if (j.gt.nh)then
          do i=nh+1,n-nh
             ! divergence of the flux
             y(j,i) = -zdx*(fx(i)-fx(i-1))-zdy*(fy(i)-fym(i))
             !swap
             fym(i)=fy(i)
          enddo
       else
          do i=nh+1,n-nh
             !simply swap
             fym(i)=fy(i)
          enddo
       endif
       
    enddo
    
  end subroutine adv_upwind


!----------------------------------------
  subroutine adv_centered(msk,x,y,u,v,cst,nh,method,order,m,n)

    ! method 0: min/max, 1: parabolic, 2: logcosh

    implicit none

    integer,intent(in):: n,m,nh,method,order
    integer*1,dimension(m,n) :: msk
    real*8,dimension(m,n) :: x,y,u,v
    real*8,dimension(5):: cst

!f2py intent(inplace)::msk,x,y,u,v

    ! fluxes (fym is the flux through the bottom edge)
    real*8,dimension(n) :: fx,fy,fym
    
    integer*1:: mx6,mx4,mx2
    integer*1,dimension(n)::my6,my4,my2

    real*8::zdx,dx,zdy,dy,a,logcosh,umax,Uparab,u1,aa,bb

    real*8::KK,UU,up,um,qp,qm

    real*8::e1,e2,e3
    real*8::d1,d2
    real*8::c1
    
    integer::i,j

    if (nh.ne.3)then
       write(*,*) "NHALO = 3 is compulsory with UP5"
       write(*,*) "catastrophic ABORT!!!!"
       stop
    endif

    e1 = 1./60.
    e2 = -2./15.
    e3 = 37./60.
    
    d1 = -1./12.
    d2 =  7./12.

    c1 = 1./2.

    dx      = cst(1)
    dy      = cst(2)

    zdx = 1./dx
    zdy = 1./dy
    !
    do i=nh,n-nh
       my6(i)=msk(1,i)+msk(2,i)+msk(3,i)+msk(4,i)+msk(5,i)+msk(6,i)
       my4(i)=         msk(2,i)+msk(3,i)+msk(4,i)+msk(5,i)
       my2(i)=                  msk(3,i)+msk(4,i)
    enddo

    ! main and single double do-loop
    do j=nh,m-nh
       ! first compute the fluxes across cells
       do i=nh,n-nh
          mx6=msk(j,i-2)+msk(j,i-1)+msk(j,i)+msk(j,i+1)+msk(j,i+2)
          mx6=mx6+msk(j,i+3)
          mx4=msk(j,i-1)+msk(j,i)+msk(j,i+1)+msk(j,i+2)
          mx2=           msk(j,i)+msk(j,i+1)
          if (mx2.eq.2)then             
             ! x flux
             if ((mx6.eq.6).and.(order.eq.6)) then
                qp=e1*(x(j,i-2)+x(j,i+3))+e2*(x(j,i-1)+x(j,i+2))
                qp=qp+e3*(x(j,i)+x(j,i+1))
             else if ((mx4.eq.4).and.(order.eq.4)) then
                qp=d1*(x(j,i-1)+x(j,i+2))+d2*(x(j,i)+x(j,i+1))
             else if ((mx2.eq.2).and.(order.ge.2)) then
                qp=c1*(x(j,i)+x(j,i+1))
             endif
             fx(i)=u(j,i)*qp
             ! flux is the sum of right and left-going
          else
             fx(i)=0.
          endif
          ! same but for the along y flux
          if (my2(i).eq.2)then
             !
             if ((my6(i).eq.6).and.(order.eq.6)) then
                qp=e1*(x(j-2,i)+x(j+3,i))
                qp=qp+e2*(x(j-1,i)+x(j+2,i))+e3*(x(j,i)+x(j+1,i))
             else if ((my4(i).eq.4).and.(order.eq.4)) then
                qp=d1*(x(j-1,i)+x(j+2,i))+d2*(x(j,i)+x(j+1,i))
             else if ((my2(i).eq.2).and.(order.ge.2)) then
                qp=c1*(x(j,i)+x(j+1,i))
             endif
             fy(i)=v(j,i)*qp
          else
             fy(i)=0.
          endif
          my6(i)=my6(i)-msk(j-2,i)+msk(j+4,i)
          my4(i)=my4(i)-msk(j-1,i)+msk(j+3,i)
          my2(i)=my2(i)-msk(j  ,i)+msk(j+2,i)
       enddo
       ! fluxes are ok we can move to the divergence
       if (j.gt.nh)then
          do i=nh+1,n-nh
             ! divergence of the flux
             y(j,i) = -zdx*(fx(i)-fx(i-1))-zdy*(fy(i)-fym(i))
             !swap
             fym(i)=fy(i)
          enddo
       else
          do i=nh+1,n-nh
             !simply swap
             fym(i)=fy(i)
          enddo
       endif
       
    enddo
    
  end subroutine adv_centered

