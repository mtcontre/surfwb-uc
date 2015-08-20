subroutine bcs(fopt,Cf,Coef,pasoRK,t,dt,Fr2,caso,qn,z,&
  Nx,Ny,CB,xi,eta,Jac,dxi,deta,qt,xit,etat,zt)
  use couplingbc
  use senales
  use coords
  use global_variables,only: treal,tinit,tfinal !to communicate metrics+z
  use mpi_surf,only:myleft,myright,myback,myfront,myrank,coords2
  use mpi
  implicit none
  integer :: Nx,Ny, caso,nt1,nt2
  real (kind=8), dimension(3,Nx,Ny)	::qn
  real (kind=8), dimension(2,Nx,Ny)	::xi,eta
  real (kind=8), dimension(3,Nx+4,Ny+4)::qt
  real (kind=8), dimension(2,Nx+4,Ny+4)::xit,etat
  real (kind=8), dimension(Nx+4,Ny+4)	::zt
  real (kind=8), dimension(Nx,Ny)	::z, Coef, Jac
  
  integer,dimension(4) :: CB  
  real (kind=8)	:: t,dt,Fr2,dxi,deta,Zx1,Zxi1, Zeta1, Zx0, Zxi0,&
    Zeta0, ZxiN, ZetaN, ZxN, ZxN1, ZxiN1, ZetaN1, Zy1, Zy0, ZyN, ZyN1, hF
  integer::i,j,fopt,pasoRK, Cf
  !Paso RK sirve solo para cuando se usa GA
  real(kind=8), dimension(3,Ny):: qA10,qA20
  real(kind=8), dimension(2,Ny)::hfija

  ! 2nd Order Boundary conditions
  !
  !The four sides are numbered in the following way:
  !	1: xi=1
  !	2: xi=Nx
  !	3: eta=1
  !	4: eta=Ny
  !
  !Defination of boundary conditions type:
  !	1: solid wall
  !	2: periodic 
  !	3: exit
  !	4: generation/absorpsion condition*
  
  !0.Initialise of matrices

  DO i=1,Nx+4; DO j=1,Ny+4
    qt(:,i,j)=1.0D0
    zt(i,j)=1.0D0
    xit(1,i,j)=1.0D0
    xit(2,i,j)=1.0D0
    etat(1,i,j)=1.0D0
    etat(2,i,j)=1.0D0
  END DO; END DO
  
  !1. Complete Inside Cells
  DO i=3,Nx+2; DO j=3,Ny+2
    qt(:,i,j)=qn(:,i-2,j-2)
    xit(:,i,j)=xi(:,i-2,j-2)
    etat(:,i,j)=eta(:,i-2,j-2)
    zt(i,j)=z(i-2,j-2)
  END DO; END DO
  
!   if (myback==mpi_proc_null) then
    DO j=3,Ny+2
      xit(:,1,j)=xi(:,2,j-2)		!xi-1=xi2
      xit(:,2,j)=xi(:,1,j-2)		!xi0=xi12	    
      etat(:,1,j)=eta(:,2,j-2)		!eta-1=eta2
      etat(:,2,j)=eta(:,1,j-2)		!eta0=eta1	    
      zt(1,j)=z(2,j-2)
      zt(2,j)=z(1,j-2)
    END DO
!   end if
  
!   if (myfront==mpi_proc_null) then
    DO j=3,Ny+2	
      xit(:,Nx+4,j)=xi(:,Nx-1,j-2)	!xiNx+4=xiNx-1
      xit(:,Nx+3,j)=xi(:,Nx,j-2)		!xiNx+3=xiNx
      etat(:,Nx+4,j)=eta(:,Nx-1,j-2)	!etaNx+4=etaNx-1
      etat(:,Nx+3,j)=eta(:,Nx,j-2)	!etaNx+3=etaNx	  
      zt(Nx+4,j)=z(Nx-1,j-2)
      zt(Nx+3,j)=z(Nx,j-2)	
    END DO
!   end if
  
!   if (myleft==mpi_proc_null) then
    DO i=3,Nx+2	
      !Metricas
      xit(:,i,1)=xi(:,i-2,2)		!xi-1=xi2
      xit(:,i,2)=xi(:,i-2,1)		!xi0=xi12	  
      etat(:,i,1)=eta(:,i-2,2)		!eta-1=eta2
      etat(:,i,2)=eta(:,i-2,1)		!eta0=eta1	  
      zt(i,1)=z(i-2,2)
      zt(i,2)=z(i-2,1)	
    END DO
!   end if
  
!   if (myright==mpi_proc_null) then
    DO i=3,Nx+2	
      !Metricas
      xit(:,i,Ny+4)=xi(:,i-2,Ny-1)	!xiNy+4=xiNy-1
      xit(:,i,Ny+3)=xi(:,i-2,Ny)		!xiNy+3=xiNy
	    
      etat(:,i,Ny+4)=eta(:,i-2,Ny-1)	!etaNy+4=etaNy-1
      etat(:,i,Ny+3)=eta(:,i-2,Ny)	!etaNy+3=etaNy
	    
      zt(i,Ny+4)=z(i-2,Ny-1)
      zt(i,Ny+3)=z(i-2,Ny)
    END DO
!   end if
  
  call exchange_2d(qt,xit,etat,zt)

  SELECT CASE (CB(1))
    case(0) !Couplingb boundary condition
      !1 Interpolate
      nt1=min(int(t/dt_xi0g1),nt_xi0g1-1)!busca el menor entre el penultimo y int(t/dt)
      nt2=min(int(t/dt_xi0g2),nt_xi0g2-1)!la idea es encontrar el punto que esté más 
      !atrás e interpolar
      !two methods:
      ! piecewise constant:  f(x)=f([x])
      ! linear interpolated: f(x)=1/dx*(f(x_k)*(x-x_k)+f(x_{k+1})(x_{k+1}-x),k=0,..,nt-1
      !como indice de q (que va de 1 a nt+1), hay que sumarle uno!! 
      !(por eso qxi01(1,j-2,NT1+2))
      !nt1 va de 0 hasta nt1-1    
      if (optxi0g1.eq.0) then
	do j=3,Ny+2
	  zt(  1,j)=qxi0g1(1,j-2,nt1+1)!z
	  qt(1,1,j)=qxi0g1(2,j-2,nt1+1)!h
	  qt(2,1,j)=qxi0g1(3,j-2,nt1+1)!u
	  qt(3,1,j)=qxi0g1(4,j-2,nt1+1)!v
	end do
      else 
	do j=3,Ny+2
	  call interpj(nt1*dt_xi0g1,qxi0g1(1,j-2,nt1+1),(nt1+1)*dt_xi0g1,&
	    qxi0g1(1,j-2,nt1+2),t,zt(1,j) )
	  call interpj(nt1*dt_xi0g1,qxi0g1(2,j-2,nt1+1),(nt1+1)*dt_xi0g1,&
	    qxi0g1(2,j-2,nt1+2),t,qt(1,1,j))
	  call interpj(nt1*dt_xi0g1,qxi0g1(3,j-2,nt1+1),(nt1+1)*dt_xi0g1,&
	    qxi0g1(3,j-2,nt1+2),t,qt(2,1,j))
	  call interpj(nt1*dt_xi0g1,qxi0g1(4,j-2,nt1+1),(nt1+1)*dt_xi0g1,&
	    qxi0g1(4,j-2,nt1+2),t,qt(3,1,j))
	end do
      end if	
      if (optxi0g2.eq.0) then
	do j=3,Ny+2
	  zt(  2,j)=qxi0g2(1,j-2,nt2+1)
	  qt(1,2,j)=qxi0g2(2,j-2,nt2+1)
	  qt(2,2,j)=qxi0g2(3,j-2,nt2+1)
	  qt(3,2,j)=qxi0g2(4,j-2,nt2+1)
	end do
      else
	do j=3,Ny+2    
	  call interpj(nt2*dt_xi0g2,qxi0g2(1,j-2,nt2+1),(nt2+1)*dt_xi0g2,&
	    qxi0g2(1,j-2,nt2+2),t,zt(2,j) )
	  call interpj(nt2*dt_xi0g2,qxi0g2(2,j-2,nt2+1),(nt2+1)*dt_xi0g2,&
	    qxi0g2(2,j-2,nt2+2),t,qt(1,2,j))
	  call interpj(nt2*dt_xi0g2,qxi0g2(3,j-2,nt2+1),(nt2+1)*dt_xi0g2,&
	    qxi0g2(3,j-2,nt2+2),t,qt(2,2,j))
	  call interpj(nt2*dt_xi0g2,qxi0g2(4,j-2,nt2+1),(nt2+1)*dt_xi0g2,&
	    qxi0g2(4,j-2,nt2+2),t,qt(3,2,j))
	end do	 
      end if

    CASE(1) 
    !Solid Wall !Ojo que son las velocidades contravariantes las que se usan
    !pero si las metricas se asumen iguales entonces queda lo mismo que en cartesianas.	  
      DO j=3,Ny+2	  
	qt(1,1,j)=qn(1,2,j-2)		!h-1=h2
	qt(1,2,j)=qn(1,1,j-2)		!h0=h1
	  
	qt(2,1,j)=-qn(2,2,j-2)		!u-1=-u2
	qt(2,2,j)=-qn(2,1,j-2)		!u0=-u1
	  
	qt(3,1,j)=-qn(3,2,j-2)		!v-1=-v2
	qt(3,2,j)=-qn(3,1,j-2)		!v0=-v1	  
      END DO
      
    CASE(2) !Periodic	  
      !handled by communication
    CASE(3) !Free Exit	
      DO j=3,Ny+2	  
	qt(1,1,j)=qn(1,2,j-2)		!h-1=h2
	qt(1,2,j)=qn(1,1,j-2)		!h0=h1
	  
	qt(2,1,j)=qn(2,2,j-2)		!u-1=-u2
	qt(2,2,j)=qn(2,1,j-2)		!u0=-u1
	  
	qt(3,1,j)=qn(3,2,j-2)		!v-1=-v2
	qt(3,2,j)=qn(3,1,j-2)		!v0=-v1	  
      END DO	

	  
    CASE(4)
      IF (pasoRK==1.OR.pasoRK==3) then
	IF (GA1==1) THEN
	  call genabs0xi_1_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,etaL1,Nsenal1,&
	    h01,t,dt,qn,zt,xi,qA1,zA1)
	ELSE IF (GA1==2) THEN
	  call genabs0xi_2_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,qs1,hs1,h01,Nsenal1,&
	    t,dt,qn,zt,xi,qA1,zA1)
	ELSE IF (GA1==3) THEN
	  call genabs0xi_3_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,us1,hs1,h01,Nsenal1,&
	    t,dt,qn,zt,xi,qA1,zA1)
	ELSE IF (GA1==9) THEN
	  call genabs0xi_9_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,etaL9,timeS9,Nsenal1,&
	    h01,t,dt,qn,zt,xi,qA1,zA1)
	END IF
      END IF	
      qA10=qA1
      
      IF (pasoRK==2.OR.pasoRK==4) then
	IF (GA1==1) THEN
	  call genabs0xi_1_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,etaL1,Nsenal1,&
	    h01,t,dt,qn,zt,xi,qA10,qA1,zA1)
	ELSE IF (GA1==2) THEN
	  call genabs0xi_2_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,qs1,hs1,h01,Nsenal1,&
	    t,dt,qn,zt,xi,qA10,qA1,zA1)
	ELSE IF (GA1==3) THEN
	  call genabs0xi_3_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,us1,hs1,h01,Nsenal1,&
	    t,dt,qn,zt,xi,qA10,qA1,zA1)			
  	ELSE IF (GA1==9) THEN
	  if ((t+dt).gt.maxval(timeS9)) then
	    do j=3,Ny+2	
	      qt(1,1,j)=qn(1,2,j-2)		!h-1=h2
	      qt(1,2,j)=qn(1,1,j-2)		!h0=h1
	      qt(2,1,j)=qn(2,2,j-2)		!u-1=-u2
	      qt(2,2,j)=qn(2,1,j-2)		!u0=-u1	    
	      qt(3,1,j)=qn(3,2,j-2)		!v-1=-v2
	      qt(3,2,j)=qn(3,1,j-2)		!v0=-v1
	    end do	  
	  else
	    call genabs0xi_9_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,etaL9,timeS9,Nsenal1,h01,t,dt,qn,zt,xi,qA10,qA1,zA1)
	    DO j=3,Ny+2
	      qt(1,1,j)=qA1(1,j-2)		!h-1=h2
	      qt(1,2,j)=qA1(1,j-2)		!h0=h1
	      qt(2,1,j)=qA1(2,j-2)		!u-1=-u2
	      qt(2,2,j)=qA1(2,j-2)		!u0=-u1	
	      qt(3,1,j)=qA1(3,j-2)		!v-1=-v2
	      qt(3,2,j)=qA1(3,j-2)		!v0=-v1
	    END DO	
	  end if
  	END IF		
      END IF    	
  END SELECT
  
!---------------------------------------------------------------------------

  
  SELECT CASE (CB(2))
    case(0) !Customized boundary condition
      !1 Interpolate
      nt1=min(int(t/dt_xiNg1),nt_xiNg1-1)
      nt2=min(int(t/dt_xiNg2),nt_xiNg2-1)
      if (optxiNg1.eq.0) then
	do j=3,Ny+2
	  zt(Nx+3  ,j)=qxiNg1(1,j-2,nt1+1)
	  qt(1,Nx+3,j)=qxiNg1(2,j-2,nt1+1)
	  qt(2,Nx+3,j)=qxiNg1(3,j-2,nt1+1)  
	  qt(3,Nx+3,j)=qxiNg1(4,j-2,nt1+1)
	end do
      else
	do j=3,Ny+2
	  call interpj(nt1*dt_xiNg1,qxiNg1(1,j-2,nt1+1),(nt1+1)*dt_xiNg1,&
	    qxiNg1(1,j-2,nt1+2),t,zt(Nx+3,j) )
	  call interpj(nt1*dt_xing1,qxiNg1(2,j-2,nt1+1),(nt1+1)*dt_xiNg1,&
	    qxiNg1(2,j-2,nt1+2),t,qt(1,Nx+3,j))
	  call interpj(nt1*dt_xiNg1,qxiNg1(3,j-2,nt1+1),(nt1+1)*dt_xiNg1,&
	    qxiNg1(3,j-2,nt1+2),t,qt(2,Nx+3,j))
	  call interpj(nt1*dt_xiNg1,qxiNg1(4,j-2,nt1+1),(nt1+1)*dt_xiNg1,&
	    qxiNg1(4,j-2,nt1+2),t,qt(3,Nx+3,j))
	end do
      end if	  
      if (optxiNg2.eq.0) then
	do j=3,Ny+2
	  zt(Nx+4,j)=qxiNg2(1,j-2,nt1+1)
	  qt(1,Nx+4,j)=qxiNg2(2,j-2,nt1+1)
	  qt(2,Nx+4,j)=qxiNg2(3,j-2,nt1+1)  
	  qt(3,Nx+4,j)=qxiNg2(4,j-2,nt1+1)
	end do
      else
	do j=3,Ny+2
	  call interpj(nt2*dt_xiNg2,qxiNg2(1,j-2,nt2+1),(nt2+1)*dt_xiNg2,&
	    qxiNg2(1,j-2,nt2+2),t,zt(Nx+4,j) )
	  call interpj(nt2*dt_xiNg2,qxiNg2(2,j-2,nt2+1),(nt2+1)*dt_xiNg2,&
	    qxiNg2(2,j-2,nt2+2),t,qt(1,Nx+4,j))
	  call interpj(nt2*dt_xiNg2,qxiNg2(3,j-2,nt2+1),(nt2+1)*dt_xiNg2,&
	    qxiNg2(3,j-2,nt2+2),t,qt(2,Nx+4,j))
	  call interpj(nt2*dt_xiNg2,qxiNg2(4,j-2,nt2+1),(nt2+1)*dt_xiNg2,&
	    qxiNg2(4,j-2,nt2+2),t,qt(3,Nx+4,j))
	end do
      end if

    CASE(1) 
    !Solid Wall 
    !Ojo que son las velocidades contravariantes las que se usan
    !pero si las metricas se asumen iguales entonces queda lo mismo que en cartesianas.    
      DO j=3,Ny+2
	qt(1,Nx+4,j)=qn(1,Nx-1,j-2)	!hNx+4=hNx-1
	qt(1,Nx+3,j)=qn(1,Nx,j-2)	!hNx+3=hNx
	
	qt(2,Nx+4,j)=-qn(2,Nx-1,j-2)	!uNx+4=-uNx-1
	qt(2,Nx+3,j)=-qn(2,Nx,j-2)	!uNx+3=-uNx
	
	qt(3,Nx+4,j)=-qn(3,Nx-1,j-2)	!vNx+4=-vNx-1
	qt(3,Nx+3,j)=-qn(3,Nx,j-2)	!vNx+3=-vNx	
      END DO
		
    CASE(2) 
      !Periodic
      !handled by communication (exchange_2d)
	
    CASE(3) 
      !Free Exit	
      DO j=3,Ny+2	
	qt(1,Nx+4,j)=qn(1,Nx-1,j-2)	!hNx+4=hNx-1
	qt(1,Nx+3,j)=qn(1,Nx,j-2)	!hNx+3=hNx
	
	qt(2,Nx+4,j)=qn(2,Nx-1,j-2)	!uNx+4=uNx-1
	qt(2,Nx+3,j)=qn(2,Nx,j-2)	!uNx+3=uNx
	
	qt(3,Nx+4,j)=qn(3,Nx-1,j-2)	!vNx+4=vNx-1
	qt(3,Nx+3,j)=qn(3,Nx,j-2)	!vNx+3=vNx	
      END DO		
      
    CASE(4) 
      !HAY QUE REVISAR ESTAS FUNCIONES		
      IF (pasoRK==1.OR.pasoRK==3) then
	IF (GA2==1) THEN
	  call genabsNxi_1_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,etaR2,Nsenal2,&
	    h02,t,dt,qn,zt,xi,qA2,zA2)
	ELSE IF (GA2==2) THEN
	  call genabsNxi_2_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,qs2,hs2,h02,Nsenal2,&
	    t,dt,qn,zt,xi,qA2,zA2)
	ELSE IF (GA2==3) THEN
	  call genabsNxi_3_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,us2,hs2,h02,Nsenal2,&
	    t,dt,qn,zt,xi,qA2,zA2)   
	END IF		
	qA20=qA2    		
      ELSE IF (pasoRK==2.OR.pasoRK==4) then
	IF (GA2==1) THEN
	  call genabsNxi_1_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,etaR2,Nsenal2,&
	    h02,t,dt,qn,zt,xi,qA20,qA2,zA2)
	ELSE IF (GA2==2) THEN
	  call genabsNxi_2_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,qs2,hs2,h02,Nsenal2,&
	    t,dt,qn,zt,xi,qA20,qA2,zA2)
	ELSE IF (GA2==3) THEN
	  call genabsNxi_3_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,us2,hs2,h02,Nsenal2,&
	    t,dt,qn,zt,xi,qA20,qA2,zA2)
	END IF
      END IF
	  
      DO j=3,Ny+2		
	qt(1,Nx+4,j)=qA2(1,j-2)		!h-1=h2
	qt(1,Nx+3,j)=qA2(1,j-2)		!h0=h1
	  
	qt(2,Nx+4,j)=qA2(2,j-2)		!u-1=-u2
	qt(2,Nx+3,j)=qA2(2,j-2)		!u0=-u1
	  
	qt(3,Nx+4,j)=qA2(3,j-2)		!v-1=-v2
	qt(3,Nx+3,j)=qA2(3,j-2)		!v0=-v1	
      END DO
		
  END SELECT
!---------------------------------------------------------------------------

  SELECT CASE (CB(3)) 
    CASE(0) !Customized boundary condition
      !1 Interpolate
      nt1=min(int(t/dt_eta0g1),nt_eta0g1-1)
      nt2=min(int(t/dt_eta0g2),nt_eta0g2-1)
      if (opteta0g1.eq.0) then
	do i=3,Nx+2
	  zt(i,1)=qeta0g1(1,i-2,nt1+1)
	  qt(1,i,1)=qeta0g1(2,i-2,nt1+1)
	  qt(2,i,1)=qeta0g1(3,i-2,nt1+1)	    
	  qt(3,i,1)=qeta0g1(4,i-2,nt1+1) 
	end do
      else
	do i=3,Nx+2
	  call interpj(nt1*dt_eta0g1,qeta0g1(1,i-2,nt2+1),(nt1+1)*dt_eta0g1,qeta0g1(1,i-2,nt1+2),t,zt(i,1) )
	  call interpj(nt1*dt_eta0g1,qeta0g1(2,i-2,nt2+1),(nt1+1)*dt_eta0g1,qeta0g1(2,i-2,nt1+2),t,qt(1,i,1))
	  call interpj(nt1*dt_eta0g1,qeta0g1(3,i-2,nt2+1),(nt1+1)*dt_eta0g1,qeta0g1(3,i-2,nt1+2),t,qt(2,i,1))
	  call interpj(nt1*dt_eta0g1,qeta0g1(4,i-2,nt2+1),(nt1+1)*dt_eta0g1,qeta0g1(4,i-2,nt1+2),t,qt(3,i,1))
	end do
      end if
	  
      if (opteta0g2.eq.0) then
	do i=3,Nx+2
	  zt(i,2)=qeta0g2(1,i-2,nt2+1)
	  qt(1,i,2)=qeta0g2(2,i-2,nt2+1)
	  qt(2,i,2)=qeta0g2(3,i-2,nt2+1)
	  qt(3,i,2)=qeta0g2(4,i-2,nt2+1)	
	end do
      else
	do i=3,Nx+2
	  call interpj(nt2*dt_eta0g2,qeta0g2(1,i-2,nt2+1),(nt2+1)*dt_eta0g2,qeta0g2(1,i-2,nt2+2),t,zt(i,2) )
	  call interpj(nt2*dt_eta0g2,qeta0g2(2,i-2,nt2+1),(nt2+1)*dt_eta0g2,qeta0g2(2,i-2,nt2+2),t,qt(1,i,2))
	  call interpj(nt2*dt_eta0g2,qeta0g2(3,i-2,nt2+1),(nt2+1)*dt_eta0g2,qeta0g2(3,i-2,nt2+2),t,qt(2,i,2))
	  call interpj(nt2*dt_eta0g2,qeta0g2(4,i-2,nt2+1),(nt2+1)*dt_eta0g2,qeta0g2(4,i-2,nt2+2),t,qt(3,i,2))
	end do
      end if

    CASE(1) 
    !Solid Wall 
    !Ojo que son las velocidades contravariantes las que se usan
    !pero si las metricas se asumen iguales entonces queda lo mismo que en cartesianas.
      DO i=3,Nx+2	
	qt(1,i,1)=qn(1,i-2,2)		!h-1=h2
	qt(1,i,2)=qn(1,i-2,1)		!h0=h1
	
	qt(2,i,1)=-qn(2,i-2,2)		!u-1=-u2
	qt(2,i,2)=-qn(2,i-2,1)		!u0=-u1
	
	qt(3,i,1)=-qn(3,i-2,2)		!v-1=-v2
	qt(3,i,2)=-qn(3,i-2,1)		!v0=-v1	
      END DO
	
    CASE(2) !Periodic	
      !exchange_2d
    CASE(3) !Free Exit	
      DO i=3,Nx+2
	qt(1,i,1)=qn(1,i-2,2)		!h-1=h2
	qt(1,i,2)=qn(1,i-2,1)		!h0=h1
	
	qt(2,i,1)=qn(2,i-2,2)		!u-1=u2
	qt(2,i,2)=qn(2,i-2,1)		!u0=u1
	
	qt(3,i,1)=qn(3,i-2,2)		!v-1=v2
	qt(3,i,2)=qn(3,i-2,1)		!v0=v1	
      END DO

    CASE(4)
! 	IF (GA3==1) THEN
! 	call genabs0eta_1(Nx,Ny,Fr2,deta,etaL3,Nsenal3,h03,t,dt,qn,zt,eta,qA3,zA3)
!	ELSE IF (GA3==2) THEN
! 	call genabs0eta_2(Nx,Ny,Fr2,deta,qs3,hs3,Nsenal3,t,dt,qn,zt,eta,qA3,zA3)
! 	ELSE IF (GA3==3) THEN   
! 	call genabs0eta_3(Nx,Ny,Fr2,deta,qs3,us3,Nsenal3,t,dt,qn,zt,eta,qA3,zA3)
! 	END IF
! 	DO i=3,Nx+2
! 	
! 	qt(1,i,1)=qA3(1,i-2)		!h-1=h2
! 	qt(1,i,2)=qA3(1,i-2)		!h0=h1
! 	
! 	qt(2,i,1)=qA3(2,i-2)		!u-1=u2
! 	qt(2,i,2)=qA3(2,i-2)		!u0=-u1
! 		
! 	qt(3,i,1)=qA3(3,i-2)		!v-1=v2
! 	qt(3,i,2)=qA3(3,i-2)		!v0=-v1
! 	
! 	END DO
END SELECT

!----------------------------------------------------------------------------------

  
  SELECT CASE (CB(4))
    case(0) !Customized boundary condition
      !1 Interpolate
      nt1=min(int(t/dt_etaNg1),nt_etaNg1-1)
      nt2=min(int(t/dt_etaNg2),nt_etaNg2-1)
      if (optetaNg1.eq.0) then
	do i=3,Nx+2
	  zt(i,Ny+3)=qetaNg1(1,i-2,nt1+1)
	  qt(1,i,Ny+3)=qetaNg1(2,i-2,nt1+1)
	  qt(2,i,Ny+3)=qetaNg1(3,i-2,nt1+1)	    
	  qt(3,i,Ny+3)=qetaNg1(4,i-2,nt1+1) 
	end do
      else
	do i=3,Nx+2
	  call interpj(nt1*dt_etaNg1,qetaNg1(1,i-2,nt1+1),(nt1+1)*dt_etaNg1,qetaNg1(1,i-2,nt1+2),t,zt(i,Ny+3) )
	  call interpj(nt1*dt_etaNg1,qetaNg1(2,i-2,nt1+1),(nt1+1)*dt_etaNg1,qetaNg1(2,i-2,nt1+2),t,qt(1,i,Ny+3))
	  call interpj(nt1*dt_etaNg1,qetaNg1(3,i-2,nt1+1),(nt1+1)*dt_etaNg1,qetaNg1(3,i-2,nt1+2),t,qt(2,i,Ny+3))
	  call interpj(nt1*dt_etaNg1,qetaNg1(4,i-2,nt1+1),(nt1+1)*dt_etaNg1,qetaNg1(4,i-2,nt1+2),t,qt(3,i,Ny+3))
	end do
      end if
      
      if (optetaNg2.eq.0) then
	do i=3,Nx+2
	  zt(i,Ny+4)=qetaNg2(1,i-2,nt2+1)
	  qt(1,i,Ny+4)=qetaNg2(2,i-2,nt2+1)
	  qt(2,i,Ny+4)=qetaNg2(3,i-2,nt2+1)
	  qt(3,i,Ny+4)=qetaNg2(4,i-2,nt2+1)	  
	end do
      else
	do i=3,Nx+2  
	  call interpj(nt2*dt_etaNg2,qetaNg2(1,i-2,nt2+1),(nt2+1)*dt_etaNg2,qetaNg2(1,i-2,nt2+2),t,zt(i,Ny+4) )
	  call interpj(nt2*dt_etaNg2,qetaNg2(2,i-2,nt2+1),(nt2+1)*dt_etaNg2,qetaNg2(2,i-2,nt2+2),t,qt(1,i,Ny+4))
	  call interpj(nt2*dt_etaNg2,qetaNg2(3,i-2,nt2+1),(nt2+1)*dt_etaNg2,qetaNg2(3,i-2,nt2+2),t,qt(2,i,Ny+4))
	  call interpj(nt2*dt_etaNg2,qetaNg2(4,i-2,nt2+1),(nt2+1)*dt_etaNg2,qetaNg2(4,i-2,nt2+2),t,qt(3,i,Ny+4))
	end do
      end if
    CASE(1)	
      DO i=3,Nx+2
	qt(1,i,Ny+4)=qn(1,i-2,Ny-1)	!hNy+4=hNy-1
	qt(1,i,Ny+3)=qn(1,i-2,Ny)	!hNy+3=hNy
	
	qt(2,i,Ny+4)=-qn(2,i-2,Ny-1)	!uNy+4=-uNy-1
	qt(2,i,Ny+3)=-qn(2,i-2,Ny)	!uNy+3=-uNy
	
	qt(3,i,Ny+4)=-qn(3,i-2,Ny-1)	!vNy+4=-vNy-1
	qt(3,i,Ny+3)=-qn(3,i-2,Ny)	!vNy+3=-vNy	
      END DO
    CASE(2) !Periodic
      DO i=3,Nx+2	
	qt(1,i,Ny+4)=qn(1,i-2,2)	!hNy+4=h2
	qt(1,i,Ny+3)=qn(1,i-2,1)	!hNy+3=h1
	
	qt(2,i,Ny+4)=qn(2,i-2,2)	!uNy+4=u2
	qt(2,i,Ny+3)=qn(2,i-2,1)	!uNy+3=u1
	
	qt(3,i,Ny+4)=qn(3,i-2,2)	!vNy+4=v2
	qt(3,i,Ny+3)=qn(3,i-2,1)	!vNy+3=v1	
      END DO	
    
    CASE(3) !Free Exit	
      DO i=3,Nx+2	
	qt(1,i,Ny+4)=qn(1,i-2,Ny-1)	!hNy+4=hNy-1
	qt(1,i,Ny+3)=qn(1,i-2,Ny)	!hNy+3=hNy
	
	qt(2,i,Ny+4)=qn(2,i-2,Ny-1)	!uNy+4=-uNy-1
	qt(2,i,Ny+3)=qn(2,i-2,Ny)	!uNy+3=-uNy
	
	qt(3,i,Ny+4)=qn(3,i-2,Ny-1)	!vNy+4=-vNy-1
	qt(3,i,Ny+3)=qn(3,i-2,Ny)	!vNy+3=-vNy
      END DO
    CASE(4) !No está programada aun
! 	IF (GA4==1) THEN
! 	call genabsNeta_1(Nx,Ny,Fr2,deta,etaL3,Nsenal3,h03,t,dt,qn,zt,eta,qA4,zA4)
! 	ELSE IF (GA4==2) THEN
! 	call genabsNeta_2(Nx,Ny,Fr2,deta,qs3,hs3,Nsenal3,t,dt,qn,zt,eta,qA4,zA4)
! 	ELSE IF (GA4==3) THEN   
! 	call genabsNeta_3(Nx,Ny,Fr2,deta,qs3,us3,Nsenal3,t,dt,qn,zt,eta,qA4,zA4)
! 	END IF!
! 	DO i=3,Nx+2
! 	
! 	qt(1,i,Ny+4)=qA4(1,i-2)		!h-1=h2
! 	qt(1,i,Ny+3)=qA4(1,i-2)		!h0=h1
! 		
! 	qt(2,i,Ny+4)=qA4(2,i-2)		!u-1=u2
! 	qt(2,i,Ny+3)=qA4(2,i-2)		!u0=-u1
! 	
! 	qt(3,i,Ny+4)=qA4(3,i-2)		!v-1=v2
! 	qt(3,i,Ny+3)=qA4(3,i-2)		!v0=-v1
! 		
! 	END DO
	
END SELECT
end subroutine bcs

subroutine interpj(x1,y1,x2,y2,x,y)
  !interpolates y(x) using the line (x1,y1)-(x2,y2) if x is between x1,x2
  !extrapolates order 0 if x is outside [x1,x2]  
  real (kind=8) ::x1,y1,x2,y2,x,y
  if (x.lt.x1) then
    y=y1
  elseif (x2.lt.x) then
    y=y2
  else
    y=y1*(x2-x)/(x2-x1)+y2*(x1-x)/(x1-x2)
  end if
end subroutine