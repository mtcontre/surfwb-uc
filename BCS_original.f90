SUBROUTINE bcs(fopt,Cf,Coef,pasoRK,t,dt,Fr2,caso,qn,z,Nx,Ny,CB,xi,eta,Jac,dxi,deta,qt,xit,etat,zt)

!2nd Order Boundary Conditions!
!bcs(state variables at n time,Nbx,Nby,BoundaryConditions, xi_metrics, eta_metrics, 
!    x_xi,y_xi,x_eta,y_eta,matrix with ghost cells)

!USE global_variables
!USE geometries
USE senales
USE coords
implicit none
integer :: Nx,Ny, caso
real (kind=8), dimension(3,Nx,Ny)	::qn
real (kind=8), dimension(2,Nx,Ny)	::xi,eta
real (kind=8), dimension(3,Nx+4,Ny+4)	::qt
real (kind=8), dimension(Nx+4,Ny+4)	::zt
real (kind=8), dimension(Nx,Ny)		::z, Coef, Jac
real (kind=8), dimension(2,Nx+4,Ny+4)	::xit,etat
integer,dimension(4) :: CB  
real (kind=8)	:: t,dt,Fr2,dxi,deta,Zx1,Zxi1, Zeta1, Zx0, Zxi0, Zeta0, ZxiN, ZetaN, ZxN, ZxN1, ZxiN1, ZetaN1, Zy1, Zy0, ZyN, ZyN1, hF
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


!0.Initialization of matrices
DO i=1,Nx+4; DO j=1,Ny+4
qt(:,i,j)=0.0D0
zt(i,j)=0.0D0
xit(1,i,j)=0.0D0
xit(2,i,j)=0.0D0
etat(1,i,j)=0.0D0
etat(2,i,j)=0.0D0
END DO; END DO

!1. Complete Inside Cells
DO i=3,Nx+2; DO j=3,Ny+2
qt(:,i,j)=qn(:,i-2,j-2)
xit(:,i,j)=xi(:,i-2,j-2)
etat(:,i,j)=eta(:,i-2,j-2)
zt(i,j)=z(i-2,j-2)
END DO; END DO

!2. Complete Ghost Cells

!Xi=1, CB(1)


	!Metricas
	DO j=3,Ny+2
	xit(:,1,j)=xi(:,2,j-2)		!xi-1=xi2
	xit(:,2,j)=xi(:,1,j-2)		!xi0=xi12
	
	etat(:,1,j)=eta(:,2,j-2)	!eta-1=eta2
	etat(:,2,j)=eta(:,1,j-2)	!eta0=eta1
	
	
	!!Batimetria Ghost Cells
	!Celda 0, zt(2,j)
! 	Zxi1=(zt(4,j)-zt(3,j))/dxi
! 	Zeta1=(zt(4,j)-zt(3,j))/deta
! 	Zx1=Zxi1*xit(1,3,j)+Zeta1*etat(1,3,j)
! 	zt(2,j)=zt(3,j)-Zx1/(xit(1,2,j)/dxi+etat(1,2,j)/deta)
! 
! 	!Celda -1, zt(1,j)
! 	Zxi0=(zt(3,j)-zt(2,j))/dxi
! 	Zeta0=(zt(3,j)-zt(2,j))/deta
! 	Zx0=Zxi0*xit(1,2,j)+Zeta0*etat(1,2,j)
! 	zt(1,j)=zt(2,j)-Zx0/(xit(1,1,j)/dxi+etat(1,1,j)/deta)
	
	zt(1,j)=zt(3,j)
	zt(2,j)=zt(3,j)
	END DO


SELECT CASE (CB(1))

	CASE(1) !Solid Wall !Ojo que son las velocidades contravariantes las que se usan, pero si las metricas se asumen iguales entonces queda lo mismo que en cartesianas.
	
	DO j=3,Ny+2
	
	qt(1,1,j)=qn(1,2,j-2)		!h-1=h2
	qt(1,2,j)=qn(1,1,j-2)		!h0=h1
	
	qt(2,1,j)=-qn(2,2,j-2)		!u-1=-u2
	qt(2,2,j)=-qn(2,1,j-2)		!u0=-u1
	
	qt(3,1,j)=-qn(3,2,j-2)		!v-1=-v2
	qt(3,2,j)=-qn(3,1,j-2)		!v0=-v1
	
	END DO
		
	CASE(2) !Periodic
	
	if (caso==8) then
	
	DO j=3,Ny+2
	
	qt(1,1,j)=qn(1,Nx-2,j-2)	!h-1=hNx-1
	qt(1,2,j)=qn(1,Nx-1,j-2)		!h0=hNx
	
	qt(2,1,j)=qn(2,Nx-2,j-2)	!u-1=uNx-1
	qt(2,2,j)=qn(2,Nx-1,j-2)		!u0=uNx
	
	qt(3,1,j)=qn(3,Nx-2,j-2)	!v-1=vNx-1
	qt(3,2,j)=qn(3,Nx-1,j-2)		!v0=vNx
	END DO
	
	else
	
	DO j=3,Ny+2
	
	qt(1,1,j)=qn(1,Nx-1,j-2)	!h-1=hNx-1
	qt(1,2,j)=qn(1,Nx,j-2)		!h0=hNx
	
	qt(2,1,j)=qn(2,Nx-1,j-2)	!u-1=uNx-1
	qt(2,2,j)=qn(2,Nx,j-2)		!u0=uNx
	
	qt(3,1,j)=qn(3,Nx-1,j-2)	!v-1=vNx-1
	qt(3,2,j)=qn(3,Nx,j-2)		!v0=vNx
	END DO
	
	end if
	
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
			call genabs0xi_1_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,etaL1,Nsenal1,h01,t,dt,qn,zt,xi,qA1,zA1)
			ELSE IF (GA1==2) THEN
			call genabs0xi_2_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,qs1,hs1,h01,Nsenal1,t,dt,qn,zt,xi,qA1,zA1)
			
			ELSE IF (GA1==3) THEN
			call genabs0xi_3_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,us1,hs1,h01,Nsenal1,t,dt,qn,zt,xi,qA1,zA1)
			
! 			ELSE IF (GA1==4) THEN
! 			call genabs0xi_4_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,hs1,h01,Nsenal1,t,dt,qn,zt,xi,qA1,zA1)
 			
			! GENABS Con diferentes señales en cada nodo
			ELSE IF (GA1==9) THEN
			
			call genabs0xi_9_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,etaL9,timeS9,Nsenal1,h01,t,dt,qn,zt,xi,qA1,zA1)
			
			END IF
			
		END IF	
 		qA10=qA1

		
		IF (pasoRK==2.OR.pasoRK==4) then

			IF (GA1==1) THEN
			call genabs0xi_1_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,etaL1,Nsenal1,h01,t,dt,qn,zt,xi,qA10,qA1,zA1)
			
			ELSE IF (GA1==2) THEN
			call genabs0xi_2_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,qs1,hs1,h01,Nsenal1,t,dt,qn,zt,xi,qA10,qA1,zA1)
			
			ELSE IF (GA1==3) THEN
			call genabs0xi_3_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,us1,hs1,h01,Nsenal1,t,dt,qn,zt,xi,qA10,qA1,zA1)
			
! 			ELSE IF (GA1==4) THEN
! 			call genabs0xi_4_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,hs1,h01,Nsenal1,t,dt,qn,zt,xi,qA10,qA1,zA1)
			ELSE IF (GA1==9) THEN
			call genabs0xi_9_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,etaL9,timeS9,Nsenal1,h01,t,dt,qn,zt,xi,qA10,qA1,zA1)
			
			END IF
		
		END IF


		DO j=3,Ny+2

		qt(1,1,j)=qA1(1,j-2)		!h-1=h2
		qt(1,2,j)=qA1(1,j-2)		!h0=h1
		
		qt(2,1,j)=qA1(2,j-2)		!u-1=-u2
		qt(2,2,j)=qA1(2,j-2)		!u0=-u1
		
		qt(3,1,j)=qA1(3,j-2)		!v-1=-v2
		qt(3,2,j)=qA1(3,j-2)		!v0=-v1
	
		END DO	


	 CASE(5) 
	 		
		if (IO1==2) then !Outflow, !Cota fija aguas afuera 
		 !OUTFLOW0_xi(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,hs1,Nsenal1,t,dt,qn,zt,xi,qA1,zA1)
		!OUTFLOW0_xi(pasoRK,fopt,tipo,MC,Nx,Ny,Fr2,dep,etas1,timeS,Ns,t,dt,q,zt,ep_x,ep2_x,qA,zA)
		call OUTFLOW0_xi(pasoRK,fopt,Cf,Coef,Nx,Ny,Fr2,dxi,etas1,timeS1,Nsenal1,t,dt,qn,zt,xi,eta,qA1,zA1)
		else
		call INFLOW0_xi(pasoRK,fopt,Cf,Coef,Nx,Ny,Fr2,dxi,qsx1,qsy1,etas1,timeS1,Nsenal1,t,dt,qn,zt,xi,eta,qA1,zA1)
		end if
		
		
		DO j=3,Ny+2
		
		qt(1,1,j)=qA1(1,j-2)		!h-1=h2
		qt(1,2,j)=qA1(1,j-2)		!h0=h1
		
		qt(2,1,j)=qA1(2,j-2)		!u-1=-u2
		qt(2,2,j)=qA1(2,j-2)		!u0=-u1
		
		qt(3,1,j)=qA1(3,j-2)		!v-1=-v2
		qt(3,2,j)=qA1(3,j-2)		!v0=-v1
	
		END DO
	
END SELECT
!---------------------------------------------------------------------------
!Xi=Nx, CB(2)
	DO j=3,Ny+2
	
	xit(:,Nx+4,j)=xi(:,Nx-1,j-2)	!xiNx+4=xiNx-1
	xit(:,Nx+3,j)=xi(:,Nx,j-2)	!xiNx+3=xiNx
	
	etat(:,Nx+4,j)=eta(:,Nx-1,j-2)	!etaNx+4=etaNx-1
	etat(:,Nx+3,j)=eta(:,Nx,j-2)	!etaNx+3=etaNx

	
	!Batimetria Celdas Fantasma
! 	!Celda N+1, zt(N+3,j)
! 	ZxiN=(zt(Nx+2,j)-zt(Nx+1,j))/dxi
! 	ZetaN=(zt(Nx+2,j)-zt(Nx+1,j))/deta
! 	ZxN=ZxiN*xit(1,Nx+2,j)+ZetaN*etat(1,Nx+2,j)
! 	zt(Nx+3,j)=zt(Nx+2,j)+ZxN/(xit(1,Nx+3,j)/dxi+etat(1,Nx+3,j)/deta)
! 
! 	!Celda N+2, zt(N+4,j)
! 	ZxiN1=(zt(Nx+3,j)-zt(Nx+2,j))/dxi
! 	ZetaN1=(zt(Nx+3,j)-zt(Nx+2,j))/deta
! 	ZxN1=ZxiN1*xit(1,Nx+3,j)+ZetaN1*etat(1,Nx+3,j)
! 	zt(Nx+4,j)=zt(Nx+3,j)+ZxN/(xit(1,Nx+4,j)/dxi+etat(1,Nx+4,j)/deta)
	
	zt(Nx+4,j)=zt(Nx+2,j)
	zt(Nx+3,j)=zt(Nx+2,j)
	
	END DO
	
SELECT CASE (CB(2))


	
	CASE(1) !Solid Wall !Ojo que son las velocidades contravariantes las que se usan, pero si las metricas se asumen iguales entonces queda lo mismo que en cartesianas.
	
	DO j=3,Ny+2
	

	qt(1,Nx+4,j)=qn(1,Nx-1,j-2)	!hNx+4=hNx-1
	qt(1,Nx+3,j)=qn(1,Nx,j-2)	!hNx+3=hNx
	
	qt(2,Nx+4,j)=-qn(2,Nx-1,j-2)	!uNx+4=-uNx-1
	qt(2,Nx+3,j)=-qn(2,Nx,j-2)	!uNx+3=-uNx
	
	qt(3,Nx+4,j)=-qn(3,Nx-1,j-2)	!vNx+4=-vNx-1
	qt(3,Nx+3,j)=-qn(3,Nx,j-2)	!vNx+3=-vNx
	
	END DO
		
	CASE(2) !Periodic
	if (caso==8) then
	
	DO j=3,Ny+2
	xit(:,Nx+4,j)=xi(:,3,j-2)	!xiNx+4=xi2
	xit(:,Nx+3,j)=xi(:,2,j-2)	!xiNx+3=xi1
	
	etat(:,Nx+4,j)=eta(:,3,j-2)	!etaNx+4=eta2
	etat(:,Nx+3,j)=eta(:,2,j-2)	!etaNx+3=eta1
		
	qt(1,Nx+4,j)=qn(1,3,j-2)	!hNx+4=h2
	qt(1,Nx+3,j)=qn(1,2,j-2)	!hNx+3=h1
	
	qt(2,Nx+4,j)=qn(2,3,j-2)	!uNx+4=u2
	qt(2,Nx+3,j)=qn(2,2,j-2)	!uNx+3=u1
	
	qt(3,Nx+4,j)=qn(3,3,j-2)	!vNx+4=v2
	qt(3,Nx+3,j)=qn(3,2,j-2)	!vNx+3=v1
	END DO
	
	else
	
	DO j=3,Ny+2
	
	qt(1,Nx+4,j)=qn(1,2,j-2)	!hNx+4=h2
	qt(1,Nx+3,j)=qn(1,1,j-2)	!hNx+3=h1
	
	qt(2,Nx+4,j)=qn(2,2,j-2)	!uNx+4=u2
	qt(2,Nx+3,j)=qn(2,1,j-2)	!uNx+3=u1
	
	qt(3,Nx+4,j)=qn(3,2,j-2)	!vNx+4=v2
	qt(3,Nx+3,j)=qn(3,1,j-2)	!vNx+3=v1
	
	END DO
	end if
	
	CASE(3) !Free Exit
	
	DO j=3,Ny+2
	
	qt(1,Nx+4,j)=qn(1,Nx-1,j-2)	!hNx+4=hNx-1
	qt(1,Nx+3,j)=qn(1,Nx,j-2)	!hNx+3=hNx
	
	qt(2,Nx+4,j)=qn(2,Nx-1,j-2)	!uNx+4=uNx-1
	qt(2,Nx+3,j)=qn(2,Nx,j-2)	!uNx+3=uNx
	
	qt(3,Nx+4,j)=qn(3,Nx-1,j-2)	!vNx+4=vNx-1
	qt(3,Nx+3,j)=qn(3,Nx,j-2)	!vNx+3=vNx
	
	END DO	
	
	
	CASE(4) !HAY QUE REVISAR ESTAS FUNCIONES
		IF (pasoRK==1.OR.pasoRK==3) then
		    IF (GA2==1) THEN
		    call genabsNxi_1_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,etaR2,Nsenal2,h02,t,dt,qn,zt,xi,qA2,zA2)
		    
		    ELSE IF (GA2==2) THEN
		    call genabsNxi_2_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,qs2,hs2,h02,Nsenal2,t,dt,qn,zt,xi,qA2,zA2)
		    
		    ELSE IF (GA2==3) THEN
		    call genabsNxi_3_1(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,us2,hs2,h02,Nsenal2,t,dt,qn,zt,xi,qA2,zA2)
		    
	    
		    	    
		    END IF
		
		qA20=qA2    
		
		ELSE IF (pasoRK==2.OR.pasoRK==4) then
		    IF (GA2==1) THEN
		    call genabsNxi_1_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,etaR2,Nsenal2,h02,t,dt,qn,zt,xi,qA20,qA2,zA2)
		    ELSE IF (GA2==2) THEN
		    call genabsNxi_2_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,qs2,hs2,h02,Nsenal2,t,dt,qn,zt,xi,qA20,qA2,zA2)
		    
		    ELSE IF (GA2==3) THEN
		    call genabsNxi_3_2(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,us2,hs2,h02,Nsenal2,t,dt,qn,zt,xi,qA20,qA2,zA2)
		        
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
	 
	 CASE(5) 
	 		
		if (IO2==2) then !Outflow, !Cota fija aguas afuera 
		call OUTFLOWN_xi(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,hs2,Nsenal2,t,dt,qn,zt,xi,qA2,zA2)
		
		else if (IO2==1) then
		call INFLOWN_xi(pasoRK,fopt,Cf,Coef,Nx,Ny,Fr2,dxi,qsx2,qsy2,etas2,timeS2,Nsenal2,t,dt,qn,zt,xi,eta,qA2,zA2)
		
		
		else
		
		call OUTFLOWN_V(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,us2,Nsenal2,t,dt,qn,zt,xi,qA2,zA2)
		
		end if
		
		
		
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
!Eta=1, CB(3)

	DO i=3,Nx+2
	
	!Metricas
	xit(:,i,1)=xi(:,i-2,2)		!xi-1=xi2
	xit(:,i,2)=xi(:,i-2,1)		!xi0=xi12
	
	etat(:,i,1)=eta(:,i-2,2)	!eta-1=eta2
	etat(:,i,2)=eta(:,i-2,1)	!eta0=eta1
	
	!Batimetria celdas ficticias
! 	!Celda 0, zt(i,2)
! 	Zxi1=(zt(i,4)-zt(i,3))/dxi
! 	Zeta1=(zt(i,4)-zt(i,3))/deta
! 	Zy1=Zxi1*xit(2,i,3)+Zeta1*etat(2,i,3)
! 	zt(i,2)=zt(i,3)-Zy1/(xit(2,i,2)/dxi+etat(2,i,2)/deta)
! 
! 	!Celda -1, zt(i,1)
! 	Zxi0=(zt(i,3)-zt(i,2))/dxi
! 	Zeta0=(zt(i,3)-zt(i,2))/deta
! 	Zy0=Zxi0*xit(2,i,2)+Zeta0*etat(2,i,2)
! 	zt(i,1)=zt(i,2)-Zy0/(xit(2,i,1)/dxi+etat(2,i,1)/deta)
	zt(i,1)=zt(i,3)
	zt(i,2)=zt(i,3)
	
	END DO
	
SELECT CASE (CB(3)) 

	CASE(1) !Solid Wall !Ojo que son las velocidades contravariantes las que se usan, pero si las metricas se asumen iguales entonces queda lo mismo que en cartesianas.
	
	DO i=3,Nx+2
	
	qt(1,i,1)=qn(1,i-2,2)		!h-1=h2
	qt(1,i,2)=qn(1,i-2,1)		!h0=h1
	
	qt(2,i,1)=-qn(2,i-2,2)		!u-1=-u2
	qt(2,i,2)=-qn(2,i-2,1)		!u0=-u1
	
	qt(3,i,1)=-qn(3,i-2,2)		!v-1=-v2
	qt(3,i,2)=-qn(3,i-2,1)		!v0=-v1
	
	END DO
	
	CASE(2) !Periodic
	
	DO i=3,Nx+2
	
	qt(1,i,1)=qn(1,i-2,Ny-1)	!h-1=hNy-1
	qt(1,i,2)=qn(1,i-2,Ny)		!h0=hNy
	
	qt(2,i,1)=qn(2,i-2,Ny-1)	!u-1=uNy-1
	qt(2,i,2)=qn(2,i-2,Ny)		!u0=uNy
	
	qt(3,i,1)=qn(3,i-2,Ny-1)	!v-1=vNy-1
	qt(3,i,2)=qn(3,i-2,Ny)		!v0=vNy
	
	END DO	
	

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
! 		IF (GA3==1) THEN
! 		call genabs0eta_1(Nx,Ny,Fr2,deta,etaL3,Nsenal3,h03,t,dt,qn,zt,eta,qA3,zA3)
! 		ELSE IF (GA3==2) THEN
! 		call genabs0eta_2(Nx,Ny,Fr2,deta,qs3,hs3,Nsenal3,t,dt,qn,zt,eta,qA3,zA3)
! 		ELSE IF (GA3==3) THEN   
! 		call genabs0eta_3(Nx,Ny,Fr2,deta,qs3,us3,Nsenal3,t,dt,qn,zt,eta,qA3,zA3)
! 		END IF
! 		DO i=3,Nx+2
! 	
! 		qt(1,i,1)=qA3(1,i-2)		!h-1=h2
! 		qt(1,i,2)=qA3(1,i-2)		!h0=h1
! 		
! 		qt(2,i,1)=qA3(2,i-2)		!u-1=u2
! 		qt(2,i,2)=qA3(2,i-2)		!u0=-u1
! 		
! 		qt(3,i,1)=qA3(3,i-2)		!v-1=v2
! 		qt(3,i,2)=qA3(3,i-2)		!v0=-v1
! 		
! 		END DO
	CASE(5) 
	 		
		if (IO3==2) then !Outflow, !Cota fija aguas afuera 
		!call OUTFLOWN(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,hs2,Nsenal2,t,dt,qn,zt,xi,qA2,zA2)
		!ERROR aun no está programada
		
		else
		
		call INFLOW0_eta(pasoRK,fopt,Cf,Coef,Nx,Ny,Fr2,deta,qsx3,qsy3,etas3,timeS3,Nsenal3,t,dt,qn,zt,eta,xi,qA3,zA3)
		
		end if
		
	
		DO i=3,Nx+2

		qt(1,i,1)=qA3(1,i-2)		!h-1=h2
		qt(1,i,2)=qA3(1,i-2)		!h0=h1
		
		qt(2,i,1)=qA3(2,i-2)		!u-1=-u2
		qt(2,i,2)=qA3(2,i-2)		!u0=-u1
		
		qt(3,i,1)=qA3(3,i-2)		!v-1=-v2
		qt(3,i,2)=qA3(3,i-2)		!v0=-v1
	
		END DO	

END SELECT
!----------------------------------------------------------------------------------
!Eta=Ny, CB(4)
	DO i=3,Nx+2
	
	!Metricas
	xit(:,i,Ny+4)=xi(:,i-2,Ny-1)	!xiNy+4=xiNy-1
	xit(:,i,Ny+3)=xi(:,i-2,Ny)	!xiNy+3=xiNy
	
	etat(:,i,Ny+4)=eta(:,i-2,Ny-1)	!etaNy+4=etaNy-1
	etat(:,i,Ny+3)=eta(:,i-2,Ny)	!etaNy+3=etaNy
	
	!Batimetria Celdas ficticias
! 	!Celda N+1, zt(N+3,j)
! 	ZxiN=(zt(i,Ny+2)-zt(i,Ny+1))/dxi
! 	ZetaN=(zt(i,Ny+2)-zt(i,Ny+1))/deta
! 	ZyN=ZxiN*xit(2,i,Ny+2)+ZetaN*etat(2,i,Ny+2)
! 	zt(i,Ny+3)=zt(i,Ny+2)+ZyN/(xit(2,i,Ny+3)/dxi+etat(2,i,Ny+3)/deta)
! 
! 	!Celda N+2, zt(N+4,j)
! 	ZxiN1=(zt(i,Ny+3)-zt(i,Ny+2))/dxi
! 	ZetaN1=(zt(i,Ny+3)-zt(i,Ny+2))/deta
! 	ZyN1=ZxiN1*xit(2,i,Ny+3)+ZetaN1*etat(2,i,Ny+3)
! 	zt(i,Ny+4)=zt(i,Ny+3)+ZyN/(xit(2,i,Ny+4)/dxi+etat(2,i,Ny+4)/deta)
	zt(i,Ny+4)=zt(i,Ny+2)
	zt(i,Ny+3)=zt(i,Ny+2)
	END DO

SELECT CASE (CB(4))
	
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
! 		IF (GA4==1) THEN
! 		call genabsNeta_1(Nx,Ny,Fr2,deta,etaL3,Nsenal3,h03,t,dt,qn,zt,eta,qA4,zA4)
! 		ELSE IF (GA4==2) THEN
! 		call genabsNeta_2(Nx,Ny,Fr2,deta,qs3,hs3,Nsenal3,t,dt,qn,zt,eta,qA4,zA4)
! 		ELSE IF (GA4==3) THEN   
! 		call genabsNeta_3(Nx,Ny,Fr2,deta,qs3,us3,Nsenal3,t,dt,qn,zt,eta,qA4,zA4)
! 		END IF
! 
! 		DO i=3,Nx+2
! 	
! 		qt(1,i,Ny+4)=qA4(1,i-2)		!h-1=h2
! 		qt(1,i,Ny+3)=qA4(1,i-2)		!h0=h1
! 		
! 		qt(2,i,Ny+4)=qA4(2,i-2)		!u-1=u2
! 		qt(2,i,Ny+3)=qA4(2,i-2)		!u0=-u1
! 		
! 		qt(3,i,Ny+4)=qA4(3,i-2)		!v-1=v2
! 		qt(3,i,Ny+3)=qA4(3,i-2)		!v0=-v1
! 		
! 		END DO
	CASE(5) 
	 		
		if (IO3==2) then !Outflow, !Cota fija aguas afuera 
		!call OUTFLOWN(fopt,Cf,Coef,Nx,Ny,Fr2,dxi,hs2,Nsenal2,t,dt,qn,zt,xi,qA2,zA2)
		!NO ESTA PROGRAMADA
		
		else
		call INFLOWN_eta(pasoRK,fopt,Cf,Coef,Nx,Ny,Fr2,deta,qsx4,qsy4,etas4,timeS4,Nsenal4,t,dt,qn,zt,eta,xi,qA4,zA4)
		end if
		
		
		
		DO i=3,Nx+2

		qt(1,i,Ny+4)=qA4(1,i-2)		!h-1=h2
		qt(1,i,Ny+3)=qA4(1,i-2)		!h0=h1
		
		qt(2,i,Ny+4)=qA4(2,i-2)		!u-1=-u2
		qt(2,i,Ny+3)=qA4(2,i-2)		!u0=-u1
		
		qt(3,i,Ny+4)=qA4(3,i-2)		!v-1=-v2
		qt(3,i,Ny+3)=qA4(3,i-2)		!v0=-v1
	
		END DO	
END SELECT


END SUBROUTINE bcs