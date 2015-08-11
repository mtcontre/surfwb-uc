SUBROUTINE bcs(fopt,Cf,Coef,pasoRK,t,dt,Fr2,caso,qn,z,Nx,Ny,CB,xi,eta,Jac,dxi,deta,qt,xit,etat,zt)

  !2nd Order Boundary Conditions!
  !bcs(state variables at n time,Nbx,Nby,BoundaryConditions, xi_metrics, eta_metrics, 
  !    x_xi,y_xi,x_eta,y_eta,matrix with ghost cells)

  ! USE global_variables
  !USE geometries
  USE senales
  USE coords
  implicit none
  integer :: Nx,Ny, caso,nt1,nt2
  real (kind=8), dimension(3,Nx,Ny)	::qn
  real (kind=8), dimension(2,Nx,Ny)	::xi,eta
  real (kind=8), dimension(3,Nx+4,Ny+4)	::qt
  real (kind=8), dimension(Nx+4,Ny+4)	::zt
  real (kind=8), dimension(Nx,Ny)		::z, Coef, Jac
  real (kind=8), dimension(2,Nx+4,Ny+4)	::xit,etat
  integer,dimension(4) :: CB  
  real (kind=8)	:: t,dt,Fr2,dxi,deta,Zx1,Zxi1, Zeta1, Zx0, Zxi0, &
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


  !Metricas simetricas a las de adentro
  !Batimetria tambien simetrica?? o igual a la ultima??
  DO j=3,Ny+2
    xit(:,1,j)=xi(:,2,j-2)		!xi-1=xi2
    xit(:,2,j)=xi(:,1,j-2)		!xi0=xi12

    etat(:,1,j)=eta(:,2,j-2)	!eta-1=eta2
    etat(:,2,j)=eta(:,1,j-2)	!eta0=eta1

    zt(1,j)=zt(4,j)
    zt(2,j)=zt(3,j)
  END DO


  SELECT CASE (CB(1))
    CASE(1) !Solid Wall
      DO j=3,Ny+2
	qt(1,1,j)=qn(1,2,j-2)		!h-1=h2
	qt(1,2,j)=qn(1,1,j-2)		!h0=h1

	qt(2,1,j)=-qn(2,2,j-2)		!u-1=-u2
	qt(2,2,j)=-qn(2,1,j-2)		!u0=-u1

	qt(3,1,j)=-qn(3,2,j-2)		!v-1=-v2
	qt(3,2,j)=-qn(3,1,j-2)		!v0=-v1
      END DO
      
    CASE(2) !Periodic
      DO j=3,Ny+2
	qt(1,1,j)=qn(1,Nx-1,j-2)	!h-1=hNx-1
	qt(1,2,j)=qn(1,Nx,j-2)		!h0=hNx

	qt(2,1,j)=qn(2,Nx-1,j-2)	!u-1=uNx-1
	qt(2,2,j)=qn(2,Nx,j-2)		!u0=uNx

	qt(3,1,j)=qn(3,Nx-1,j-2)	!v-1=vNx-1
	qt(3,2,j)=qn(3,Nx,j-2)		!v0=vNx
      END DO

    CASE(3) !Free Exit	
      DO j=3,Ny+2
	qt(1,1,j)=qn(1,2,j-2)		!h-1=h2
	qt(1,2,j)=qn(1,1,j-2)		!h0=h1

	qt(2,1,j)=qn(2,2,j-2)		!u-1=-u2
	qt(2,2,j)=qn(2,1,j-2)		!u0=-u1

	qt(3,1,j)=qn(3,2,j-2)		!v-1=-v2
	qt(3,2,j)=qn(3,1,j-2)		!v0=-v1
      END DO	
  END SELECT
  !---------------------------------------------------------------------------

  DO j=3,Ny+2
    xit(:,Nx+4,j)=xi(:,Nx-1,j-2)	!xiNx+4=xiNx-1
    xit(:,Nx+3,j)=xi(:,Nx,j-2)		!xiNx+3=xiNx

    etat(:,Nx+4,j)=eta(:,Nx-1,j-2)	!etaNx+4=etaNx-1
    etat(:,Nx+3,j)=eta(:,Nx,j-2)	!etaNx+3=etaNx

    zt(Nx+4,j)=zt(Nx+2,j)
    zt(Nx+3,j)=zt(Nx+2,j)
  END DO

  SELECT CASE (CB(2))
    CASE(1) !Solid Wall 
      DO j=3,Ny+2
	qt(1,Nx+4,j)=qn(1,Nx-1,j-2)	!hNx+4=hNx-1
	qt(1,Nx+3,j)=qn(1,Nx,j-2)	!hNx+3=hNx

	qt(2,Nx+4,j)=-qn(2,Nx-1,j-2)	!uNx+4=-uNx-1
	qt(2,Nx+3,j)=-qn(2,Nx,j-2)	!uNx+3=-uNx

	qt(3,Nx+4,j)=-qn(3,Nx-1,j-2)	!vNx+4=-vNx-1
	qt(3,Nx+3,j)=-qn(3,Nx,j-2)	!vNx+3=-vNx
      END DO

    CASE(2) !Periodic
      DO j=3,Ny+2
	qt(1,Nx+4,j)=qn(1,2,j-2)	!hNx+4=h2
	qt(1,Nx+3,j)=qn(1,1,j-2)	!hNx+3=h1

	qt(2,Nx+4,j)=qn(2,2,j-2)	!uNx+4=u2
	qt(2,Nx+3,j)=qn(2,1,j-2)	!uNx+3=u1

	qt(3,Nx+4,j)=qn(3,2,j-2)	!vNx+4=v2
	qt(3,Nx+3,j)=qn(3,1,j-2)	!vNx+3=v1
      END DO
    CASE(3) !Free Exit
      DO j=3,Ny+2
	qt(1,Nx+4,j)=qn(1,Nx-1,j-2)	!hNx+4=hNx-1
	qt(1,Nx+3,j)=qn(1,Nx,j-2)		!hNx+3=hNx

	qt(2,Nx+4,j)=qn(2,Nx-1,j-2)	!uNx+4=uNx-1
	qt(2,Nx+3,j)=qn(2,Nx,j-2)		!uNx+3=uNx

	qt(3,Nx+4,j)=qn(3,Nx-1,j-2)	!vNx+4=vNx-1
	qt(3,Nx+3,j)=qn(3,Nx,j-2)		!vNx+3=vNx
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
    zt(i,1)=zt(i,3)
    zt(i,2)=zt(i,3)	
  END DO

  SELECT CASE (CB(3)) 
    CASE(1) !Solid Wall
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
  END SELECT
  !----------------------------------------------------------------------------------
  !Eta=Ny, CB(4)
  DO i=3,Nx+2
    !Metricas
    xit(:,i,Ny+4)=xi(:,i-2,Ny-1)	!xiNy+4=xiNy-1
    xit(:,i,Ny+3)=xi(:,i-2,Ny)	!xiNy+3=xiNy

    etat(:,i,Ny+4)=eta(:,i-2,Ny-1)	!etaNy+4=etaNy-1
    etat(:,i,Ny+3)=eta(:,i-2,Ny)	!etaNy+3=etaNy

    zt(i,Ny+4)=zt(i,Ny+1)
    zt(i,Ny+3)=zt(i,Ny+2)
  END DO

  SELECT CASE (CB(4))
    CASE(1)
      DO i=3,Nx+2
	qt(1,i,Ny+4)=qn(1,i-2,Ny-1)		!hNy+4=hNy-1
	qt(1,i,Ny+3)=qn(1,i-2,Ny)		!hNy+3=hNy

	qt(2,i,Ny+4)=-qn(2,i-2,Ny-1)	!uNy+4=-uNy-1
	qt(2,i,Ny+3)=-qn(2,i-2,Ny)		!uNy+3=-uNy

	qt(3,i,Ny+4)=-qn(3,i-2,Ny-1)	!vNy+4=-vNy-1
	qt(3,i,Ny+3)=-qn(3,i-2,Ny)		!vNy+3=-vNy
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
  END SELECT
END SUBROUTINE bcs

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