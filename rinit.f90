!Rutina que inicia las variables,lee datos de entrada, lee batimetria, condiciones iniciales, condiciones de borde y asigna estos datos a variables
!Llama a tranformacion de coordenadas
SUBROUTINE init

USE global_variables
USE geometries
USE senales
USE coords

implicit none

integer :: i,j 
real (kind=8) :: U1, U2

  
!Leer informacion de la simulacion
call input_control

FR2=U**2.0D0/(g*H)
print*, 'Fr2= ', FR2

call input_geom		!Leer batimetria
call input_ic		!Leer condición inicial
call init_TS ! Time Series

print*,'InitFlowfieldOK'
!Coeficiente de friccion: debe venir adimensionalizado
allocate(MCoef(Nbx,Nby))

IF (fopt==1) THEN
  
  IF (fM==1) THEN !Completa matriz de friccion con el valor ingresado en input.dat
    DO i=1,Nbx 
    DO j=1,Nby
    MCoef(i,j)=Coef
    END DO
    END DO
  print*, 'friccionOK'
  ELSE !Lee matriz de friccion de un archivo binario
    !open	(unit=99, file ='friccion.dat', form='unformatted')
    open	(unit=99, file ='data/friction_run31.dat')
    read(99,*) ((MCoef(i,j),i=1,Nbx),j=1,Nby)
    !read	(unit=99) ((MCoef(i,j),i=1,Nbx),j=1,Nby)
    close(unit=99)
    print*, 'friccionOK'
  END IF

ELSE
  Cf=0
  DO i=1,Nbx 
  DO j=1,Nby
  MCoef(i,j)=0.0D0
  END DO
  END DO
 
END IF
      

!Adimensionalizar CI y Bathy
call adimension

!Calcular las metricas
call metrics
!Crea Coordenadas Xi,Eta

allocate(coordxi(Nbx),coordeta(Nby))
 coordxi(1)=dxi/2.0D0
 coordeta(1)=deta/2.0D0
do i=2,Nbx
 coordxi(i)=coordxi(i-1)+dxi
end do
do i=2,Nby
 coordeta(i)=coordeta(i-1)+deta
end do

!Angulo normales a los bordes con respecto a un eje

call angulo

!Velocidades Contravariantes para CFL
do i=1,Nbx; do j=1,Nby

	U1=qold_global(2,i,j)*xi_global(1,i,j)+qold_global(3,i,j)*xi_global(2,i,j)
	U2=qold_global(2,i,j)*eta_global(1,i,j)+qold_global(3,i,j)*eta_global(2,i,j)
	S1_global(i,j)=U1+C_global(i,j)*sqrt(xi_global(1,i,j)**2+xi_global(2,i,j)**2)
	S2_global(i,j)=U2+C_global(i,j)*sqrt(eta_global(1,i,j)**2+eta_global(2,i,j)**2)


end do; end do


print*, 'Simulacion Incializada'

END SUBROUTINE init


SUBROUTINE init_flowfield
!Condiciones Iniciales: Diferentes para cada caso
!Se podrian bajar de archivos, pero por mientras mejor hacerlas aqui

USE global_variables
USE geometries
USE senales
USE time0
implicit none

integer :: mitad, i, j, error, Ntot
real (kind=8), dimension(:), allocatable :: hin, uin, vin
real (kind=8)::	zaux, hz, xmed, xo, yo, hc, hs,r,d, Haux, sigma, hmean, hestanque, &
		p3, c3, yaux, Am, a, ro, ho, eta0, omega,tau,p,S, B, uo, Dsyn, gama, x1, m, z85
real (kind=8) :: D1, D2, D3, D4, D5, D6
integer :: io

allocate (qnew_global(3,Nbx,Nby), qold_global(3,Nbx,Nby), &
	  qreal_global(3,Nbx,Nby),q0_global(3,Nbx,Nby), V_global(Nbx,Nby), C_global(Nbx,Nby), &
	  VC(Nbx,Nby),S1_global(Nbx,Nby),S2_global(Nbx,Nby), STAT = error)
!GA
allocate(qA1(3,Nby),qA2(3,Nby),qA3(3,Nbx),qA4(3,Nbx),zA1(Nby),zA2(Nby),zA3(Nbx),zA4(Nbx))


SELECT CASE (caso)
	CASE(999)
!  	  allocate (qold_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(hin(Nbx*Nby),uin(Nbx*Nby), vin(Nbx*Nby))
 	  Ntot=Nbx*Nby
 	  print*,'Reading initq.dat . . .'
	  !Lee batimetria Leandro
	  open(unit=99,file='data/initq.dat')
	  
	  Do i=1,Ntot
	    read(99,*) hin(i), uin(i), vin(i)
	  End Do
	  close(unit=99)
	  do i=1,Nbx,1; do j=1,Nby,1
	  qold_global(1,i,j)=hin(j+(i-1)*Nby)+0.0D0
	  qold_global(2,i,j)=uin(j+(i-1)*Nby)+0.0D0
	  qold_global(3,i,j)=vin(j+(i-1)*Nby)+0.0D0
	  end do; end do
	  print *, 'Cond Inicial'
END SELECT

!q0_global=qold_global

END SUBROUTINE init_flowfield


SUBROUTINE ADIMENSION
!Function that aplies the adimensionalization to the initial conditions
!Funcion que aplica la adimensionalizacion a las condiciones inciales y a todo

USE global_variables
USE geometries
implicit none
integer :: i,j
real (kind=8)::U1,U2
	do i=1,Nbx; do j=1,Nby
			x_global(i,j)=x_global(i,j)/L
			y_global(i,j)=y_global(i,j)/L
			z_global(i,j)=z_global(i,j)/H
			qold_global(1,i,j)=qold_global(1,i,j)/H
			qold_global(2,i,j)=qold_global(2,i,j)/U
			qold_global(3,i,j)=qold_global(3,i,j)/U
			
			V_global(i,j)=sqrt((qold_global(2,i,j))**2.0D0+(qold_global(3,i,j))**2.0D0)
			
			C_global(i,j)=sqrt(qold_global(1,i,j)/FR2)
			VC(i,j)=V_global(i,j)+C_global(i,j)

			
	end do; end do

END SUBROUTINE ADIMENSION
