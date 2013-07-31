SUBROUTINE input_ic
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

END SUBROUTINE input_ic