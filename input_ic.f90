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


SELECT CASE (int(initqopt))
  CASE(0)
    open(unit=2,file=initqname(1),form='unformatted')
    read(2) ((qold_global(1,i,j),j=1,Nby),i=1,Nbx)
    close(unit=2)
    
    open(unit=2,file=initqname(2),form='unformatted')
    read(2) ((qold_global(2,i,j),j=1,Nby),i=1,Nbx)
    close(unit=2)
    
    open(unit=2,file=initqname(3),form='unformatted')
    read(2) ((qold_global(3,i,j),j=1,Nby),i=1,Nbx)
    close(unit=2)
  CASE(1)
    open(unit=2,file=initqname(1))
    read(2,*) ((qold_global(1,i,j),j=1,Nby),i=1,Nbx)
    close(unit=2)
    
    open(unit=2,file=initqname(2))
    read(2,*) ((qold_global(2,i,j),j=1,Nby),i=1,Nbx)
    close(unit=2)
    
    open(unit=2,file=initqname(3))
    read(2,*) ((qold_global(3,i,j),j=1,Nby),i=1,Nbx)
    close(unit=2)
  CASE (2)
    allocate(hin(Nbx*Nby),uin(Nbx*Nby), vin(Nbx*Nby))
    open(unit=99,file='data/initq.dat')
    Do i=1,Nbx*Nby
      read(99,*) hin(i), uin(i), vin(i)
    End Do
    close(unit=99)
    do i=1,Nbx,1; do j=1,Nby,1
      qold_global(1,i,j)=hin(j+(i-1)*Nby)
      qold_global(2,i,j)=uin(j+(i-1)*Nby)
      qold_global(3,i,j)=vin(j+(i-1)*Nby)
    end do; end do
    deallocate(hin,uin,vin)  
  CASE (3)
    allocate(hin(Nbx*Nby),uin(Nbx*Nby), vin(Nbx*Nby))
    open(unit=99,file='data/initq.dat')
    Do i=1,Nbx*Nby
      read(99,*) hin(i), uin(i), vin(i)
    End Do
    close(unit=99)
    do i=1,Nby,1; do j=1,Nbx,1
      qold_global(1,i,j)=hin(i+(j-1)*Nbx)
      qold_global(2,i,j)=uin(i+(j-1)*Nbx)
      qold_global(3,i,j)=vin(i+(j-1)*Nbx)
    end do; end do
    deallocate(hin,uin,vin)

END SELECT
! print*,'h------------'
! do i =1,3
!   print*,qold_global(1,i,1),qold_global(1,i,2),qold_global(1,i,3)
! end do
! pause
! 
! print*,'u------------'
! do i =1,3
!   print*,qold_global(2,i,1),qold_global(2,i,2),qold_global(2,i,3)
! end do
! pause
! 
! print*,'v------------'
! do i =1,3
!   print*,qold_global(3,i,1),qold_global(3,i,2),qold_global(3,i,3)
! end do
! pause


END SUBROUTINE input_ic