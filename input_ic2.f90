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

real (kind=8), dimension(:), allocatable:: told !vector de tiempo de la simulacion que ya estaba
integer ::st,oldcase,oldnx,oldny,olda1,olda2,oldnit,olddit,it0fromold
integer :: temp
character (len=200) :: sol2dI !elindice I del ultimo archivo por cargar como condicion inicial

allocate (qnew_global(3,Nbx,Nby), qold_global(3,Nbx,Nby), &
	  qreal_global(3,Nbx,Nby),q0_global(3,Nbx,Nby), V_global(Nbx,Nby), C_global(Nbx,Nby), &
	  VC(Nbx,Nby),S1_global(Nbx,Nby),S2_global(Nbx,Nby), STAT = error)
!GA
allocate(qA1(3,Nby),qA2(3,Nby),qA3(3,Nbx),qA4(3,Nbx),zA1(Nby),zA2(Nby),zA3(Nbx),zA4(Nbx))
allocate(hin(Nbx*Nby),uin(Nbx*Nby), vin(Nbx*Nby))

inquire(file='results/param.dat',exist=st)
print*,'st=',st
pause
print *, .not. st


if (.not. st) go to 10
pause
  open(unit=97,file='results/param.dat',iostat=st)
  read(97,*), oldcase,oldnx,oldny,olda1,olda2,oldnit,olddit
  open(unit=98,file='results/Time999.dat',iostat=st)  
if (st>0) go to 10  
  it0fromold=int(oldnit/olddit)*olddit !calculo la iteracion desde la que parto, para guardar consistentemente
  allocate(told(it0fromold/olddit))
  read(98,*),(told(temp),temp=1,it0fromold/olddit)
  
if (maxval(told) .lt. tfinal) go to 11
  write(sol2dI,*) it0fromold!sol2dI=trim(it0fromold)
  sol2dI='results/SOL2D.'//sol2dI//'.dat'
  open(unit=99,file=sol2dI)
  Ntot=Nbx*Nby  
  do i=1,Ntot
    read(99,*),temp,temp,temp,hin(i),uin(i),vin(i) 
  end do
  go to 11
  
10  Ntot=Nbx*Nby
    pause
    print*,'Reading initq.dat . . .'
    !Lee batimetria Leandro
    open(unit=99,file='data/initq.dat')
    Do i=1,Ntot
      read(99,*) temp,temp,temp, hin(i), uin(i), vin(i)
    End Do
    close(unit=99)
    do i=1,Nbx,1; do j=1,Nby,1
    qold_global(1,i,j)=hin(j+(i-1)*Nby)+0.0D0
    qold_global(2,i,j)=uin(j+(i-1)*Nby)+0.0D0
    qold_global(3,i,j)=vin(j+(i-1)*Nby)+0.0D0
    end do; end do

11 print*, 'SOL2D in results/, delete these to update the solution'
12 print *,'Cond Inicial'

 
!primero leo el vector tiempo
!si el ultimo valor es menor q tfinal (del inputcontrol) entonces 
!la condicion inicial va a ser el ultimo valor
!sino termina nomas
	
	
	

! if st 	
! print*,'Reading from last solution'
! read(99,*)( t(i),i
! 	
! END SELECT

!q0_global=qold_global

END SUBROUTINE input_ic