!Rutina que lee la batimetria y la guarda en x_global, y_global, z_global


SUBROUTINE input_geom

!Programa que crea batimetria para ir probando el codigo, 
!la idea es usar una malla regular para las pruebas
USE global_variables
USE geometries
USE senales
implicit none

real (kind=8), dimension(:), allocatable :: xaux, yaux, zaux

real (kind=8):: dx,dy,zf,Lx, Ly, dx1, dx2, ri, rext, ds, xo, yo, sigma,r, ho, a, beta

integer	:: i,j,Ntot

SELECT CASE (caso)
CASE(999)	!Ensayo para un perfil de playa
	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx*Nby),yaux(Nbx*Nby), zaux(Nbx*Nby))
	  Ntot=Nbx*Nby
	  !Lee batimetria Leandro
	  open(unit=99,file='data/bathy.dat')
	  Do i=1,Ntot
	    read(99,*) xaux(i), yaux(i), zaux(i)
	  End Do
	  close(unit=99)
	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(j+(i-1)*Nby)+0.00
	  y_global(i,j)=yaux(j+(i-1)*Nby)
	  z_global(i,j)=zaux(j+(i-1)*Nby)+0.00
	  end do; end do
	  print *, 'bathy', x_global(1,1)
	  !Ensayo para caso 2D	  
END SELECT
 


END SUBROUTINE input_geom
! SUBROUTINE input_geom
! 
! USE global_variables
! USE geometries
! 
! implicit none
! 
! integer::i,j
! 
! allocate (x_global(Nbx,Nby),y_global(Nbx,Nby),z_global(Nbx,Nby))
! 
! print*, Nbx, Nby
! open	(unit=2, file ='gridX.dat', form='unformatted')
! read	(unit=2) ((x_global(i,j),i=1,Nbx),j=1,Nby)
! close(unit=2)
! 
! open	(unit=3, file ='gridY.dat', form='unformatted')
! read	(unit=3) ((y_global(i,j),i=1,Nbx),j=1,Nby)
! close(unit=3)
! 
! open	(unit=4,file='gridZ.dat', form='unformatted')
! read	(unit=4) ((z_global(i,j),i=1,Nbx),j=1,Nby)
! close(unit=4)
! 
! END SUBROUTINE input_geom