!Rutina que lee la batimetria y la guarda en x_global, y_global, z_global
SUBROUTINE input_geom

!Programa que crea batimetria para ir probando el codigo, 
!la idea es usar una malla regular para las pruebas
USE global_variables
USE geometries
USE senales
implicit none

real (kind=8), dimension(:), allocatable :: xaux, yaux, zaux
integer	:: i,j
allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))

SELECT CASE (int(batiopt))
  CASE(0)
    open(unit=2,file=batiname(1),form='unformatted')
    read(2) ((x_global(i,j),j=1,Nby),i=1,Nbx)
    close(unit=2)
    
    open(unit=2,file=batiname(2),form='unformatted')
    read(2) ((y_global(i,j),j=1,Nby),i=1,Nbx)
    close(unit=2)
    
    open(unit=2,file=batiname(3),form='unformatted')
    read(2) ((z_global(i,j),j=1,Nby),i=1,Nbx)
    close(unit=2)
  CASE(1)
    open(unit=2,file=batiname(1))
    read(2,*) ((x_global(i,j),j=1,Nby),i=1,Nbx)
    close(unit=2)
    
    open(unit=2,file=batiname(2))
    read(2,*) ((y_global(i,j),j=1,Nby),i=1,Nbx)
    close(unit=2)
    
    open(unit=2,file=batiname(3))
    read(2,*) ((z_global(i,j),j=1,Nby),i=1,Nbx)
    close(unit=2)
  CASE(2)	
    allocate(xaux(Nbx*Nby),yaux(Nbx*Nby), zaux(Nbx*Nby))
    open(unit=2,file=batiname(1))
    Do i=1,Nbx*Nby
      read(2,*) xaux(i), yaux(i), zaux(i)
    End Do
    close(unit=2)
    do i=1,Nbx,1; do j=1,Nby,1
      x_global(i,j)=xaux(j+(i-1)*Nby)
      y_global(i,j)=yaux(j+(i-1)*Nby)
      z_global(i,j)=zaux(j+(i-1)*Nby)
    end do; end do
    deallocate(xaux,yaux,zaux)
  CASE(3)
    allocate(xaux(Nbx*Nby),yaux(Nbx*Nby), zaux(Nbx*Nby))
    open(unit=2,file=batiname(1))
    Do i=1,Nbx*Nby
      read(2,*) xaux(i), yaux(i), zaux(i)
    End Do
    close(unit=2)
    do j=1,Nby,1; do i=1,Nbx,1
      x_global(i,j)=xaux(i+(j-1)*Nbx)
      y_global(i,j)=yaux(i+(j-1)*Nbx)
      z_global(i,j)=zaux(i+(j-1)*Nbx)
    end do; end do
    deallocate(xaux,yaux,zaux)
END SELECT

END SUBROUTINE input_geom
