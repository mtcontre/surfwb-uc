SUBROUTINE pascua(Nx,Ny)

!Lee la batimetría del pascua y la guarda en un archivo no formateado


real (kind=8), allocatable, dimension(:,:)::x, y, z
integer:: Nx, Ny, i,j


allocate(x(Nx,Ny),y(Nx,Ny),z(Nx,Ny))

open(unit=1,file='gridX_mat.dat')!, form='formatted')
open(unit=2,file='gridY_mat.dat')!, form='formatted')
open(unit=3,file='gridZ_mat.dat')!, form='formatted')


     read(1,*) ((x(i,j),i=1,Nx),j=1,Ny)
     read(2,*) ((y(i,j),i=1,Nx),j=1,Ny)
     read(3,*) ((z(i,j),i=1,Nx),j=1,Ny)

!END DO; END DO;


close(unit=1)
close(unit=2)
close(unit=3)

!Writing Formatted files

print*, x(:,1)

open  (unit = 41, file = 'gridX.dat', form='unformatted')
          write  (unit = 41) ((x(i,j),i=1,Nx),j=1,Ny)
close(unit = 41)

open  (unit = 42, file = 'gridY.dat', form='unformatted')
          write  (unit = 42) ((y(i,j),i=1,Nx),j=1,Ny)
close(unit = 42)

open  (unit = 43, file = 'gridZ.dat', form='unformatted')
          write  (unit = 43) ((z(i,j),i=1,Nx),j=1,Ny)
close(unit = 43)

print*, Nx, Ny

call input_geom

END SUBROUTINE pascua

