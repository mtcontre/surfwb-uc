SUBROUTINE angulo
!Subrutina que calcula el angulo de la normal al borde de las celdas con respecto a un eje
!Borde 1  y Borde 2 (xi=0, xi=Nx) calcula alfa con respecto a eje x
!Borde 3  y Borde 4 (eta=0, eta=Ny) calcula alfa con respecto a eje y
USE coords
USE global_variables
USE geometries

integer:: i,j
real (kind=8):: dx, dy, pi

allocate(angulo1(Nby),angulo2(Nby),angulo3(Nbx),angulo4(Nbx))

pi=3.14159265D0

Do j=2,(Nby-1)
!Borde 1
dx=x_global(1,j+1)-x_global(1,j)
dy=y_global(1,j+1)-y_global(1,j)

angulo1(j)=atan(dx/dy)

!Borde Nbx
dx=x_global(Nbx,j+1)-x_global(Nbx,j)
dy=y_global(Nbx,j+1)-y_global(Nbx,j)
angulo2(j)=atan(dx/dy)
End Do
angulo1(1)=angulo1(2)
angulo1(Nby)=angulo1(Nby-1)
angulo2(1)=angulo2(2)
angulo2(Nby)=angulo2(Nby-1)

Do i=2,(Nbx-1)
!Borde 3
dx=x_global(i+1,1)-x_global(i,1)
dy=y_global(i+1,1)-y_global(i,1)
angulo3(i)=atan(dy/dx)

!Borde Nby
dx=x_global(i+1,Nby)-x_global(i,Nby)
dy=y_global(i+1,Nby)-y_global(i,Nby)
angulo4(i)=atan(dy/dx)+pi
End Do
angulo3(1)=angulo3(2)
angulo3(Nbx)=angulo3(Nbx-1)
angulo4(1)=angulo4(2)
angulo4(Nbx)=angulo4(Nbx-1)

END SUBROUTINE angulo
