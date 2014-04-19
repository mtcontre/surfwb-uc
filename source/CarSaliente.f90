SUBROUTINE Saliente(zmin,jmin,hminL,AL,TL,qL)

!Funcion que calcula la altura de aguas en un punto en la seccion
!resuelve implicitamente: Rmas-(uA+2C)epx=0
USE global_variables
USE geometries
real (kind=8),dimension(Nby,2)::seccion
real (kind=8):: zmin,jmin,hminL,hL,AL,TL,PL,Rh,Dh,qL

seccion(:,1)=y_global(Nbx,:)
seccion(:,2)=z_global(Nbx,:)

!hmin es la altura de agua en el punto mas bajo de la seccion

call PropSecNat(seccion,Nbx,hminL,AL,PL,TL,Rh,Dh)

!qL
qL=0.0D0
Do i=1,Nby
qL=qL+qold_global(1,Nbx,j)*sqrt(qold_global(2,Nbx,j)**2.0D0+qold_global(3,Nbx,j)**2.0D0)
end Do
!hL
END SUBROUTINE Saliente