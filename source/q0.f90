! Funcion que completa la condicion inicial U (undisturbed state) que se usa en inflow_xi.f90 e inflow_eta.f90.
! Esta condicion corresponde al primer estado calculado desde dentro del dominio  (RR (oRL, dependiendo del borde)) que se calcula a partir de la invariante de riemann que parte desde dentro del domino y llega justo al borde.

SUBROUTINE q0(borde,paso,epxR,epyR,hR,uR,vR,i)

USE global_variables !iteraciones
USE time0 !q0

integer :: borde, i, paso
real (kind=8):: hR, uR, vR, epxR, epyR



IF (it==0 .AND.paso==1) THEN

  IF (borde==1) THEN
  q0_global(1,1,i)=hR
  q0_global(2,1,i)=uR
  q0_global(3,1,i)=vR
  epxR0=epxR
  epyR0=epyR
  END IF
  
  IF (borde==2) THEN
  q0_global(1,Nbx,i)=hR
  q0_global(2,Nbx,i)=uR
  q0_global(3,Nbx,i)=vR
  END IF
  
  
  IF (borde==3) THEN
  q0_global(1,i,1)=hR
  q0_global(2,i,1)=uR
  q0_global(3,i,1)=vR
  END IF
  
  IF (borde==4) THEN
  q0_global(1,i,Nby)=hR
  q0_global(2,i,Nby)=uR
  q0_global(3,i,Nby)=vR
  END IF

END IF

END SUBROUTINE q0
  
