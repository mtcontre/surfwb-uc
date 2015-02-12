SUBROUTINE tau(tipo,C,Fr2,h,u,taub)

implicit none
integer:: tipo
real (kind=8):: C, Fr2,h,u,taub

IF (tipo==1) THEN !Manning
taub=1.0D0/Fr2*C**2.0D0*u*abs(u)/h**(1.0D0/3.0D0)

ELSE IF (tipo==2) THEN !Chezy

taub=1.0D0/Fr2*1.0D0/(C**2.0D0)*u*abs(u)

ELSE IF (tipo==3) THEN !Sampson
taub=h*u*C

ELSE IF (tipo==4) THEN
taub = C*u*abs(u)
ELSE
print*, 'No existe tipo de friccion'
stop

END IF



END SUBROUTINE tau

SUBROUTINE tauU(tipo,C,Fr2,h,u,v,taub)
  implicit none
  integer:: tipo
  real (kind=8):: C, Fr2,h,u,v,taub

  IF (tipo==1) THEN !Manning
  taub=1.0D0/Fr2*C**2.0D0*u*sqrt(u**2.0D0+v**2.0D0)/h**(1.0D0/3.0D0)

  ELSE IF (tipo==2) THEN !Chezy

  taub=1.0D0/Fr2*1.0D0/(C**2.0D0)*u*sqrt(u**2.0D0+v**2.0D0)

  ELSE IF (tipo==3) THEN !Sampson
  taub=h*sqrt(u**2.0D0+v**2.0D0)*C

  ElSE IF (tipo==4) THEN  
  taub = C*u*sqrt(u*u + v*v)
  ELSE
  print*, 'No existe tipo de friccion'
  stop

  END IF

END SUBROUTINE tauU