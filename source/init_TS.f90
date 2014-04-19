!Rutina que inicia las coordenadas de los puntos en los cuales se necesitan las series de tiempo

SUBROUTINE init_TS
!Condiciones Iniciales: Diferentes para cada caso
!Se podrian bajar de archivos, pero por mientras mejor hacerlas aqui

USE TimeSeries
!USE global_variables
!USE geometries
!USE senales
!USE time0
implicit none
integer ::error,i,j

! Nts=1 ! Number of points for time series output
open(unit=60,file='data/gauges.dat')
!coordx-coordy
read(60,*) Nts
! dt_TS=0.05 ! timestep to save Time Series and LPT   
!allocate (x0(Nts), y0(Nts), i0(Nts), j0(Nts), r(Nts), s(Nts), H0(Nts), Uts(Nts), Vts(Nts), &
!		x00(Nts), x01(Nts), x10(Nts), x11(Nts), STAT = error)
allocate (id0(Nts),x0(Nts), y0(Nts), i0(Nts), j0(Nts), r(Nts), s(Nts), H0(Nts), Uts(Nts), Vts(Nts), &
		x00(3,Nts), x01(3,Nts), x10(3,Nts), x11(3,Nts), STAT = error)

DO i=1,Nts    
      read(60,*) id0(i),x0(i),y0(i)
END DO
print *, 'Init Time Series'

END SUBROUTINE init_TS

