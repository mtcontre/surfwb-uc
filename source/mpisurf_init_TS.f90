!Rutina que inicia las coordenadas de los puntos en los cuales se necesitan las series de tiempo

SUBROUTINE init_TS
  !Condiciones Iniciales: Diferentes para cada caso
  !Se podrian bajar de archivos, pero por mientras mejor hacerlas aqui

  USE TimeSeries
  USE geometries, only: x_global, y_global
  use global_variables, only: outdir,indir
  use mpi_surf
  use multigrid_surf
  !USE geometries
  !USE senales
  !USE time0
  implicit none
  integer ::error,i,j,l,m,level
  logical :: f_exists

  ! Nts=1 ! Number of points for time series output
  open(unit=600,file=trim(indir)//'/gauges.dat')
  !coordx-coordy
  read(600,*) Nts
  ! dt_TS=0.05 ! timestep to save Time Series and LPT   
  allocate (id0(Nts),x0(Nts), y0(Nts), STAT = error)
  allocate(m1_temp(2,Nts))

  level = 1
  DO i=1,Nts        
    read(600,*) id0(i),x0(i),y0(i)
    m1_temp(:,i) = minloc((geom(level)%X-x0(i))*(geom(level)%X-x0(i))+ &
			  (geom(level)%Y-y0(i))*(geom(level)%Y-y0(i)))
			
  END DO
  close(unit=600)
  

END SUBROUTINE init_TS

