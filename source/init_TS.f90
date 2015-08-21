!Rutina que inicia las coordenadas de los puntos en los cuales se necesitan las series de tiempo

SUBROUTINE init_TS
  !Condiciones Iniciales: Diferentes para cada caso
  !Se podrian bajar de archivos, pero por mientras mejor hacerlas aqui

  USE TimeSeries
  USE geometries, only: x_global, y_global
  use global_variables, only: outdir,indir
  !USE geometries
  !USE senales
  !USE time0
  implicit none
  integer ::error,i,j,l,m
  logical :: f_exists

  ! Nts=1 ! Number of points for time series output
  open(unit=60,file=trim(indir)//'/gauges.dat')
  !coordx-coordy
  read(60,*) Nts
  ! dt_TS=0.05 ! timestep to save Time Series and LPT   
  !allocate (x0(Nts), y0(Nts), i0(Nts), j0(Nts), r(Nts), s(Nts), H0(Nts), Uts(Nts), Vts(Nts), &
  !		x00(Nts), x01(Nts), x10(Nts), x11(Nts), STAT = error)
  allocate (id0(Nts),x0(Nts), y0(Nts), i0(Nts), j0(Nts), r(Nts), s(Nts), H0(Nts), Uts(Nts), Vts(Nts), &
		  x00(3,Nts), x01(3,Nts), x10(3,Nts), x11(3,Nts), STAT = error)
  allocate(m1(2,Nts))

  DO i=1,Nts    
    read(60,*) id0(i),x0(i),y0(i)
    m1(:,i) = minloc((x_global-x0(i))*(x_global-x0(i))+(y_global-y0(i))*(y_global-y0(i)))    
    l = m1(1,i)
    m = m1(2,i)
  END DO
  close(unit=60)
  
  
  !guardar los indices encontrados para revisar despues.
  
  inquire(file=trim(outdir)//'/timeseries/indices.dat',exist=f_exists)
  if (.not. f_exists) then
    open(unit=60, file=trim(outdir)//'/timeseries/indices.dat',status='new',action='write')
  else
    open(unit=60, file=trim(outdir)//'/timeseries/indices.dat',status='replace')
  end if
  
  do i=1,Nts
    write(unit=60,fmt=*) m1(1,i), m1(2,i)
  end do
  close(unit=60)

END SUBROUTINE init_TS

