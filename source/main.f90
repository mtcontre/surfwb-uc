!NSWE Finite Volume Model
!Curvilinear coordinate system
!Well-Balanced, Shock-Capturing
!Approximate Riemann Solver VFRoe-ncv plus 2nd order Hydrostatic Reconstruction
!Semi-implicit friction step

!Developed by Maricarmen Guerra P.

PROGRAM MAIN !Ponerle un nombre decente

  USE global_variables !Use the global variables module
  USE geometries	!Use geometries module
  USE senales	!Use senales module for inflow or outflow BC
  USE TimeSeries
  use custombc
  
  implicit none 

  !Declare variables that are used here and only here if it is necessary

  integer::i,j,k
  real (kind=8) :: mindxdy,minxieta,maxUC, maxV, maxC, zmax, maxS1, maxS2, pVol1,dtreal
  real (kind=8), dimension(2) :: xieta, loc
  real (kind=8), dimension(:), allocatable:: dxdy
  real (kind=8) :: time_start,time_finish, time_estim
  integer :: clock_start, clock_rate, clock_finish
  integer :: omp_nthreads, omp_threadid, omp_get_num_threads, omp_get_thread_num
  character (len=255) :: formatstring
  logical :: fexists
  
  g=9.812D0
  hmin=1.0E-7
  it=0.0D0
  dt=0.0D0
  t=0.0D0!0.0D0
  treal=0.0D0!960.0D0 !0.0D0

  ! Initialize calculation

  call init !init function: reads input data and bathymethry, initialize flowfield, calculates metrics
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	call LPT () ! Lagrangian Particle Tracking - init

  if (outopt==1) then
    call outputmat	!write initial conditions into a file

    call outputgauges(treal)
  else
	  !call outputtec(treal)
  ! 	call outputparaview(treal)
  end if
  
  call massbalance	!Calculates the initial volume and mass of water

  write(*,100)
  100 FORMAT ('Start Solver')
  print*,'------------------------------'

  zmax=maxVal(z_global(:,:))

  allocate(dxdy((Nbx-1)*Nby+Nbx*(Nby-1)))
  k=1
  do i=1,Nbx-1; do j=1,Nby
    dxdy(k)=x_global(i+1,j)-x_global(i,j)
    k=k+1
  end do; end do

  do i=1,Nbx; do j=1,Nby-1
    dxdy(k)=y_global(i,j+1)-y_global(i,j)
    k=k+1
  end do; end do

  !MAIN CYCLE!

  !Do while (it<=1) !to find errors
!   call cpu_time(time_start)
  
  call system_clock(clock_start,clock_rate)
  
  time_start = real(clock_start,kind=8) /real(clock_rate,kind=8)
  
  DO while(treal<=tfinal+1.d-10)
    !1. Calculate time step with the CFL condition
    !    and boundary timestep
    call tstep(dtreal)    
  
    !2.Calling main_solver, which solves the 4 stages of RK method, and calculates qnew_global(h,u,v)
    if (fopt==0) then
      if (rk==1) then	
	
	call solver1 !Solver RK4, 	sin fricción
	
      else
	call solver2 !Solver RK2, 	sin fricción
      end if
    else
      if (rk==1) then
    	call solverf4 !Solver RK4, 	con fricción
      else
	call solverf2 !Solver RK2, 	con fricción
      end if
    end if
    treal=treal+dtreal
    t=treal*U/L
!     t=t+dt
!     treal=t*L/U
      it=it+1
!       if (it==528) then
! 	goout=.true.
!       end if
    if (print_out) then
      !Print iteration information on the screen

      print*, 'dt= ', dt
      write(*,199) it
      199 FORMAT ('Iteracion= ',I3)
      print*, 'Time= ', t

      !Adimensional Mass Balance Verification
      call massbalance
      write(*,170) Pvol
      170 FORMAT ('%Volumen= ',F7.3 ,'%')
      
      print*,'nxi=',Nbx,'neta=',Nby
      !Pvol1=Pvol
! 	      call cpu_time(time_finish)
      call system_clock(clock_finish, clock_rate)
      time_finish = real(clock_finish ,kind=8) / real(clock_rate,kind=8)
      write(*,171) time_finish-time_start
      171 format ('Time Elapsed = ',1f10.1,' seconds.')

      time_estim=(time_finish-time_start)*tfinal/treal
      write(*,172) time_estim
      172 format ('Time Estimated = ',f10.1,' seconds.')
      print*,'------------------------------'
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	call LPT () ! Lagrangian Particle Tracking
    !3.Dimensionalize results and write results into a file
    if (outopt==1) then
	  call outputmat
	  call outputgauges(treal)
    end if
    
    !Actualization of the global variables (adimensionalized)
    qold_global=qnew_global
    
  END DO
! call cpu_time(time_finish)
call system_clock(clock_finish, clock_rate)
      time_finish = real(clock_finish ,kind=8) / real(clock_rate,kind=8)
  print *, 'Simulation Ended'
  print *, 'Final Time of Computation= ', t
  print *, 'tfinal', tfinal
  print *, 'treal+dtreal',treal+dtreal
  print *, 'Iteraciones', it
  print *, 'Time Elapsed = ',time_finish-time_start,' seconds.'
  
  !save execution times
!$omp parallel private(omp_threadid, omp_nthreads)
  omp_nthreads = omp_get_num_threads()
  omp_threadid = omp_get_thread_num()
  
  if (omp_threadid==0) then
    print *, 'Num threads used =', omp_nthreads
    inquire(file='times.txt',exist=fexists)
    
    if (fexists) then
      open(unit=0, file='times.txt', status='old', position='append', action='write')
    else
      open(unit=0, file='times.txt', status='new', action='write')
    end if
    
    formatstring = '(i4,",",i4,",",E20.10,",",E20.10,",",E20.10)'
    write(unit=0,fmt=formatstring) Nbx, omp_nthreads, time_finish-time_start, time_start, time_finish
  end if
!$omp end parallel
END PROGRAM MAIN


