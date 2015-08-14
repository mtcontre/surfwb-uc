PROGRAM MAIN !Ponerle un nombre decente

  USE global_variables !Use the global variables module
  USE geometries	!Use geometries module
  USE senales	!Use senales module for inflow or outflow BC
  USE TimeSeries
  USE MPI
  USE MPI_surf
  
  implicit none 

  
  !Declare variables that are used here and only here IF it is necessary

  integer::i,j,k
  real (kind=8) :: mindxdy,minxieta,maxUC, maxV, maxC, zmax, maxS1, maxS2, pVol1,dtreal
  real (kind=8), dimension(2) :: xieta, loc
  real (kind=8), dimension(:), allocatable:: dxdy
  real (kind=8) :: time_start,time_finish, time_estim
  integer :: clock_start, clock_rate, clock_finish
  character (len=255) :: formatstring
  logical :: fexists
  
  CALL MPI_INIT(ierror) 
  
  g=9.812D0
  hmin=1.0E-7
  it=0.0D0
  dt=0.0D0
  t=0.0D0
  treal=0.0D0

  call init
  
  IF (outopt==1) THEN
    CALL outputmat
    CALL outputgauges(treal)
  END IF
  
  CALL massbalance	!Calculates the initial volume and mass of water

  CALL system_clock(clock_start,clock_rate)  
  time_start = real(clock_start,kind=8) /real(clock_rate,kind=8)    
  
  
  
  DO WHILE(treal<=tfinal+1.d-10)
    !1. Calculate time step with the CFL condition
    CALL tstep(dtreal)    

    !2.Calling main_solver, which solves the 4 stages of RK method, and calculates qnew_global(h,u,v)
    IF (fopt==0) THEN
      IF (rk==1) THEN		
	CALL solver1 !Solver RK4, 	sin fricción	
      else
	CALL solver2 !Solver RK2, 	sin fricción
      END IF
    else
      IF (rk==1) THEN
    	CALL solverf4 !Solver RK4, 	con fricción
      else
	CALL solverf2 !Solver RK2, 	con fricción
      END IF
    END IF
    
    treal=treal+dtreal
    t=treal*U/L
    it=it+1
    IF (print_out.and.myrank==0) THEN
      CALL massbalance
      
      !Print iteration information on the screen
      print*, 'dt= ', dt
      write(*,199) it
      199 FORMAT ('Iteracion= ',I3)
      print*, 'Time= ', t
      
      !Adimensional Mass Balance VerIFication      
      write(*,170) Pvol
      170 FORMAT ('%Volumen= ',F7.3 ,'%')      
      print*,'nxi=',Nbx,'neta=',Nby
      
      CALL system_clock(clock_finish, clock_rate)
      time_finish = real(clock_finish ,kind=8) / real(clock_rate,kind=8)
      write(*,171) time_finish-time_start
      171 format ('Time Elapsed = ',1f10.1,' seconds.')

      time_estim=(time_finish-time_start)*tfinal/treal
      write(*,172) time_estim
      172 format ('Time Estimated = ',f10.1,' seconds.')
      print*,'------------------------------'
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	CALL LPT () ! Lagrangian Particle Tracking
    !3.Dimensionalize results and write results into a file
    IF (outopt==1) THEN
	  CALL outputmat
	  CALL outputgauges(treal)
    END IF
    
    !Update global variables (adimensionalized)
    qold_global=qnew_global    
  END DO
  
  ! CALL cpu_time(time_finish)
  CALL MPI_Barrier(MPI_COMM_WORLD,ierror)
  IF (myrank==0) THEN
    CALL system_clock(clock_finish, clock_rate)
	time_finish = real(clock_finish ,kind=8) / real(clock_rate,kind=8)
    print *, 'Simulation Ended'
    print *, 'Final Time of Computation= ', t
    print *, 'tfinal', tfinal
    print *, 'treal+dtreal',treal+dtreal
    print *, 'Iteraciones', it
    print *, 'Time Elapsed = ',time_finish-time_start,' seconds.'
  END IF
  CALL MPI_finalize(ierror)
END PROGRAM MAIN


