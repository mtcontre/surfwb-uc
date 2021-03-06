program main
  use mpi
  use mpi_surf
  use multigrid_surf
  use global_variables
  use performance_stats
  implicit none
  logical :: fexists
  
  call mpi_init(ierror)
    call init
    
    call massbalance
    
    !write initial condition
    if (outopt==1) then      
      call outputmat_par
      call outputgauges(treal)
    end if
   
    time_start=mpi_wtime()
    do while(treal<=tfinal)
      !set a stable delta t
      call setdt      
      !pick one solver from input_control params
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
      
      t=t+dt
      treal=treal+dtreal
      it=it+1	
      !print screen information for this iteration
      call it_verbose
      
      !Lagrangian Particle Tracking
      !call LPT () 
      !print current new results
      if (outopt==1) then
	call outputmat_par
	call outputgauges(treal)
      end if
      
      !Actualization of the global variables (adimensionalized)
      qold_global=qnew_global
    end do
    
    !print the last iteration
    
    if (myrank==master) then
      print *, 'Simulation Ended'
      print *, 'Final Time of Computation= ', t
      print *, 'tfinal', tfinal
      print *, 'treal+dtreal',treal+dtreal
      print *, 'Iteraciones', it
      print *, 'Time Elapsed = ',time_finish-time_start,' seconds.'
      inquire(file='stats.txt',exist=fexists)
      if (fexists) then
	open(unit=100,file='stats.txt',status='old',position='append')
      else
	open(unit=100,file='stats.txt',status='new',action='write')
      end if
      write(unit=100,fmt=*) nproc,nxi,neta,time_finish-time_start
      close(unit=100)
    end if
    
  call mpi_finalize(ierror)
  
end program main

