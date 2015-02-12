program main
  use mpi
  use mpi_surf
  use global_variables
  implicit none
  
  real(kind=8)::time_start, time_finish, time_estim
  logical :: dir_exists
  call mpi_init(ierror)
    call init_par
    
    call init_seq
    
    call massbalance
    
    call outputmat_par
    
    print*, 'Simulacion inicializada'
    
    time_start=mpi_wtime()
    do while (treal<= tfinal)
      call setdt
      
      call solver1
      
      t = t +dt
      treal = treal+dt*L/U
      it = it + 1
      
      call massbalance
      if (print_out.and.myrank==master) then
	!Print iteration information on the screen
	print*, 'dt= ', dt
	write(*,199) it
	199 FORMAT ('Iteracion= ',I10)
	print*, 'Time= ', t

	!Adimensional Mass Balance Verification
	
	write(*,170) Pvol
	170 FORMAT ('%Volumen= ',F7.3 ,'%')
	
	print*,'nxi=',Nbx,'neta=',Nby
	!Pvol1=Pvol
	time_finish = mpi_wtime()
	write(*,171) time_finish-time_start
	171 format ('Time Elapsed = ',f10.1,' seconds.')

	time_estim=(time_finish-time_start)*tfinal/treal
	write(*,172) time_estim
	172 format ('Time Estimated = ',f10.1,' seconds.')
	print*,'------------------------------'	
      end if
      
      call outputmat_par
      
      qold_global = qnew_global
    end do
    
  !save stats data
  if (myrank==master) then
    inquire(file='stats.txt', exist=dir_exists)
    if (.not. dir_exists) then
      open(1,file='stats.txt', status='new', action='write')
    else
      open(1,file='stats.txt',status='old',action='write',position='append')
    end if
    write(1,'(I5,",",E20.14)') nproc, time_finish-time_start
  end if
    
  call mpi_finalize(ierror)
  
end program