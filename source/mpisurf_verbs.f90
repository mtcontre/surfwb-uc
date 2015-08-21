subroutine it_verbose
  use mpi_surf
  use mpi
  use global_variables
  use performance_stats
  
  call massbalance
  if (mod(it,10)==0.and.myrank==master) then
    !Print iteration information on the screen
    print*, 'dt= ', dt
    write(*,199) it
    199 FORMAT ('Iteracion= ',I10)
    print*, 'Time= ', t
    !Adimensional Mass Balance Verification    
    write(*,170) Pvol
    170 FORMAT ('%Volumen= ',F7.3 ,'%')	
    time_finish= mpi_wtime()
    write(*,171) time_finish-time_start
    171 format ('Time Elapsed = ',f10.1,' seconds.')
    time_estim=(time_finish-time_start)*tfinal/treal
    write(*,172) time_estim
    172 format ('Time Estimated = ',f10.1,' seconds.')
    print*,'------------------------------'
  end if
end subroutine