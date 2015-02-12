subroutine setdt
  !(JGalazM,2013) calculates dtreal satisfying the cfl condition, 
  !so if there are boundary conditions with a given time step, it always 
  !passes through t_i, for instance, considering the xi0 boundary
  !and given t_i=0+i*dtxi0 as the timesteps for that boundary,
  !then tnext=min{t+dt,0+k*dtxi0}, withk = min{k : k*dtxi0>t}
  
  ! min grid size, max celerity inside the domain
  
  use mpi
  use mpi_surf
  use global_variables
  implicit none
  real (kind=8) ::minxieta,maxS1,maxS2,maxUC,dtreal, dtbuf
  integer,dimension(2) :: nt1,nt2
  
  minxieta=minval((/dxi,deta/))
  maxS1=maxval(S1_global)
  maxS2=maxval(S2_global)  
  maxUC=maxval((/maxS1,maxS2/))
  if (myrank==0) then

  end if
  !CFL condition
  dt=CFL*minxieta/maxUC		!This dt is adimensional
  dtreal=dt*L/U

!   print*,dtreal,myrank
  call mpi_allreduce(dtreal,dtbuf,1,mpi_double_precision,mpi_min,comm2d,ierror)
  dtreal = dtbuf
  dt = dtreal*U/L
!   print*,dtreal,myrank
  !decide wether write/step to output or not
  if (dit==-1) then
    print_out = .False.
    if (treal+dtreal>=(nitout+1)*dtout+tinit) then
      dtreal = (nitout+1)*dtout+tinit-treal
      call mpi_allreduce(dtreal,dtbuf,1,mpi_double_precision,mpi_min,comm2d,ierror)
      dtreal = dtbuf
      print_out = .True.    
      nitout=nitout+1
    else
      print_out = .False.
    end if    
  else if (mod(it,dit)==0) then
    print_out = .True.
    nitout = it
  end if
 
  dt=dtreal*U/L
end subroutine
