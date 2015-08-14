subroutine tstep(dtreal)
  !(JGalazM,2013) calculates dtreal satisfying the cfl condition, 
  !so if there are boundary conditions with a given time step, it always 
  !passes through t_i, for instance, considering the xi0 boundary
  !and given t_i=0+i*dtxi0 as the timesteps for that boundary,
  !then tnext=min{t+dt,0+k*dtxi0}, withk = min{k : k*dtxi0>t}
  
  ! min grid size, max celerity inside the domain
  
  USE custombc
  USE global_variables
  USE MPI
  USE MPI_SURF
  implicit none
  real (kind=8) ::minxieta,maxS1,maxS2,maxUC,dtreal
  integer :: nt1,nt2
  
  
  minxieta=minval((/dxi,deta/))
  maxS1=maxval(S1_global)
  maxS2=maxval(S2_global)  
  maxUC=maxval((/maxS1,maxS2/))


  
  ! check celerities in ghost cells
  if (flagxi0.eq.1) then    
    maxUC=maxval((/maxUC,maxval(Sxi0)/))
  end if  
  if (flagxiN.eq.1) then  
    maxUC=maxval((/maxUC,maxval(SxiN)/))
  end if
  if (flageta0.eq.1) then
    maxUC=maxval((/maxUC,maxval(Seta0)/))
  end if
  if (flagetaN.eq.1) then
    maxUC=maxval((/maxUC,maxval(SetaN)/))
  end if
  !reduce maxUC over all procs
  call mpi_allreduce(maxUC,maxUC,1,mpi_double_precision,mpi_max,comm2d,ierror)
    
  
  !CFL condition
  dt=CFL*minxieta/maxUC		!This dt is adimensional
  dtreal=dt*L/U
  ! now fix dtreal so it satisfies t+dtreal<=nt*dtboundary
!   print*,'\n------treal=',treal
!   print*,'------dtreal=',dtreal
  if (flagxi0.eq.1) then
    nt1=int(t/dt_xi0g1+1e-6)+1!satisfies nt1=min{n : n*dt_xi0g1>t}
    nt2=int(t/dt_xi0g2+1e-6)+1
    dtreal=minval((/dtreal,nt1*dt_xi0g1-treal,nt2*dt_xi0g2-treal/))
!     print*,'------xi0dtreal=',dtreal,nt1*dt_xi0g1-treal
  end if
  if (flagxiN.eq.1) then
    nt1=int(t/dt_xiNg1+1e-6)+1!satisfies nt1=min{n : n*dt_xi0g1>t}
    nt2=int(t/dt_xiNg2+1e-6)+1
    dtreal=minval((/dtreal,nt1*dt_xiNg1-treal,nt2*dt_xiNg2-treal/))
!     print*,'------xiNdtreal=',nt1*dt_xiNg1-treal,nt2*dt_xiNg2-treal
  end if
  if (flageta0.eq.1) then
    nt1=int(t/dt_eta0g1+1e-6)+1!satisfies nt1=min{n : n*dt_xi0g1>t}
    nt2=int(t/dt_eta0g2+1e-6)+1
    dtreal=minval((/dtreal,nt1*dt_eta0g1-treal,nt2*dt_eta0g2-treal/))
!     print*,'------eta0dtreal=',nt1*dt_eta0g1-treal
  end if
  if (flagetaN.eq.1) then
    nt1=int(t/dt_etaNg1+1e-6)+1!satisfies nt1=min{n : n*dt_xi0g1>t}
    nt2=int(t/dt_etaNg2+1e-6)+1
    dtreal=minval((/dtreal,nt1*dt_etaNg1-treal,nt2*dt_etaNg2-treal/))
!     print*,'------etaNdtreal=',dtreal
  end if
  
  if (dit==-1) then
    print_out = .False.
    if (treal+dtreal>=(nitout+1)*dtout+tinit) then
      dtreal = (nitout+1)*dtout+tinit-treal
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
