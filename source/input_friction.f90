subroutine input_friction
  use mpi
  use mpi_surf
  use global_variables
  
  integer ::i,j
  
  allocate(fric_buff(nxi,neta))
  if (fopt==0) then !no friction
    Cf = 0
    Coef =0.0D0
    fric_buff = Coef
  elseif (fopt==1) then
    if (fM==1) then
      fric_buff = Coef
    elseif (fM==2) then
      open(unit=99, file='data/fricmat.dat')
      read(99,*) ((fric_buff(i,j),j=1,neta),i=1,nxi)
      close(unit=99)
    end if    
  end if


end subroutine