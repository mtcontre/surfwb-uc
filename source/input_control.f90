subroutine input_control
  use global_variables
  use mpi_surf
  
  implicit none
  
  open(1,file='input.dat')
  read(1,*) caso
  
  !grids
  read(1,*) tinit,tfinal,CFL
  read(1,*) nxi, neta
  read(1,*) batiopt				!solo 0 o 1 por ahora
  if ((batiopt==0).or.(batiopt==1)) then
    read(1,'(A)') batiname(1)
    read(1,'(A)') batiname(2)
    read(1,'(A)') batiname(3)
  elseif ((batiopt==2).or.(batiopt==3)) then
    read(1,'(A)') batiname(1)
  end if
  
  !initial conditions
  read(1,*) initqopt				!solo 0 o 1 por ahora
  if ((initqopt==0).or.(initqopt==1)) then
    read(1,'(A)') initqname(1)
    read(1,'(A)') initqname(2)
    read(1,'(A)') initqname(3)
  elseif ((initqopt==2).or.(initqopt==3)) then
    read(1,'(A)') initqname(1)
  endif
  
  !dimensionalization parameters
  read(1,*) dxi
  read(1,*) deta
  read(1,*) L
  read(1,*) H
  read(1,*) U
  
  !boundary conditions
  read(1,*) CB_real(1)
  read(1,*) CB_real(2)
  read(1,*) CB_real(3)
  read(1,*) CB_real(4)
  
  !numerical method options
  read(1,*) kappa
  read(1,*) rk
  read(1,*) mmopt
  read(1,*) fopt  
  if (fopt==0) then
    Cf = 0
    Coef =0.0d0
  elseif (fopt==1) then
    read(1,*) fM !Si fM=1, un coeficiente de friccion, si fM=2, se usará la matrzi de coeficientes
    read(1,*) Cf !Tipo Fricción: Manning ==1, Chezy ==2, Sampson==3, 0=ninguno
    if (fM==1) then !Un coeficiente para todo el dominio
      read(1,*) Coef !Coeficiente friccion segun sea el caso ya adimensionalizado!!!!
    end if
  end if
  
  !output options
  read(1,*) dit
  if (dit==-1) then
    read(1,*) dtout
  end if
  read(1,*) outopt
  read(1,*) outdir
  close(unit=1)  
  
end subroutine

