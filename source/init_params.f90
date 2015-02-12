subroutine init_params
  !initialize parameters like t,treal, or boundaries parameters
  
  use mpi
  use mpi_surf
  use global_variables
  
  implicit none
  
  logical :: dir_exists
  character(256) :: command, intchar
  
  !time
  treal = tinit
  t = treal*U/L
  
  !adimensional
  FR2 = U*U/(g*H)
  
  !outputdirectory
  write(intchar,*) nproc
  outdir=trim(outdir)//trim(adjustl(intchar))//'/'  
  if (myrank == master .and. .not. dir_exists) then
    
    command='mkdir '//trim(outdir)
    call system(trim(command))
    command='mkdir '//trim(outdir)//'/grids'
    call system(trim(command))
  end if
  
 
  !print to screen
  if (myrank==master) then
    write(*,100) caso
    100 FORMAT ('Caso: ', T25, I4)  
    write(*,'("T_init (s)",T25,F8.2)') tinit  
    
    write(*,110) tfinal
    110 FORMAT ('T_final (s): ', T25, F8.2)

    write(*,120) CFL
    120 FORMAT ('CFL: ', T25, F5.2)

    write(*,130) Nxi
    130 FORMAT ('Nbx: ', T25,  I6)

    write(*,140) Neta
    140 FORMAT ('Nby: ', T25, I6)


    write(*,150)  L, H, U
    150 FORMAT ('Escalas (L, H, U): ', T25, 3F5.2)

    write(*,160) CB_real(1), CB_real(2), CB_real(3), CB_real(4)
    160 FORMAT ('Condiciones de Borde: ', T25, 4I2)

    if (mmopt==1) then
    write(*,161) mmopt
    161 FORMAT ('Usando limitador MINMOD, mmopt= ', T30, 4I2)
    else if (mmopt==2) then
    write(*,162) mmopt
    162 FORMAT ('Usando limitador SUPERBEE, mmopt= ', T30, 4I2)
    else
    write(*,163) mmopt
    163 FORMAT ('Usando limitador MC, mmopt= ', T30, 4I2)
    end if
    
    print*, 'Fr2= ', FR2
  end if
end subroutine