subroutine input_control
  use global_variables
  use senales
  use couplingbc
  use multigrid_surf
  use mpi_surf
  implicit none
  integer ::i
  character(len=256) ::command,intchar
  logical::dir_exists

  !start reading input data
  open(1,file=trim(indir)//'/input.dat')

  read(1,*) caso
  read(1,*) tinit
  read(1,*) tfinal
  read(1,*) CFL	  
  ngrids=1			
  !number of overlapped grids (1fornow)
  allocate(Nxi(ngrids),Neta(ngrids),&
    batinames(ngrids,3),initqnames(ngrids,3),batiopts(ngrids),initqopts(ngrids)) 
  do i=1,ngrids
    read(1,*) Nxi(i)
    read(1,*) Neta(i)
    read(1,*) batiopts(i)	
    if ((batiopts(i)==0).or.(batiopts(i)==1)) then	
      read(1,'(A)') batinames(i,1)!Xmatrix	
      read(1,'(A)') batinames(i,2)!Ymatrix
      read(1,'(A)') batinames(i,3)!Zmatrix
    else if ((batiopts(i)==2).or.(batiopts(i)==3)) then
      read(1,'(A)') batinames(i,1)!3columns X,Y,Z
    end if
    
    read(1,*) initqopts(i)
    if ((initqopts(i)==0).or.(initqopts(i)==1)) then
      read(1,'(A)') initqnames(i,1)!h-matrix
      read(1,'(A)') initqnames(i,2)!u-matrix
      read(1,'(A)') initqnames(i,3)!v-matrix
    else if ((initqopts(i)==2).or.(initqopts(i)==3)) then
      read(1,'(A)') initqnames(i,1)!recolumns, h,u,v
    end if
  end do
  read(1,*) dxi
  read(1,*) deta
  read(1,*) L
  read(1,*) H
  read(1,*) U  
  !now i can assign time
  treal=tinit
  !adimensionalize
  t=treal*U/L
  read(1,*) CB_real(1)	!Boundary Condition for xi_0,j

  !*Boundary conditions for the coarsest level
  !Diferentes tipos de CB
  !Boundary Condition for xi=0
  if (CB_real(1)==0) then 
  !user-defined values for boundary ghost cells
    read(1,*) fnamexi0(1)
    read(1,*) nt_xi0g1
    read(1,*) dt_xi0g1
    read(1,*) optxi0g1
    read(1,*) fnamexi0(2)
    read(1,*) nt_xi0g2
    read(1,*) dt_xi0g2
    read(1,*) optxi0g2
    flagxi0=1
    call couplingbc_xi0    
  else
    flagxi0=0
  end if
    
  if (CB_real(1)==4) then !GENABS 1, 2, 3 o 9
    !Cienfuegos Generation-Absorption boundary condition
    read(1,*) GA1
    read(1,*) Nsenal1
  end if
      
  !Boundary Condition for xi=nbx
  read(1,*) CB_real(2)	
  if (CB_real(2)==0) then 
    !Customized values for boundary ghost cells
    read(1,*) fnamexiN(1)
    read(1,*) nt_xiNg1
    read(1,*) dt_xiNg1
    read(1,*) optxiNg1
    read(1,*) fnamexiN(2)
    read(1,*) nt_xiNg2
    read(1,*) dt_xiNg2
    read(1,*) optxiNg2
    call couplingbc_xiN
    flagxiN=1
  else
    flagxiN=0
  end if  
  IF (CB_real(2)==4) THEN
    read(1,*) GA2
    read(1,*) Nsenal2
    call readGA(2,GA2,Nsenal2)
  END IF
   
  !Boundary Condition for eta_i,0
  read(1,*) CB_real(3)
  if (CB_real(3)==0) then 
    !Customized values for boundary ghost cells
    read(1,'(A)') fnameeta0(1)
    read(1,*) nt_eta0g1
    read(1,*) dt_eta0g1
    read(1,*) opteta0g1
    read(1,'(A)') fnameeta0(2)
    read(1,*) nt_eta0g2
    read(1,*) dt_eta0g2
    read(1,*) opteta0g2    
    call couplingbc_eta0
    flageta0=1
  else
    flageta0=0
  end if
  IF (CB_real(3)==4) THEN
    read(1,*) GA3
    read(1,*) Nsenal3
  call readGA(3,GA3,Nsenal3)
  END IF

  !Boundary Condition for eta_i,Nby
  read(1,*) CB_real(4)
  if (CB_real(4)==0) then 
    !Customized values for boundary ghost cells
    read(1,'(A)') fnameetaN(1)
    read(1,*) nt_etaNg1
    read(1,*) dt_etaNg1
    read(1,*) optetaNg1
    read(1,'(A)') fnameetaN(2)
    read(1,*) nt_etaNg2
    read(1,*) dt_etaNg2
    read(1,*) optetaNg2
    print*,trim(fnameetaN(1)),trim(fnameetaN(2))
    call couplingbc_etaN
    flagetaN=1
  else
    flagetaN=0
  end if
  IF (CB_real(4)==4) THEN
    read(1,*) GA4
    read(1,*) Nsenal4
    call readGA(4,GA4,Nsenal4)
  END IF
  read(1,*) dit	!dit to print results, write files every pdt iterations
  if (dit == -1) then
    read(1,*) dtout
  end if
  read(1,*) kappa !To consider a 0.0 value
  read(1,*) rk !Runge Kutta method 1=Rk4, 2=Rk2
  read(1,*) mmopt !1=Minmod, 2=Superbee Limiters
  read(1,*) fopt !Con o sin friccion, No == 0, Si ==1
  !for now, same choice with every grid
  IF (fopt==0) THEN
    Cf=0
    Coef=0.0D0
  END IF
  IF (fopt==1) THEN
    read(1,*) fM !Si fM=1, un coeficiente de friccion, si fM=2, se usará la matrzi de coeficientes
    read(1,*) Cf !Tipo Fricción: Manning ==1, Chezy ==2, Sampson==3, 0=ninguno
    if (fM==1) then !Un coeficiente para todo el dominio
    read(1,*) Coef !Coeficiente friccion segun sea el caso ya adimensionalizado!!!!
    end if
  END IF
  read(1,*) outopt !1=Matlab, 2=Tecplot files
  read(1,*) outdir
  
  !create outdir if necessary
  inquire(file=outdir,exist=dir_exists)
  if (.not. dir_exists) then
    command='mkdir '//trim(outdir)
    call system(trim(command))
    command='mkdir '//trim(outdir)//'/grids'
    call system(trim(command))
    command='mkdir '//trim(outdir)//'/timeseries'
    call system(trim(command))
  end if
  close(1) 
  
  
  FR2=U**2.0D0/(g*H)
  
  !this print should be somewhere else
  write(*,100) caso
  100 FORMAT ('Caso: ', T25, I4)
  
  write(*,'("T_init (s)",T25,F8.2)') tinit
  
  write(*,110) tfinal
  110 FORMAT ('T_final (s): ', T25, F8.2)

  write(*,120) CFL
  120 FORMAT ('CFL: ', T25, F5.2)

  write(*,130) Nxi(1)
  130 FORMAT ('Nbx: ', T25,  I6)

  write(*,140) Neta(1)
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

end subroutine input_control

