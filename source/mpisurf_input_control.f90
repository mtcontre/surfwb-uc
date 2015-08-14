subroutine input_control

use global_variables
use senales
use custombc
USE MPI_SURF

implicit none

integer :: omp_nthreads, omp_threadid, omp_get_num_threads, omp_get_thread_num
  
!get input data directory from environment
call get_environment_variable('INDIR',indir)

!start reading input data
open(1,file=trim(indir)//'/input.dat')
read(1,*) caso  
read(1,*) tinit
treal=tinit
read(1,*) tfinal
read(1,*) CFL
read(1,*) Nbx
read(1,*) Nby
read(1,*) batiopt
read(1,'(A)') batiname(1)!Xmatrix
read(1,'(A)') batiname(2)!Ymatrix
read(1,'(A)') batiname(3)!Zmatrix
read(1,*) initqopt
read(1,'(A)') initqname(1)!h-matrix
read(1,'(A)') initqname(2)!u-matrix
read(1,'(A)') initqname(3)!v-matrix
read(1,*) dxi
read(1,*) deta
read(1,*) L
read(1,*) H
read(1,*) U
t=treal*U/L
read(1,*) CB(1)	!Boundary Condition for xi_0,j
read(1,*) CB(2)	!Boundary Condition for xi=nbx
read(1,*) CB(3) !Boundary Condition for eta_i,0
read(1,*) CB(4) !Boundary Condition for eta_i,Nby
read(1,*) dit	!dit to print results, write files every dit  iterations
if (dit == -1) then	
  read(1,*) dtout	!or every dtout seconds
end if
read(1,*) kappa 	!To consider a 0.0... value where is it used?
read(1,*) rk 		!Runge Kutta method 1=Rk4, 2=Rk2
read(1,*) mmopt 	!1=Minmod, 2=Superbee Limiters
read(1,*) fopt 		!Con o sin friccion, No == 0, Si ==1
IF (fopt==0) THEN
 Cf=0
 Coef=0.0D0
END IF
IF (fopt==1) THEN
  !Si fM=1, un coeficiente de friccion
  !si fM=2, se usará la matriz de coeficientes
  read(1,*) fM 
  read(1,*) Cf !Tipo Fricción: Manning ==1, Chezy ==2, Sampson==3, 0=ninguno
    if (fM==1) then !Un coeficiente para todo el dominio
    read(1,*) Coef !Coeficiente friccion segun sea el caso ya adimensionalizado!!!!
    end if
END IF
read(1,*) outopt !1=Matlab, 2=Tecplot files
read(1,*) outdir
print*,'-------------',outdir
close(1)

if (myrank==0) then
  write(*,100) caso
  100 FORMAT ('Caso: ', T25, I4)

  write(*,110) tfinal
  110 FORMAT ('T_final (s): ', T25, F8.2)

  write(*,120) CFL
  120 FORMAT ('CFL: ', T25, F5.2)

  write(*,130) Nbx
  130 FORMAT ('Nbx: ', T25,  I6)

  write(*,140) Nby
  140 FORMAT ('Nby: ', T25, I6)


  write(*,150)  L, H, U
  150 FORMAT ('Escalas (L, H, U): ', T25, 3F5.2)

  write(*,160) CB(1), CB(2), CB(3), CB(4)
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
end if

end subroutine input_control
