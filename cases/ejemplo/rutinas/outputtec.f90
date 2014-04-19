SUBROUTINE outputtec(time)

!ifort -o outputtec outputtec.f90 /usr/local/tecplot/lib/tecio64.a -lstdc++

USE global_variables
USE geometries

implicit none

integer	:: i,j, ent, itent
real (kind=8)	:: ind, dec, time
logical	:: lexist, lexistT
character(len=256):: filename, filenameT
character(len=10)::number, numbercaso
character(len=1000)::intchar, ncaso
real (kind=8), dimension(:,:),allocatable:: xtec,ytec,ztec,htec,hztec,utec,vtec

! variables to enable writing of TecPlot binary (*.plt) files

  integer (kind = 4)           :: TecIni, TecDat, TecZne
  integer (kind = 4)           :: TecEnd
  integer (kind = 4)           :: VIsDouble = 0
  integer (kind = 4)           :: Debug = 1
  
  integer (kind = 4)           :: III

  character (len = 1)          :: nullchr = char(0)
!--------------------------------------------------------------

if (time==0.0D0) then
qnew_global=qold_global
end if

!Dimensionalization of the variables
do i=1,Nbx; do j=1,Nby
		qreal_global(1,i,j)=qnew_global(1,i,j)*H
		qreal_global(2,i,j)=qnew_global(2,i,j)*U
		qreal_global(3,i,j)=qnew_global(3,i,j)*U
end do; end do



!Writing files
!Write the initial conditions in the first results file
!Escribir tb las condiciones iniciales en un archivo
ind=it/dit
ent=ind
dec=ind-ent


IF (dec==0.0D0.OR.it==0.0D0) THEN
  
  itent=it
  
  write(intchar,*) itent
  intchar=adjustl(intchar)
  number=intchar
  write(ncaso,*) caso
  ncaso=adjustl(ncaso)
  numbercaso=ncaso

!Tiempo
  
  filenameT='Time'//trim(numbercaso)//'.dat'
  inquire(FILE=filenameT, EXIST=lexistT)
  
  IF (.NOT. lexistT) THEN
  open(20,file=filenameT,status='new',action='write')
  ELSE
    
    IF (it==0.0D0) THEN
    open(20,file=filenameT,status='replace',action='write')
    ELSE
    open(20,file=filenameT,status='old',action='write', position='append')
    END IF
    
  END IF

  write(20,180) treal
  180 format(F10.3) !Formato para el vector tiempo
  
  close(20)

!Resultados TECPLOT
  
      !filename='SOL2D.'//trim(number)
      
      write(filename, fmt = '(a,i6.6)')'SOL2D'
      
      I = TecIni('Caso'//trim(ncaso)//NULLCHR,			&   ! title of file
		'X, Y, Z, H, Hz, U, V'//NULLCHR,	&   ! list of variables
		trim(filename)//trim(number)//'.plt'//NULLCHR,	&   ! output file name
		'.'//NULLCHR,						&
		Debug,								&
		VIsDouble)
  

allocate(xtec(Nbx,Nby),ytec(Nbx,Nby),ztec(Nbx,Nby),htec(Nbx,Nby),hztec(Nbx,Nby),utec(Nbx,Nby),vtec(Nbx,Nby))

  DO j=1,Nby; DO i=1,Nbx
    
    if (abs(qreal_global(1,i,j))<=kappa) then
    qreal_global(1,i,j)=0.0D0 
    end if
    if (abs(qreal_global(2,i,j))<=kappa) then
    qreal_global(2,i,j)=0.0D0 
    end if
    if (abs(qreal_global(3,i,j))<=kappa) then
    qreal_global(3,i,j)=0.0D0 
    end if
    xtec(i,j)=x_global(i,j)
    ytec(i,j)=y_global(i,j)
    ztec(i,j)=z_global(i,j)
    htec(i,j)=qreal_global(1,i,j)
    hztec(i,j)=qreal_global(1,i,j)+z_global(i,j)
    utec(i,j)=qreal_global(2,i,j)
    vtec(i,j)=qreal_global(3,i,j)
  END DO; END DO;

	I = TecZne('Zone'//NULLCHR,    &     
			Nbx,			           &
			Nby,			           &
			1, &
			'BLOCK'//NULLCHR,		   &
			NULLCHR//NULLCHR)
			! total number of points
		III = Nby * Nbx
			! write each variable
			! the last argument in the following calls indicates
			! the precision of the the variable,
			! 0 = single, 1 = double
		I   = TecDat(III,xtec,1)
		I   = TecDat(III,ytec,1)
		I   = TecDat(III,ztec,1)
		I   = TecDat(III,htec,1)
		I   = TecDat(III,hztec,1)
		I   = TecDat(III,utec,1)
		I   = TecDat(III,vtec,1)

deallocate(xtec,ytec,ztec,htec,hztec,utec,vtec)
    I   = TecEnd()
close(1)

END IF


END SUBROUTINE outputtec
