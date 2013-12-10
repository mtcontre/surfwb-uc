Program Convert

!Programa que convierte resultados .dat -ASCII a .plt (tecplot)
!ifort -o convert Convert.f90 /usr/local/tecplot/lib/tecio64.a -lstdc++

implicit none

integer	:: i,j, ent, itent, it, Nx, Ny, caso
real (kind=8), dimension(:,:),allocatable:: x, y, z, hz
real (kind=8), dimension(:,:,:),allocatable:: q
logical	:: lexist
character(len=1100):: filename
character(len=10)::number
character(len=1000)::intchar, ncaso
! variables to enable writing of TecPlot binary (*.plt) files

  integer (kind = 4)           :: TecIni, TecDat, TecZne
  integer (kind = 4)           :: TecEnd
  integer (kind = 4)           :: VIsDouble = 0
  integer (kind = 4)           :: Debug = 1
  
  integer (kind = 4)           :: III

  character (len = 1)          :: nullchr = char(0)
!--------------------------------------------------------------

!Writing files
!Write the initial conditions in the first results file
!Escribir tb las condiciones iniciales en un archivo
ind=it/dit
ent=ind
dec=ind-ent

Nf=4250  !Numero final archivos SOL2D.Nf.dat
dn=10 	 !Intervalo archivos, cada cuantas iteraciones fue grabado un archivo
it=0
Nx=201
Ny=41
 caso=24
allocate(x(Nx,Ny),y(Nx,Ny),z(Nx,Ny),q(3,Nx,Ny),hz(i,j))


write(ncaso,*) caso
ncaso=adjustl(ncaso)
  
  
Do i=0,dn,Nf

  itent=it
  write(intchar,*) itent
  intchar=adjustl(intchar)
  number=intchar
  
  filename='SOL2D.'//trim(number)//'.dat'
  !filename='/Resultados/Case'//trim(numbercaso)//'/SOL2D.'//trim(number)//'.dat'
  inquire(FILE=filename, EXIST=lexist)
  
  IF (.NOT. lexist) THEN
  !open(10,file=filename,status='new',action='write')
  print*, 'Error, file does not exist'
  stop
  
  ELSE
  open(10,file=filename)
  !Leo Resultados
  Do i=1,Nx;  Do j=1,Ny
  
  read(10,170) x(i,j),y(i,j),z(i,j),q(1,i,j), q(2,i,j),q(3,i,j)
  170 format(F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5 / )
  
  hz(i,j)=z(i,j)+q(1,i,j)
  
  End Do; End Do;
  close(10)
  END IF
  
!Resultados TECPLOT

      !filename='SOL2D.'//trim(number)
      
      write(filename, fmt = '(a,i6.6)')'SOL2Dtec.'
      
      I = TecIni('Caso'//trim(ncaso)//NULLCHR,			&   ! title of file
		'X, Y, Z, H, Hz, U, V'//NULLCHR,	&   ! list of variables
		trim(filename)//trim(number)//'.plt'//NULLCHR,	&   ! output file name
		'.'//NULLCHR,						&
		Debug,								&
		VIsDouble)
  

   
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
		I   = TecDat(III,x,1)
		I   = TecDat(III,y,1)
		I   = TecDat(III,z,1)
		I   = TecDat(III,q(1,:,:),1)
		I   = TecDat(III,hz,1)
		I   = TecDat(III,q(2,:,:),1)
		I   = TecDat(III,q(3,:,:),1)

    I   = TecEnd()
close(1)

End Do




End Program Convert