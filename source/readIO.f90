SUBROUTINE readIO(borde,IO,Ns)
!Rutina que carga las señales de entrada o salida en alguno de los bordes para usar Outflow or Inflow

use senales
use global_variables

!Cambiar nombre de archivos si es necesario

implicit none
integer::IO,Ns,borde,i,j

select case(borde)

    CASE(1) !Borde xi=0
   
    if (IO==1) then !Inflow
    allocate(qsx1(Nby,Ns),qsy1(Nby,Ns),timeS1(Ns))
    allocate(etas1(Nby,Ns))
    open(unit=50,file='qsx1.dat',form='unformatted') !En binario
    read (unit=50) ((qsx1(i,j),i=1,Nby),j=1,Ns)
    close(unit=50)
    
     
    open(unit=51,file='qsy1.dat',form='unformatted')
    read (unit=51) ((qsy1(i,j),i=1,Nby),j=1,Ns)
    close(unit=51)

    open(unit=51,file='etas1.dat',form='unformatted')
    read (unit=51) ((etas1(i,j),i=1,Nby),j=1,Ns)
    close(unit=51)

    open(unit=52,file='timeS1.dat',form='unformatted') !En binario
    read (unit=52) (timeS1(i),i=1,Ns)
    close(unit=52)
    end if
       
    if (IO==2) then !Outflow h fija
    allocate(etas1(Nby,Ns),timeS1(Ns))

    open(unit=60,file='etas1H.dat',form='unformatted')
    read (unit=60) ((etas1(i,j),i=1,Nby),j=1,Ns)
    close(unit=60)
        
    open(unit=53,file='timeS1.dat',form='unformatted') !En binario
    read (unit=53) (timeS1(i),i=1,Ns)
    close(unit=53)
    
    end if
    
    CASE(2) !Borde xi=Nbx
   
    if (IO==1) then !Inflow

    allocate(qsx2(Nby,Ns),qsy2(Nby,Ns),etas2(Nby,Ns),timeS2(Ns))
    

    
    open(unit=60,file='qsx2.dat',form='unformatted')
    read (unit=60) ((qsx2(i,j),i=1,Nby),j=1,Ns)
    close(unit=60)

    open(unit=61,file='qsy2.dat',form='unformatted')
    read (unit=61) ((qsy2(i,j),i=1,Nby),j=1,Ns)
    close(unit=61)

     open(unit=611,file='etas2.dat',form='unformatted')
    read (unit=611) ((etas2(i,j),i=1,Nby),j=1,Ns)
    close(unit=611)

    open(unit=62,file='timeS2.dat',form='unformatted')
    read (unit=62) (timeS2(i),i=1,Ns)
    close(unit=62)

    
    end if
    
    
    
    if (IO==2) then !Outflow h fija
    allocate(hs2(Ns,2))
    open(unit=80,file='hs2.dat')
    DO i=1,Ns
    read(80,*) hs2(i,1), hs2(i,2)
    END DO
    close(unit=80)
    end if
    
    if (IO==3) then !Outflow U fija (Solo para U fija, V=0)
    allocate(us2(Ns,2))
    open(unit=81,file='us2.dat')
    DO i=1,Ns
    read(81,*) us2(i,1), us2(i,2)
    END DO
    close(unit=81)
    end if
    
    
    CASE(3) !Borde eta=0
   
    if (IO==1) then !Inflow
    allocate(qsx3(Nbx,Ns),qsy3(Nbx,Ns),timeS3(Ns))
    allocate(etas3(Nbx,Ns))
    open(unit=80,file='qsx3.dat',form='unformatted')
    read (unit=80) ((qsx3(i,j),i=1,Nbx),j=1,Ns)
    close(unit=80)
  
    open(unit=81,file='qsy3.dat',form='unformatted')
    read (unit=81) ((qsy3(i,j),i=1,Nbx),j=1,Ns)
    close(unit=81)
    
    open(unit=81,file='etas3.dat',form='unformatted')
    read (unit=81) ((etas3(i,j),i=1,Nbx),j=1,Ns)
    close(unit=81)
    
    open(unit=82,file='timeS3.dat',form='unformatted')
    read (unit=82) (timeS3(i),i=1,Ns)
    close(unit=82)
   
    end if
    
    if (IO==2) then !Outflow h fija
    allocate(hs3(Ns,2))
    open(unit=100,file='hs3.dat')
    DO i=1,Ns
    read(100,*) hs3(i,1), hs3(i,2)
    END DO
    close(unit=100)
    end if
    
    CASE(4) !Borde eta=Nby
   
    if (IO==1) then !Inflow
    allocate(qsx4(Nbx,Ns),qsy4(Nbx,Ns),timeS4(Ns))
    allocate(etas4(Nbx,Ns))
    open(unit=80,file='qsx4.dat',form='unformatted')
    read (unit=80) ((qsx4(i,j),i=1,Nbx),j=1,Ns)
    close(unit=80)
    
    open(unit=81,file='qsy4.dat',form='unformatted')
    read (unit=81) ((qsy4(i,j),i=1,Nbx),j=1,Ns)
    close(unit=81)
    
    open(unit=81,file='etas4.dat',form='unformatted')
    read (unit=81) ((etas4(i,j),i=1,Nbx),j=1,Ns)
    close(unit=81)
    
    open(unit=82,file='timeS4.dat',form='unformatted')
    read (unit=82) (timeS4(i),i=1,Ns)
    close(unit=82)
    end if
    
    if (IO==2) then !Outflow h fija
    allocate(hs4(Ns,2))
    open(unit=120,file='hs4.dat')
    DO i=1,Ns
    read(120,*) hs4(i,1), hs4(i,2)
    END DO
    close(unit=120)
    end if

end Select


END SUBROUTINE readIO