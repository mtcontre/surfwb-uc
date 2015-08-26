SUBROUTINE readGA(borde,GA,Ns)

use senales
use global_variables
!use Pich
!Cambiar nombre de archivos si es necesario
implicit none
integer::GA,Ns,borde,i,j

select case(borde)

    CASE(1) !Borde xi=0
    !allocate(qA1(3,Nby),zA1(Nby))
    if (GA==1) then
    allocate(etaL1(Ns,2))
    open(unit=50,file='data/etaL.dat')
    !read(50,*) h01
    h01=0.2D0
    DO i=1,Ns
    read(50,*) etaL1(i,1), etaL1(i,2)
    END DO
    close(unit=50)
    end if
    
    if (GA==2) then
    allocate(qs1(Ns,2),hs1(Ns,2))
    h01=0.41D0
    open(unit=60,file='qs1.dat')
    DO i=1,Ns
    read(60,*) qs1(i,1), qs1(i,2)
    END DO
    close(unit=60)
    open(unit=61,file='hs1.dat')
    DO i=1,Ns
    read(61,*) hs1(i,1), hs1(i,2)
    END DO
    close(unit=61)
    end if
    
    if (GA==3) then
    allocate(us1(Ns,2),hs1(Ns,2))
    h01=0.41D0
    open(unit=70,file='hs1.dat')
    DO i=1,Ns
    read(70,*) hs1(i,1), hs1(i,2)
    END DO
    close(unit=70)
    
    open(unit=71,file='us1.dat')
    DO i=1,Ns
    read(71,*) us1(i,1), us1(i,2)
    END DO
    close(unit=71)
    end if
    
!     if (GA==4) then
!     allocate(hs1(Ns,2))
!     h01=0.9D0
!     open(unit=102,file='Hpich09.dat')
!     DO i=1,Ns
!     read(102,*) hs1(i,1), hs1(i,2)
!     END DO
!     close(unit=102)
!     end if
    
    !Gen Abs Leandro, diferente señal en cada nodo

    if (GA==9) then
 
    allocate(etaL9(Ns,Nby),timeS9(Ns))
    
    !Cambiar nombre para otros casos
    open(unit=90,file='data/etaxi0.dat')
    h01=0.0D0
    read(90,*) ((etaL9(i,j),i=1,Ns),j=1,Nby)
	!etaL9=etaL9+0.001D0
	!etaL9=etaL9*1.5D0
    close(unit=50)
    print*,'fdsa',shape(etaL9)
    print*,'Ns',Ns
    open(unit=91,file='data/timexi0.dat')
    read (91,*) (timeS9(i),i=1,Ns)
    close(unit=91)
    end if
        
!     
    CASE(2) !Borde xi=N
    !allocate(qA2(3,Nby),zA2(Nby))
    if (GA==1) then
    allocate(etaR2(Ns,2))
    open(unit=80,file='etaR.dat')
    read(80,*) h02
    DO i=1,Ns
    read(80,*) etaR2(i,1), etaR2(i,2)
    END DO
    close(unit=80)
    end if
    
    if (GA==2) then
    allocate(qs2(Ns,2),hs2(Ns,2))
    h02=0.33D0
    open(unit=90,file='qs2.dat')
    DO i=1,Ns
    read(90,*) qs2(i,1), qs2(i,2)
    END DO
    close(unit=90)
    
    open(unit=91,file='hs2.dat')
    DO i=1,Ns
    read(91,*) hs2(i,1), hs2(i,2)
    END DO
    close(unit=91)
    end if
    
    if (GA==3) then
    allocate(hs2(Ns,2),us2(Ns,2))
    
    open(unit=100,file='hs2.dat')
    h02=0.33D0
    DO i=1,Ns
    read(100,*) hs2(i,1), hs2(i,2)
    END DO
    close(unit=100)
    
    open(unit=101,file='us2.dat')
    DO i=1,Ns
    read(101,*) us2(i,1), us2(i,2)
    END DO
    close(unit=101)
    end if
    
!     if (GA==4) then
!     allocate(qs2(Ns,2))
!     
!     h02=0.9D0
!     open(unit=102,file='Q5pich.dat')
!     DO i=1,Ns
!     read(102,*) qs2(i,1), qs2(i,2)
!     END DO
!     close(unit=102)
!     
!     end if
!     
!     if (GA==5) then
!     allocate(qs2(Ns,2))
!     
!     h02=0.9D0
!     open(unit=102,file='Q5pich.dat')
!     DO i=1,Ns
!     read(102,*) qs2(i,1), qs2(i,2)
!     END DO
!     close(unit=102)
!        
!     end if
    
    CASE(3) !Borde eta=0
    !allocate(qA3(3,Nbx),zA3(Nbx))
    if (GA==1) then
    allocate(etaL3(Ns,2))
    open(unit=110,file='etaL3.dat')
    read(110,*) h03
    DO i=1,Ns
    read(110,*) etaL3(i,1), etaL3(i,2)
    END DO
    close(unit=110)
    end if
    
    if (GA==2) then
    allocate(qs3(Ns,2),hs3(Ns,2))
    open(unit=111,file='qs3.dat')
    DO i=1,Ns
    read(111,*) qs3(i,1), qs3(i,2)
    END DO
    close(unit=111)
    open(unit=120,file='hs3.dat')
    DO i=1,Ns
    read(120,*) hs3(i,1), hs3(i,2)
    END DO
    close(unit=120)
    end if
    
    if (GA==3) then
    allocate(qs3(Ns,2),us3(Ns,2))
    open(unit=121,file='hs3.dat')
    DO i=1,Ns
    read(121,*) qs3(i,1), qs3(i,2)
    END DO
    close(unit=121)
    open(unit=130,file='us3.dat')
    DO i=1,Ns
    read(130,*) us3(i,1), us3(i,2)
    END DO
    close(unit=130)
    end if
    
    CASE(4) !Borde eta=N
    !allocate(qA4(3,Nbx),zA4(Nbx))
    if (GA==1) then
    allocate(etaR4(Ns,2))
    open(unit=140,file='etaR4.dat')
    read(140,*) h04
    DO i=1,Ns
    read(140,*) etaR4(i,1), etaR4(i,2)
    END DO
    close(unit=140)
    end if
    
    if (GA==2) then
    allocate(qs4(Ns,2),hs4(Ns,2))
    open(unit=150,file='qs4.dat')
    DO i=1,Ns
    read(150,*) qs4(i,1), qs4(i,2)
    END DO
    close(unit=150)
    open(unit=151,file='hs4.dat')
    DO i=1,Ns
    read(151,*) hs4(i,1), hs4(i,2)
    END DO
    close(unit=151)
    end if
    
    if (GA==3) then
    allocate(qs4(Ns,2),us4(Ns,2))
    open(unit=160,file='hs4.dat')
    DO i=1,Ns
    read(160,*) qs4(i,1), qs4(i,2)
    END DO
    close(unit=160)
    open(unit=161,file='us.dat')
    DO i=1,Ns
    read(161,*) us4(i,1), us4(i,2)
    END DO
    close(unit=161)
    end if

END SELECT

END SUBROUTINE readGA
