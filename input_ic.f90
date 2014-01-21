SUBROUTINE input_ic
  use global_variables
  use multigrid_surf
  use mpi_surf
  implicit none

  integer :: i,j,level
  real (kind=8), dimension(:), allocatable :: hin, uin, vin
  allocate(initq(ngrids))
  

  do level=1,ngrids
    allocate(initq(level)%H(nxi(level),neta(level)), &
	initq(level)%U(nxi(level),neta(level)),&
	initq(level)%V(nxi(level),neta(level)))
    SELECT CASE (int(initqopt(level)))
      CASE(0)
	open(unit=2,file=initqname(level, 1),form='unformatted')
	read(2) ((initq(level)%H(i,j),j=1,neta(level)), i=1,nxi(level))
	close(unit=2)
	
	open(unit=2,file=initqname(level,2),form='unformatted')
	read(2) ((initq(level)%U(i,j), j=1,neta(level)), i=1,nxi(level))
	close(unit=2)
	
	open(unit=2,file=initqname(level,3),form='unformatted')
	read(2) ((initq(level)%V(i,j),j=1,neta(level)),i=1,nxi(level))
	close(unit=2)
      CASE(1)
	open(unit=2,file=initqname(level,1))
	read(2,*) ((initq(level)%H(i,j),j=1,neta(level)),i=1,nxi(level))
	close(unit=2)
	
	open(unit=2,file=initqname(level,2))
	read(2,*) ((initq(level)%U(i,j),j=1,neta(level)),i=1,nxi(level))
	close(unit=2)
	
	open(unit=2,file=initqname(level,3))
	read(2,*) ((initq(level)%V(i,j),j=1,neta(level)),i=1,nxi(level))
	close(unit=2)
      CASE (2)
	allocate(hin(nxi(level)*neta(level)),&
	  uin(nxi(level)*neta(level)), vin(nxi(level)*neta(level)))
	open(unit=99,file=initqname(level,1))
	Do i=1,nxi(level)*neta(level)
	  read(99,*) hin(i), uin(i), vin(i)
	End Do
	close(unit=99)
	do i=1,nxi(level),1; do j=1,neta(level),1
	  initq(level)%H(i,j)=hin(j+(i-1)*neta(level))
	  initq(level)%U(i,j)=uin(j+(i-1)*neta(level))
	  initq(level)%V(i,j)=vin(j+(i-1)*neta(level))
	end do; end do
	deallocate(hin,uin,vin)  
      CASE (3)
	allocate(hin(nxi(level)*neta(level)),&
	  uin(nxi(level)*neta(level)), vin(nxi(level)*neta(level)))
	open(unit=99,file=initqname(level,1))
	Do i=1,nxi(level)*neta(level)
	  read(99,*) hin(i), uin(i), vin(i)
	End Do
	close(unit=99)
	do i=1,neta(level),1; do j=1,nxi(level),1
	  initq(level)%H(i,j)=hin(i+(j-1)*nxi(level))
	  initq(level)%U(i,j)=uin(i+(j-1)*nxi(level))
	  initq(level)%V(i,j)=vin(i+(j-1)*nxi(level))
	end do; end do
	deallocate(hin,uin,vin)

    END SELECT
  end do
  !for debugging..show first 3row/columns
  level=1
  print*,'h------------'
  do i =1,3
    print*,initq(level)%H(i,1),initq(level)%H(i,2),initq(level)%H(i,3)
  end do
!   pause
  
  print*,'u------------'
  do i =1,3
    print*,initq(level)%U(i,1),initq(level)%U(i,2),initq(level)%U(i,3)
  end do
!   pause
  
  print*,'v------------'
  do i =1,3
    print*,initq(level)%V(i,1),initq(level)%V(i,2),initq(level)%V(i,3)
  end do

END SUBROUTINE input_ic