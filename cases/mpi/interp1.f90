SUBROUTINE interp1(N,x,y,xn,yn)
!PROGRAM interp1

implicit none

real (kind=8):: xn, yn, maxdo0, minup0, x1, y1, x2, y2
real (kind=8), dimension(N)::x, y, xaux !up0,do0,index1,index2
integer:: N, i, j, k, indexmin, indexmax, sigue

! N=10

!allocate(x(N),y(N),xaux(N))
!PRUEBA:

! x=(/0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,5.0D0,6.0D0,7.0D0,8.0D0,9.0D0/)
! y=(/0.0D0,0.50D0,1.0D0,1.50D0,2.0D0,2.50D0,3.0D0,3.50D0,4.0D0,4.50D0/)
! xn=4.2D0

!Encontrar entre que puntos se encuentra xn
j=1
k=1
minup0=100000000D0
maxdo0=-1000000000D0

IF (xn>=x(1).AND.xn<=x(N)) THEN !Interpolar
  DO i=1,N
    xaux(i)=xn-x(i)

    if (xaux(i)==0.0D0) then
      sigue=0
      yn=y(i)
    else
      if(xaux(i)>0.0D0) then

	if(xaux(i)<minup0) then
	minup0=xaux(i)
	indexmin=i
	end if
	
      else if (xaux(i)<0.0D0) then
	if(xaux(i)>maxdo0) then
	maxdo0=xaux(i)
	indexmax=i
	end if
      end if
      sigue=1
    end if
    if (sigue==0) EXIT
  END DO
    
  IF (sigue==1) THEN
    
    x1=x(indexmin)
    y1=y(indexmin)
      
    x2=x(indexmax)
    y2=y(indexmax)
      
    yn=y1+(xn-x1)*(y1-y2)/(x1-x2)
    
  END IF

ELSE IF (xn<x(1)) THEN !Extrapolar hacia abajo
yn=y(1)+(xn-x(1))*(y(1)-y(2))/(x(1)-x(2))
ELSE IF (xn>x(N)) THEN !Extrapolar hacia afuera
yn=y(N)+(xn-x(N))*(y(N)-y(N-1))/(x(N)-x(N-1))
END IF

! print*,'interp1', xn, yn
! pause




  
!END PROGRAM interp1
END SUBROUTINE interp1
