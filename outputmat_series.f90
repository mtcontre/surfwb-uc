SUBROUTINE outputmat(time)

USE global_variables
USE geometries
USE LagrangianPT
USE TimeSeries
implicit none

integer	:: i,j, ent, itent, ent1
real (kind=8)	:: ind, dec, time, dit1
real (kind=8)	:: ind1, dec1
logical	:: lexist, lexistT
 character(len=1100):: filename, filenameT
 character(len=10)::number, numbercaso
 character(len=1000)::intchar, ncaso
 character(len=1100):: filename_LPTx, filename_LPTy, filename_LPTu, filename_LPTv
!-------------------------
logical	:: lexistPh, lexistPu, lexistPv, lexistT_eta, lexist_LPTx, lexist_LPTy, lexist_LPTu, lexist_LPTv 
 character(len=1100):: filenamePh,filenamePu,filenamePv,filenameT_eta, path
logical	:: lexist_param
 character(len=1100):: filename_param
logical	:: lexist_TSxy
 character(len=1100):: filename_TSxy
logical	:: lexist_H, lexist_U, lexist_V
 character(len=1100):: filename_H, filename_U, filename_V
!! Blend Interpolation for times series at a point
!real (kind=8), dimension(135)::x0, y0, i0, j0, H0
!real (kind=8), dimension(1)::m1,m2
!real (kind=8), dimension(135)::r, s, x, x00, x01, x10, x11 ! Interpolation by blending (see blend.f90, function blend_102)
!! Coordenadas Puntos
! 
do i=1,Nts
		m1=minloc((x_global-x0(i))*(x_global-x0(i))+(y_global-y0(i))*(y_global-x0(i)),1, &
			  mask=(x_global-x0(i)).ge.0 .and. &
			       (y_global-y0(i)).ge.0 )
		m2=minloc((x_global-x0(i))*(x_global-x0(i))+(y_global-y0(i))*(y_global-x0(i)),2, &
			  mask=(x_global-x0(i)).ge.0 .and. &
			       (y_global-y0(i)).ge.0 )		
			       		!m=x0(i)
		m1=minloc(x_global(:,1)-x0(i),1,mask=(x_global(:,1)-x0(i)).ge.0)
		i0(i)=m1(1)-1
		r(i)=(x0(i)-x_global(i0(i),1))/(x_global(2,1)-x_global(1,1))		
		m2=minloc(y_global(1,:)-y0(i),1,mask=(y_global(1,:)-y0(i)).ge.0)
		j0(i)=m2(1)-1
		s(i)=(y0(i)-y_global(1,j0(i)))/(y_global(1,2)-y_global(1,1))		
		x00(1,i)=qreal_global(1,i0(i),j0(i))    +z_global(i0(i),j0(i))
		x10(1,i)=qreal_global(1,i0(i)+1,j0(i))  +z_global(i0(i)+1,j0(i))
		x01(1,i)=qreal_global(1,i0(i),j0(i)+1)  +z_global(i0(i),j0(i)+1)
		x11(1,i)=qreal_global(1,i0(i)+1,j0(i)+1)+z_global(i0(i)+1,j0(i)+1)
		x00(2,i)=qreal_global(2,i0(i),j0(i))
		x10(2,i)=qreal_global(2,i0(i)+1,j0(i))
		x01(2,i)=qreal_global(2,i0(i),j0(i)+1)
		x11(2,i)=qreal_global(2,i0(i)+1,j0(i)+1)
		x00(3,i)=qreal_global(3,i0(i),j0(i))    
		x10(3,i)=qreal_global(3,i0(i)+1,j0(i))  
		x01(3,i)=qreal_global(3,i0(i),j0(i)+1)  
		x11(3,i)=qreal_global(3,i0(i)+1,j0(i)+1)
		!x00(i)=qreal_global(1,i0(i),j0(i))    +z_global(i0(i),j0(i))
		!x10(i)=qreal_global(1,i0(i)+1,j0(i))  +z_global(i0(i)+1,j0(i))
		!x01(i)=qreal_global(1,i0(i),j0(i)+1)  +z_global(i0(i),j0(i)+1)
		!x11(i)=qreal_global(1,i0(i)+1,j0(i)+1)+z_global(i0(i)+1,j0(i)+1)
		
		!print *, '----------------------------'
		!print *, 'r i0 ', r(i) , i0(i)
		!print *, 's j0', s(i) , j0(i)
		!print *, 'x0 xg ', x0(:,i), x_global(i0(i),1)		
		!print *, 'y0 yg', y0(i), y_global(1,j0(i))
		!!print *, 'y0-yg', y_global(1,j0(i))-y0(i)
		!print *,'x00',x00(:,i)
		!print *,'qreal',qreal_global(:,i0(i),j0(i))
		!print *,'z_g',z_global(i0(i),j0(i))
		!print *, 'iteration', i
		
end do
!pause
!print *, 'r s', r , s

!		print *,'x00',x00(:,1)
!		print *,'qreal',qreal_global(:,i0(1),j0(1))
!		print *,'z_g',z_global(i0(1),j0(1))
!------------------------
if (time==0.0D0) then
qnew_global=qold_global
!~ print *, 'x0 y0', x0 , y0
!~ print *, 'x1 y1', x_global(i0,1), y_global(1,j0) 
!~ write(*,'("Paused, Pres ENTER to Continue")') !BP
!~ read(*,*) 
end if

!Dimensionalization of the variables
do i=1,Nbx; do j=1,Nby
		qreal_global(1,i,j)=qnew_global(1,i,j)*H
		qreal_global(2,i,j)=qnew_global(2,i,j)*U
		qreal_global(3,i,j)=qnew_global(3,i,j)*U
end do; end do

do i=1,Nts
		!H0(i)=qreal_global(1,i0(i),j0(i))+z_global(i0(i),j0(i))
		! function blend_102 by John Burkardt (see blend.f90)
		H0(i)=            	          + x00(1,i) &
				+ r(i) *        ( - x00(1,i) + x10(1,i) ) & 
				+ s(i) *        ( - x00(1,i)              + x01(1,i) ) &
				+ r(i) * s(i) * ( + x00(1,i) - x10(1,i)   - x01(1,i) + x11(1,i) )

		Uts(i)=            	          + x00(2,i) &
				+ r(i) *        ( - x00(2,i) + x10(2,i) ) & 
				+ s(i) *        ( - x00(2,i)              + x01(2,i) ) &
				+ r(i) * s(i) * ( + x00(2,i) - x10(2,i)   - x01(2,i) + x11(2,i) )

		Vts(i)=            	          + x00(3,i) &
				+ r(i) *        ( - x00(3,i) + x10(3,i) ) & 
				+ s(i) *        ( - x00(3,i)              + x01(3,i) ) &
				+ r(i) * s(i) * ( + x00(3,i) - x10(3,i)   - x01(3,i) + x11(3,i) )
			
end do

!print *, 'r s', r(3) , s(3)
!print *, 'x00 x10', x00(3) , x10(3)
!print *, 'x01 x11', x01(3) , x11(3)
!pause

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
  path='results/'
  filename='results/SOL2D.'//trim(number)//'.dat'
  filenameT='results/Time'//trim(numbercaso)//'.dat'
  inquire(FILE=filename, EXIST=lexist)
  inquire(FILE=filenameT, EXIST=lexistT) 

  IF (.NOT. lexist) THEN
  open(10,file=filename,status='new',action='write')
  ELSE
  open(10,file=filename,status='replace',action='write')
  END IF
  
  IF (.NOT. lexistT) THEN
  open(20,file=filenameT,status='new',action='write')
  ELSE
    
    IF (it==0.0D0) THEN
    open(20,file=filenameT,status='replace',action='write')
    ELSE
    open(20,file=filenameT,status='old',action='write', position='append')
    END IF
    
  END IF
!Escribo
!Tiempo
  write(20,180) treal
  180 format(F10.3) !Formato para el vector tiempo
!Resultados

  DO j=1,Nby; DO i=1,Nbx
!  DO j=1,Nby,4; DO i=1,Nbx,4
    
    if (abs(qreal_global(1,i,j))<=kappa) then
    qreal_global(1,i,j)=0.0D0 
    end if
    if (abs(qreal_global(2,i,j))<=kappa) then
    qreal_global(2,i,j)=0.0D0 
    end if
    if (abs(qreal_global(3,i,j))<=kappa) then
    qreal_global(3,i,j)=0.0D0 
    end if
    
    write(10,170) x_global(i,j),y_global(i,j),z_global(i,j),qreal_global(1,i,j), qreal_global(2,i,j),qreal_global(3,i,j)
    170 format(F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5 / )
    
  END DO; END DO;

	

  close(10)
  close(20)
  
!   call system('gzip -f '//filename)

END IF


!--------------------------------------------------------- PARAMETRES
filename_param=trim(path)//'param.dat'  
inquire(FILE=filename_param, EXIST=lexist_param)
	IF (.NOT. lexist_param) THEN
  		open(25,file=filename_param,status='new',action='write')
  	ELSE
  		open(25,file=filename_param,status='replace',action='write')
  	END IF

!write(25,250) caso, Nbx, Nby, L, H, U, Nts, num
    	!250 format(I5, TR1, I5, TR1, I5, TR1, F15.5, TR1, F15.5, TR1, F15.5, TR1, I5, TR1, I5/ )
	write(25,250) caso, Nbx, Nby, Nts, num, it, dit
    	250 format(I5, TR1, I5, TR1, I5, TR1, I5, TR1, I5, TR1, F15.5, TR1, I10/ )

	!print *,'caso', caso, Nbx, Nby, L, H, U, Nts, num
!pause

 close(25)
!Scope falso!!!
if (1==0) then


filenamePh='results/timeseries/Etah_'//trim(number)//'.dat'
filenamePu='results/timeseries/Etau_'//trim(number)//'.dat'
filenamePv='results/timeseries/Etav_'//trim(number)//'.dat'
inquire(FILE=filenamePh, EXIST=lexistPh)
inquire(FILE=filenamePu, EXIST=lexistPu)
inquire(FILE=filenamePv, EXIST=lexistPv)

	  IF (.NOT. lexistPh) THEN  
		open(30,file=filenamePh,status='new',action='write',access='stream')  
	  ELSE
	  
	      IF (dec==0.0D0.OR.it==0.0D0) THEN
			    open(30,file=filenamePh,status='replace',action='write',access='stream')
		ELSE
			    open(30,file=filenamePh,status='old',action='write', position='append',access='stream')
	      END IF

	  END IF
	  
	  IF (.NOT. lexistPu) THEN  
		open(31,file=filenamePu,status='new',action='write',access='stream')  
	  ELSE
	  
	      IF (dec==0.0D0.OR.it==0.0D0) THEN
			    open(31,file=filenamePu,status='replace',action='write',access='stream')
	      ELSE
			    open(31,file=filenamePu,status='old',action='write', position='append',access='stream')
	      END IF

	  END IF
	  IF (.NOT. lexistPv) THEN  
		open(32,file=filenamePv,status='new',action='write',access='stream')  
	  ELSE
	  
	      IF (dec==0.0D0.OR.it==0.0D0) THEN
			    open(32,file=filenamePv,status='replace',action='write',access='stream')
	      ELSE
			    open(32,file=filenamePv,status='old',action='write', position='append',access='stream')
	      END IF

	  END IF
	  
print*, H0
write(30,111) H0
write(31,111) Uts
write(32,111) Vts
111 format(100000(F15.5) / )

	  
	  
end if
END SUBROUTINE outputmat



