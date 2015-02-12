subroutine lpt_p() 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Program to perform the
! Lagrangian Particle Tracking of sediment particles
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

USE global_variables
USE LagrangianPT


implicit none

  ! local variables
real (kind = 8) :: dtp
!~ character (len = 256) :: filename

! BEGIN!

dtp = dt  / 2.0E+00  ! particle time-step
!~ print *, 'dtp lpt_p', dtp

! allocate geometry and flow variables

  q_LPT(1,:,:) = qreal_global(2,:,:)
  q_LPT(2,:,:) = qreal_global(3,:,:)
  q_LPT(3,:,:) = zero
  treeflag = 0
! print *, 'LPT_P q ', qreal_global(2,1:5,1:5)!, q_LPT(:,1,1)

!-------- midpoint predictor ---------

	do n_LPT = 1,num
                hitflag = 0

			if (status_n(n_LPT).eq.1) then

				call tracking()
! print *, 'LPT_P k1', k1(1:3)
IF (t<6.0D0) THEN ! Initializacion
k1=k1*0.0D0
END IF
				xp(n_LPT)=xpn(n_LPT)+dtp*k1(1)/two 
				yp(n_LPT)=ypn(n_LPT)+dtp*k1(2)/two 
				zp(n_LPT)=zpn(n_LPT)+dtp*k1(3)/two  !max(zero, zpn(3)+dtp*k1(3)/two)
				up(n_LPT)=upn(n_LPT)+dtp*k1(4)/two
				vp(n_LPT)=vpn(n_LPT)+dtp*k1(5)/two 
				wp(n_LPT)=wpn(n_LPT)+dtp*k1(6)/two

			end if
	end do !NUM
! print *, 'LPT_P xp xpn', xp, xpn
end subroutine lpt_p
