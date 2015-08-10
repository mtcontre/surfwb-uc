
subroutine tracking ()

USE global_variables
USE LagrangianPT
USE geometries

  implicit none

  ! local variables
integer ::nn


! START TRACKING

!~ print *, 'tracking 1'

ibck = ii(n_LPT)-3
ifwd = ii(n_LPT)+3
jbck = jj(n_LPT)-3
jfwd = jj(n_LPT)+3
kbck = kk(n_LPT)!-3
kfwd = kk(n_LPT)!+3

iflag = 0

!~ print *, 'TRACKING u_vel', u_vel

call locate(xp(n_LPT), yp(n_LPT), zp(n_LPT), ii(n_LPT), jj(n_LPT), kk(n_LPT), &
	u_vel(n_LPT), v_vel(n_LPT), w_vel(n_LPT), &
	status_n(n_LPT))
!~ print *, 'TRACKING 2 u_vel', u_vel
!~ print *, 'TRACKING xp Locate ', xp(n_LPT), yp(n_LPT), zp(n_LPT)
if (iflag.eq.0) then

!~ if(ii(n).ge.im(nz)-2) ii(n)= 4

!~ if(ii(n).le.2) ii(n)=im(nz)-5 
 
ibck = ii(n_LPT)-5
ifwd = ii(n_LPT)+5
jbck = jj(n_LPT)-5
jfwd = jj(n_LPT)+5
kbck = kk(n_LPT)!-5
kfwd = kk(n_LPT)!+5

call locate( xp(n_LPT), yp(n_LPT), zp(n_LPT), ii(n_LPT), jj(n_LPT), kk(n_LPT), &
        u_vel(n_LPT), v_vel(n_LPT), w_vel(n_LPT), &
        status_n(n_LPT))

end if



if (iflag.eq.0) then

treeflag = 1
ibck = 1; ifwd = Nbx - 1
jbck = 1; jfwd = Nby - 1
kbck = 1; kfwd = 1

	call locate(xp(n_LPT), yp(n_LPT), zp(n_LPT), ii(n_LPT), jj(n_LPT), kk(n_LPT), &
		u_vel(n_LPT), v_vel(n_LPT), w_vel(n_LPT), &
		status_n(n_LPT))

end if


! if particle is not found, it exited the domain
! and it won't be calculated any longer.
if (iflag.eq.0) then
	status_n(n_LPT) = 0
      print *,'particle=',n_LPT,'exited'
      print *,'x,y,z',xp(n_LPT),yp(n_LPT),zp(n_LPT),status_n(n_LPT)
      print *,'last known position',ii(n_LPT),jj(n_LPT),kk(n_LPT)
end if

treeflag = 0

end subroutine tracking

