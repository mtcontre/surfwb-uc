subroutine lpt_c() 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Program to perform the
! Lagrangian Particle Tracking of sediment particles
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

USE global_variables
USE LagrangianPT


  implicit none


  ! local variables
real (kind = 8) :: dtp, fr
character (len = 256) :: filename

!~ integer :: i, j, k, m, n, nv, lt

  ! interpolated variables
!~ real (kind = 8), dimension(:), allocatable :: u_vel,v_vel,w_vel

!~ integer :: iflag, tp, hitflag, np
!~ 
!~   ! variables for kdtree search:
!~ integer :: treeflag

! BEGIN!


dtp = dt / 2.0E+00  ! particle time-step

!---------------------------------------------------------------------
treeflag = 0

! allocate particle variables
!~ allocate (u_vel(num),v_vel(num),w_vel(num))

!~ u_vel=zero
!~ v_vel=zero
!~ w_vel=zero
!~ ibck=0
!~ ifwd=0
!~ jbck=0
!~ jfwd=0
!~ kbck=0
!~ kfwd=0



!-------- midpoint corrector ---------

	do n_LPT = 1,num
                hitflag = 0				
               
	          if (status_n(n_LPT).eq.1) then

				call tracking()
IF (t<6.0D0) THEN ! Initializacion
k1=k1*0.0D0
END IF
				xp(n_LPT)=xpn(n_LPT)+dtp*k1(1); 
				yp(n_LPT)=ypn(n_LPT)+dtp*k1(2); 
				zp(n_LPT)=zpn(n_LPT)+dtp*k1(3)  !max(zero, zpn(n)+dtp*k1(3))
				up(n_LPT)=upn(n_LPT)+dtp*k1(4); 
				vp(n_LPT)=vpn(n_LPT)+dtp*k1(5); 
				wp(n_LPT)=wpn(n_LPT)+dtp*k1(6)

 
              end if

                 xpn(n_LPT) = xp(n_LPT)
                 ypn(n_LPT) = yp(n_LPT)
                 zpn(n_LPT) = zp(n_LPT)

                 upn(n_LPT) = up(n_LPT)
                 vpn(n_LPT) = vp(n_LPT)
                 wpn(n_LPT) = wp(n_LPT)
	end do !NUM



!~ 	do n = 1,num
!~             hitflag = 0				
!~ 
!~ 		particle_loop: do tp = 1,3999
!~ 
!~ 	          if (status_n(n).eq.1) then
!~ 
!~ 				call tracking()
!~ 
!~ 				call rhs(xp(n),yp(n),zp(n),up(n),vp(n),wp(n),&
!~ 						u_vel(n),v_vel(n),w_vel(n),&
!~ 						Dux(n),Duy(n),Duz(n),&
!~                                     vortx(n),vorty(n),vortz(n),&
!~ 						fr,ds(n),k1(1:6))
!~ 
!~ 				xp(n)=xpn(n)+dtp*k1(1)/two; yp(n)=ypn(n)+dtp*k1(2)/two; zp(n)=zpn(n)+dtp*k1(3)/two  !max(zero, zpn(3)+dtp*k1(3)/two)
!~ 				up(n)=upn(n)+dtp*k1(4)/two; vp(n)=vpn(n)+dtp*k1(5)/two; wp(n)=wpn(n)+dtp*k1(6)/two
!~ 
!~                         call checkcollision()
!~ 
!~ 				call tracking()
!~ 
!~ 				call rhs(xp(n),yp(n),zp(n),up(n),vp(n),wp(n),&
!~ 						u_vel(n),v_vel(n),w_vel(n),&
!~ 						Dux(n),Duy(n),Duz(n),&
!~                                     vortx(n),vorty(n),vortz(n),&
!~ 						fr,ds(n),k1(1:6))
!~ 
!~ 				xp(n)=xpn(n)+dtp*k1(1); yp(n)=ypn(n)+dtp*k1(2); zp(n)=zpn(n)+dtp*k1(3)  !max(zero, zpn(n)+dtp*k1(3))
!~ 				up(n)=upn(n)+dtp*k1(4); vp(n)=vpn(n)+dtp*k1(5); wp(n)=wpn(n)+dtp*k1(6)
!~ 
!~                         call checkcollision()
!~ 
!~                  xpn(n) = xp(n)
!~                  ypn(n) = yp(n)
!~                  zpn(n) = zp(n)
!~ 
!~                  upn(n) = up(n)
!~                  vpn(n) = vp(n)
!~                  wpn(n) = wp(n)
!~ 
!~                  end if
!~ 
!~ 		end do particle_loop !TP
!~ 
!~ 	end do !NUM



!~ if (myid.eq.996) then !19) then
!~ n = 11 !30
!~ write(filename, fmt = '(a,i2.2)')'icheck',myid
!~ 
!~         			open (unit = myid, file = trim(filename), form = 'formatted', position = 'append')
!~         			write (unit = myid, fmt = *) xp(n),yp(n),zp(n),up(n),vp(n),wp(n),status_n(n)
!~ 
!~        			close (unit = myid)
!~ end if

!call mpi_barrier (mpi_comm_world, ierr)

!~ deallocate(u_vel, v_vel, w_vel)

!~ contains
!~   include 'tracking.f90'

end subroutine lpt_c
