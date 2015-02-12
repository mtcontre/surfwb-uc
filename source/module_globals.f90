MODULE global_variables
  logical::goout=.false.
  !Global Variables MGP

  !Variables de Estado
  real (kind=8), dimension(:,:),save,allocatable	:: V_global, C_global, VC, S1_global, S2_global

  !Metrics
  real (kind=8), dimension(:,:,:),save,allocatable	::xi_global, eta_global, &
					    qnew_global, qold_global, &
					    qreal_global

  !q real es la variable dimensionalizada					  
  real (kind=8),dimension(:,:),save,allocatable	::xc,yc,xe,ye, aj_global,MCoef


  !Input control
  real (kind=8)	::CFL,tfinal,tinit,L,H,U, t, iteration, dt, FR2, g,dxi, &
  deta, hmin, kappa, Pvol, vol0, treal, Coef
  integer :: batiopt, initqopt
  integer :: it,dit
  integer :: nitout=0
  real (kind=8) :: dtout
  logical :: print_out = .False.
  integer,dimension(4) :: CB
  integer	::caso,Nbx, Nby, nxi, neta, mmopt, rk, outopt, fopt, fM, Cf
  character (len=120),dimension(3) :: batiname,initqname
  character(len=256)::outdir
  !Lagrangian Particle tracking
  integer, parameter :: Nbz=1
END MODULE global_variables