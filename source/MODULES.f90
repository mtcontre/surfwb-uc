MODULE global_variables
  !Global Variables MGP

  !Variables de Estado
  real (kind=8), dimension(:,:),save,allocatable	:: &
    V_global, C_global, VC, S1_global, S2_global

  !Metrics
  real (kind=8), dimension(:,:,:),save,allocatable:: &
    xi_global, eta_global, qnew_global, qold_global,qreal_global

  !q real es la variable dimensionalizada					  
  real (kind=8),dimension(:,:),save,allocatable::xc,yc,xe,ye, aj_global,MCoef


  !Input control
  real (kind=8)	::CFL,tinit,tfinal,L,H,U, t, iteration, dt,dtreal, FR2, g,dxi,&
    deta, hmin, kappa, Pvol, vol0, treal, Coef,batiopt,initqopt
  integer :: nitout=0
  real (kind=8) :: dtout
  logical :: print_out = .False.
  integer,dimension(4) :: CB
  integer	::caso,Nbx, Nby, it,dit, mmopt, rk, outopt, fopt, fM, Cf
  character (len=120),dimension(3) :: batiname,initqname
  character( len=255) :: outdir,indir

  !Lagrangian Particle tracking
  integer, parameter :: Nbz=1
END MODULE global_variables


!(jdgalaz@,2014): once parallelized and the overlapped grid system implemented
!the idea is to define new 'global variables' for each processor, just as
!when the code was not parallelized
!routine init.f90 takes care of this
MODULE mpi_surf
  logical,dimension(2)::isperiodic  
  logical::reorder
  integer,parameter::ndim=2,master=0
  integer,dimension(2)::dims,coords,coords2
  integer::comm2d,myleft,myright,myfront,myback,shift=1
  integer::ierror,myrank,nproc,myrank2d
  real (kind=8) :: time_start,time_finish, time_estim
END MODULE
module custombc
! Customized ghost cells for boundary conditions
integer:: nt_xi0g1, nt_xi0g2,nt_xiNg1, nt_xiNg2, nt_eta0g1,nt_eta0g2,nt_etaNg1,nt_etaNg2
integer:: optxi0g1,optxi0g2,optxiNg1,optxiNg2,opteta0g1,opteta0g2,optetaNg1,optetaNg2    !if 1 then use piecewise constant, if 2 use linear interpolation
real (kind=8) :: dt_xi0g1, dt_xi0g2,dt_xiNg1, dt_xiNg2, dt_eta0g1,dt_eta0g2,dt_etaNg1,dt_etaNg2
real (kind=8), dimension(:,:,:), allocatable :: qxi0g1,qxi0g2,qxiNg1,qxiNg2,qeta0g1,qeta0g2,qetaNg1,qetaNg2
integer:: flagxi0,flagxiN,flageta0,flagetaN !revisar init.f90 la parte del cfl inicial
real (kind=8), dimension(:,:), allocatable ::Sxi0,SxiN,Seta0,SetaN
! real (kind=8), dimension(:,:), allocatable S2xi0,S2xiN,S2eta0,S2etaN
end module
MODULE multigrid_surf  
  integer::ngrids !total number of grids
  integer::si,ei,sj,ej!first and last index in global matrices
  integer::myzone=1!current zone,should vary between groups of cores
!   integer, dimension(:),allocatable::myzone!current zone
  integer,dimension(4)::CB_real
  integer, dimension(:),allocatable::nxi,neta !vectors with size for each zone 
  integer,dimension(:),allocatable::batiopts,initqopts
  character(len=255),dimension(:,:),allocatable::batinames,initqnames
  
  type gridarray
      real(kind=8),dimension(:,:),allocatable :: X,Y,Z
  endtype gridarray
  
  type initqarray
    real(kind=8),dimension(:,:),allocatable ::H,U,V
  endtype initqarray
  
  type bMcoef!for buffering the Mcoef
    real(kind=8),dimension(:,:),allocatable ::bMcoef
  end type bMcoef
  
  !send/recv buffers
  type(gridarray), dimension(:),allocatable ::geom
  type(initqarray),dimension(:),allocatable ::initq
  type(bmcoef),dimension(:),allocatable ::buffMcoef 
  
  !to group boundaries
  integer,dimension(:),allocatable::inclxi0,inclxiN,incleta0,incletaN
  integer ::allgroup,groupxi0,groupxiN,groupeta0,groupetaN
  integer ::commxi0,commxiN,commeta0,commetaN
  
  !buffers for scattering input boundaries  

  real(kind=8),dimension(:,:,:),allocatable::bufqxi0g1,bufqxi0g2,&
    bufqxiNg1,bufqxiNg2, bufqeta0g1,bufqeta0g2, bufqetaNg1,bufqetaNg2
END MODULE

MODULE couplingbc
  !Coupling boundary conditions
  !specify explicit values for h,u,v
  !interpolate in time if neccesary
  integer:: nt_xi0g1, nt_xi0g2,nt_xiNg1, nt_xiNg2, nt_eta0g1,nt_eta0g2,nt_etaNg1,nt_etaNg2
  integer:: optxi0g1,optxi0g2,optxiNg1,optxiNg2,opteta0g1,opteta0g2,optetaNg1,optetaNg2    !if 1 then use piecewise constant, if 2 use linear interpolation
  real (kind=8) :: dt_xi0g1, dt_xi0g2,dt_xiNg1, dt_xiNg2, dt_eta0g1,dt_eta0g2,dt_etaNg1,dt_etaNg2
  real (kind=8), dimension(:,:,:), allocatable :: qxi0g1,qxi0g2,qxiNg1,qxiNg2,qeta0g1,qeta0g2,qetaNg1,qetaNg2
  integer:: flagxi0=0,flagxiN=0,flageta0=0,flagetaN=0 !revisar init.f90 la parte del cfl inicial
  real (kind=8), dimension(:,:), allocatable ::Sxi0,SxiN,Seta0,SetaN
  character(len=255),dimension(2) ::fnamexi0,fnamexiN,fnameeta0,fnameetaN
  !real (kind=8), dimension(:,:), allocatable S2xi0,S2xiN,S2eta0,S2etaN
END MODULE

MODULE coords
!Coordenadas Curvilíneas
  real (kind=8),dimension(:), save, allocatable :: coordxi, coordeta
  real (kind=8),dimension(:),save,allocatable :: angulo1,angulo2,angulo3,angulo4
END MODULE coords


MODULE geometries
  real (kind=8),dimension(:,:),save,allocatable	::x_global,y_global,z_global
  real (kind=8),dimension(:,:,:),save,allocatable	:: z_LPT
END MODULE geometries

MODULE senales
  integer:: GA1,GA2,GA3,GA4,Nsenal1,Nsenal2,Nsenal3,Nsenal4, IO1,IO2,IO3,IO4
  real (kind=8):: h01,h02,h03,h04
  real (kind=8),dimension(:,:),save,allocatable:: etaL1,qs1,us1,hs1, &
						  etaR2,qs2,us2,hs2,&
						  etaL3,qs3,us3,hs3, &
						  etaR4,qs4,us4,hs4, &
						  qA1,qA2,qA3,qA4, &
						  qsx1,qsy1, qsx2, qsy2, &
						  qsx3,qsy3, qsx4, qsy4, &
						  etaL9, hs9, &
						  etas1,etas2,etas3,etas4
						  
  real (kind=8),dimension(:),save,allocatable:: zA1,zA2,zA3,zA4,timeS1,timeS2,timeS3,timeS4,timeS9
  
 
END MODULE senales


MODULE time0
!Modulo que guarda las condiciones iniciales
real (kind=8), dimension(:,:,:),save,allocatable:: q0_global
real (kind=8):: epxR0, epyR0
END MODULE time0

!-----------------------
MODULE Jacobianos
!Jacobianos en interfaces de la celda
real (kind=8), dimension(:,:),save,allocatable:: Jac_global_xi, Jac_global_eta, Jac_global

END MODULE Jacobianos
!---------------------------
! MODULE Pich
!  real (kind=8),dimension(:),save,allocatable:: AreaPich
!  real (kind=8),dimension(:),save,allocatable:: VelPich
! END MODULE Pich

MODULE LagrangianPT ! Lagrangian Particle Tracking
real (kind=8),dimension(:),save,allocatable	::xp, yp, zp, up, vp, wp
real (kind=8),dimension(:),save,allocatable	::xpn, ypn, zpn, upn, vpn, wpn
real (kind=8),dimension(:),save,allocatable	::u_vel, v_vel, w_vel ! interpolated variables
integer,dimension(:),save,allocatable	::ii, jj, kk, status_n
!~ real (kind=8),dimension(:,:),save,allocatable	:: z_LPT
real (kind=8),dimension(:,:,:),save,allocatable	:: q_LPT
integer	:: num, ibck, ifwd, jbck, jfwd, kbck, kfwd, n_LPT
real (kind=8),parameter :: zero=0.0D0, one=1.0D0, two =2.0D0
  ! Integration variables
real (kind = 8),dimension(6)::k1 ,k2!,k3,k4,xold
integer :: treeflag! variables for kdtree search:
integer :: hitflag, iflag , LPT_init
END MODULE

MODULE TimeSeries ! Time Series at different points
!check init_TS.f90
integer, dimension(:),save,allocatable ::id0
real (kind=8), dimension(:),save,allocatable ::x0, y0, i0, j0, H0, Uts, Vts
integer, dimension(:,:),save, allocatable::m1,m2,m1_temp
logical, dimension(:), allocatable :: gaugeflag
real (kind=8), dimension(:),save,allocatable ::r, s, x 
  ! Interpolation by blending (see blend.f90, function blend_102)
real (kind=8), dimension(:,:),save,allocatable ::x00, x01, x10, x11 
  ! real (kind=8), dimension(:,:),save,allocatable :: x1,x2,y1,y2

integer	:: Nts ! Number of points for the time series
real (kind=8):: dt_TS ! timestep to record the time serie
real (kind=8):: sav_TS ! timestep to record the time serie
END MODULE

MODULE performance_stats
  
END MODULE
