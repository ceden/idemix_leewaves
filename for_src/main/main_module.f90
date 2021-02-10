

module main_module   
!=======================================================================
! main module containing most important arrays and parameter
! others can be found in specific modules
!=======================================================================
      implicit none
!---------------------------------------------------------------------------------
!     constants and parameter
!---------------------------------------------------------------------------------
      real*8, parameter :: version = 2.11
      real*8, parameter :: pi      = 3.14159265358979323846264338327950588
      real*8, parameter :: radius  = 6370.0e3        ! Earth radius in m
      real*8, parameter :: degtom  = radius/180.0*pi ! conversion degrees latitude to meters
      real*8, parameter :: mtodeg  = 1/degtom        ! revers conversion 
      real*8, parameter :: omega   = pi/43082.0      ! earth rotation frequency in 1/s
      real*8, parameter :: rho_0   = 1024.0          ! Boussinesq reference density in kg/m^3
      real*8, parameter :: grav    = 9.81            ! gravitational constant in m/s^2
!---------------------------------------------------------------------------------
!     Parallel domain setup
!---------------------------------------------------------------------------------
      integer :: n_pes     ! total number of processors
      integer :: my_pe     ! index of this processor from 0 to n_pes-1
      integer :: n_pes_i   ! total number of processors in x direction
      integer :: n_pes_j   ! total number of processors in y direction
      integer :: my_blk_i  ! index of this processor in x direction from 1 to n_pes_i
      integer :: my_blk_j  ! index of this processor in y direction from 1 to n_pes_j
      integer :: i_blk     ! grid points of domain decompostion in x direction 
      integer :: j_blk     ! grid points of domain decompostion in y direction
      integer :: is_pe     ! start index of grid points in x direction of this processor
      integer :: ie_pe     ! end index of grid points in x direction of this processor
      integer :: js_pe     ! start index of grid points in y direction of this processor
      integer :: je_pe     ! end index of grid points in y direction of this processor
      integer :: onx=2     ! number of overlapping points in x and y direction
      integer :: my_comm=0 ! communicator for MPI library
!---------------------------------------------------------------------------------
!     model parameter
!---------------------------------------------------------------------------------
      integer :: nx            ! grid points in zonal (x,i) direction
      integer :: ny            ! grid points in meridional (y,j) direction
      integer :: nz            ! grid points in vertical (z,k) direction
      integer :: taum1     = 1 ! pointer to last time step  
      integer :: tau       = 2 ! pointer to current time step
      integer :: taup1     = 3 ! pointer to next time step
      real*8  :: dt_mom    = 0 ! time step in seconds for momentum
      real*8  :: dt_tracer = 0 ! time step for tracer can be larger than for momentum
      integer :: itt           ! time step number
      integer :: enditt        ! last time step of simulation
      real*8  :: runlen=0.     ! length of simulation in seconds
      real*8  :: AB_eps = 0.1  ! deviation from Adam-Bashforth weighting
!---------------------------------------------------------------------------------
!     logical switches for general model setup
!---------------------------------------------------------------------------------
      logical :: coord_degree                      = .false. ! either spherical (true) or cartesian (false) coordinates 
      logical :: enable_cyclic_x                   = .false. ! enable zonal cyclic boundary conditions 
      logical :: enable_cyclic_y                   = .false. ! enable meridional cyclic boundary conditions
      integer :: eq_of_state_type = 1                        ! equation of state: 1: linear, 3: nonlinear with comp., 5: TEOS
      logical :: enable_momentum_equation          = .true.  ! calculate momentum changes                
!---------------------------------------------------------------------------------
!     variables related to numerical grid
!---------------------------------------------------------------------------------
      real*8, allocatable, dimension(:,:,:)   :: maskT     ! mask in physical space for tracer points
      real*8, allocatable, dimension(:,:,:)   :: maskU     ! mask in physical space for U points
      real*8, allocatable, dimension(:,:,:)   :: maskV     ! mask in physical space for V points
      real*8, allocatable, dimension(:,:,:)   :: maskW     ! mask in physical space for W points
      real*8, allocatable, dimension(:,:,:)   :: maskZ     ! mask in physical space for Zeta points
      integer, allocatable, dimension(:,:)    :: kbot       ! 0 denotes land, 0<kmt<=nz denotes deepest cell zt(kmt)
      real*8, allocatable, dimension(:)       :: xt,dxt     ! zonal (x) coordinate of T-grid point in meters
      real*8, allocatable, dimension(:)       :: xu,dxu     ! zonal (x) coordinate of U-grid point in meters
      real*8, allocatable, dimension(:)       :: yt,dyt     ! meridional (y) coordinate of T-grid point in meters
      real*8, allocatable, dimension(:)       :: yu,dyu     ! meridional (y) coordinate of V-grid point in meters
      real*8                                  :: x_origin,y_origin ! origin of grid in x and y direction, located at xu_1, yu_1
      real*8, allocatable, dimension(:)       :: zt,zw      ! vertical coordinate in m
      real*8, allocatable, dimension(:)       :: dzt,dzw    ! box thickness in m
      real*8, allocatable, dimension(:,:)     :: coriolis_t ! coriolis frequency at T grid point in 1/s
      real*8, allocatable, dimension(:)       :: cost       ! metric factor for spherical coordinates on T grid
      real*8, allocatable, dimension(:)       :: cosu       ! metric factor for spherical coordinates on U grid
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
      real*8, allocatable, dimension(:,:,:) :: Nsqr                  ! Square of stability frequency in 1/s^2
      real*8, allocatable, dimension(:,:,:) :: Nsqrt                 
      real*8, allocatable, dimension(:,:,:) :: dNsqrdz                 
      real*8, allocatable, dimension(:,:,:,:) :: u,du!,v  
      real*8, allocatable, dimension(:,:,:) :: udiss                 !
      real*8, allocatable, dimension(:,:,:)   :: flux_east,flux_north,flux_top ! multi purpose fluxes

end module main_module   


subroutine allocate_main_module
!=======================================================================
! allocate all arrays within main module
!=======================================================================
 use main_module   
 implicit none



  allocate( xt(is_pe-onx:ie_pe+onx), xu(is_pe-onx:ie_pe+onx)) ; xt=0;xu=0
  allocate( yt(js_pe-onx:je_pe+onx), yu(js_pe-onx:je_pe+onx)) ; yt=0;yu=0
  allocate( dxt(is_pe-onx:ie_pe+onx), dxu(is_pe-onx:ie_pe+onx)) ; dxt=0;dxu=0
  allocate( dyt(js_pe-onx:je_pe+onx), dyu(js_pe-onx:je_pe+onx))  ; dyt=0;dyu=0

  allocate( zt(nz), dzt(nz), zw(nz), dzw(nz) ); zt=0; zw=0; dzt=0; dzw=0
  allocate( cost(js_pe-onx:je_pe+onx), cosu(js_pe-onx:je_pe+onx)); cost=1.0; cosu=1.0;
  allocate( coriolis_t(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); coriolis_t=0

  allocate( kbot(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); kbot=0

  allocate( maskT(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
  allocate( maskU(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
  allocate( maskV(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
  allocate( maskW(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
  allocate( maskZ(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
  maskW=0.; maskT=0.; maskU=0.; maskV=0.; maskZ =0.;

  allocate(Nsqr(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); Nsqr = 0
  allocate(Nsqrt(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); Nsqrt = 0
  allocate(dNsqrdz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); dNsqrdz = 0
  allocate( u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); u = 0
  allocate( du(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); du = 0
  !allocate( v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); v = 0

  allocate( udiss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); udiss = 0
  allocate( flux_east(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); flux_east = 0
  allocate(flux_north(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); flux_north = 0
  allocate(  flux_top(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); flux_top = 0


end subroutine allocate_main_module





