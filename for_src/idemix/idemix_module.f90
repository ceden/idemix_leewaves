
module idemix_module   
!=======================================================================
! module containing all relevant arrays and parameter for IDEMIX
!=======================================================================
      implicit none


      logical :: enable_idemix_lee_simple = .false.
     
      real*8, allocatable :: c0(:,:,:)                        ! mean vertical group velocity
      real*8, allocatable :: v0(:,:,:)                        ! mean lateral group velocity
      real*8, allocatable :: cstar(:,:)                       ! modal gravity wave speed
      real*8, allocatable :: alpha_c(:,:,:)                   ! dissipation parameter
  
      real*8, allocatable :: iw_diss(:,:,:)                   ! total dissipation of IW energy
       
     
      real*8 :: alpha_v=0                                      ! inverse time scale for vertical symmetrisation
      real*8 :: gamma=1.57                                    ! 
      real*8 :: jstar = 10.0                                  ! spectral bandwidth in modes
      real*8 :: mu0   = 4.0/3.0                               ! dissipation parameter
      
      real*8 :: E_gm = 1e-5
     
!---------------------------------------------------------------------------------
!     lee wave forcing
!---------------------------------------------------------------------------------
  
     
      real*8              :: nu = 0.8               ! Hurst number
      real*8              :: inv_fr_c = 0.75        ! critical inverse Froude number
      
      
      real*8, allocatable, dimension(:,:)  :: k_s   !Topographic wavenumber
      real*8, allocatable, dimension(:,:)  :: h_rms !RMS-height of topographic spectrum
      
      real*8, allocatable :: flux_lee_tot(:,:)
      real*8, allocatable :: lee_stress_u(:,:)
      real*8, allocatable :: inv_fr(:,:)            ! inverse Froude number 
      real*8, allocatable :: u_bot(:,:)             ! Bottom velocity in u-direction
      real*8, allocatable :: uz(:,:,:) 
      
      
      real*8, allocatable :: E_lee_d(:,:,:,:),dE_lee_d(:,:,:,:) 
      real*8, allocatable :: E_lee_s(:,:,:,:),dE_lee_s(:,:,:,:) 
      real*8, allocatable :: c_lee(:,:,:), Tc_lee(:,:,:) , lambda_0(:,:,:)
      real*8, allocatable :: mean_flow_to_E_lee(:,:,:) 
      real*8, allocatable :: iw_diss_lee(:,:,:) 
      
      real*8, allocatable :: E_lee_p(:,:,:,:),dE_lee_p(:,:,:,:)
      real*8, allocatable :: E_lee_m(:,:,:,:),dE_lee_m(:,:,:,:) 
      
end module idemix_module   



subroutine allocate_idemix_module
  use main_module
  use idemix_module
  implicit none
  
   allocate( cstar(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); cstar = 0
   allocate( c0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); c0 = 0
   allocate( v0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); v0 = 0
   allocate( alpha_c(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); alpha_c = 0
   allocate( iw_diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); iw_diss = 0
   allocate( uz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); uz = 0
      
     allocate( k_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) );       k_s   = 0.0
     allocate( h_rms(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) );     h_rms = 0.0
     allocate( flux_lee_tot(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); flux_lee_tot = 0
     allocate( inv_fr(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); inv_fr = 0
     allocate( u_bot(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) ); u_bot = 0
     allocate( lee_stress_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) );lee_stress_u=0
     
     allocate(  E_lee_d(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) );  E_lee_d = 0
     allocate( dE_lee_d(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dE_lee_d = 0
     allocate(  E_lee_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) );  E_lee_s = 0
     allocate( dE_lee_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dE_lee_s = 0
     
     allocate(  E_lee_p(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) );  E_lee_p = 0
     allocate(  E_lee_m(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) );  E_lee_m = 0
     allocate( dE_lee_p(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dE_lee_p = 0
     allocate( dE_lee_m(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3) ); dE_lee_m = 0
     
     allocate( c_lee(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); c_lee = 0
     allocate( lambda_0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); lambda_0 = 0
    
     allocate( Tc_lee(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); Tc_lee = 0
     allocate( mean_flow_to_E_lee(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); mean_flow_to_E_lee = 0
     allocate( iw_diss_lee(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); iw_diss_lee = 0
 
   
   
   
end subroutine allocate_idemix_module
