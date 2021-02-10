
!=======================================================================
! 1 D test case 
!======================================================================= 

module config_module
 ! use this module only locally in this file
 implicit none
 real*8, parameter :: N_0 = 30e-4,z_b = 800.0  !N_0 = 0.008
 real*8, parameter :: zu0=200.0,zu1=750
 real*8, parameter :: h0=2000.
 real*8, allocatable :: u_target(:,:,:)
end module config_module


subroutine set_parameter
 ! ----------------------------------
 !       set here main parameter
 ! ----------------------------------
 use main_module   
 use config_module
 use diagnostics_module   
 use idemix_module   
 implicit none
 
 nx   = 1; nz   = 200; ny   = 1

 coord_degree           = .false.
 enable_cyclic_x        = .true.
 enable_cyclic_y        = .true.
 eq_of_state_type       = 5
      
 !enable_idemix_lee_simple = .true.
 
 alpha_v = 0.
 mu0 = 0.
 
 dt_mom     = 360
 dt_tracer  = dt_mom
 E_gm = 1e-4
 runlen =  86400.*120*3
 enable_diag_ts_monitor = .true.; ts_monint =  6*3600! dt_mom
 enable_diag_snapshots  = .true.; snapint  =   6*3600! dt_mom
 
 enable_momentum_equation = .false.
end subroutine set_parameter



subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 dxt    = 5e3
 dyt    = 5e3
 dzt    = h0/nz
end subroutine set_grid

subroutine set_coriolis
 use main_module   
 use config_module
 implicit none
 coriolis_t = 1e-4
end subroutine set_coriolis


subroutine set_initial_conditions
 use main_module
 use idemix_module   
 use config_module
 implicit none
 integer :: i,j,k


 do k=1,nz
   do j=js_pe-onx,je_pe+onx
      do i=is_pe-onx,ie_pe+onx    
       Nsqr(i,j,k)  = ( 2e-4 + N_0*exp(zw(k)/z_b)  )**2
       Nsqrt(i,j,k) = ( 2e-4 + N_0*exp(zt(k)/z_b)  )**2      
     enddo
   enddo
 enddo

 allocate( u_target(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 do k=1,nz
   u_target(:,:,k) =  0.1+0.1*exp( -(zt(k)+zu1)**2/zu0**2)
   u(:,:,k,tau)    =  u_target(:,:,k)
   u(:,:,k,taum1)  =  u_target(:,:,k)
 enddo
 
 h_rms = 10
 k_s = 1./10e3

end subroutine set_initial_conditions


subroutine set_forcing
 use main_module   
 use config_module
  use idemix_module   
 implicit none
 !u(:,:,:,tau) = u(:,:,:,tau) + dt_tracer*(.02/86400.)*(u_target-u(:,:,:,taum1))
 
                                        
  E_lee_p(:,:,:,tau) = E_lee_p(:,:,:,tau) -  maskW*dt_tracer*0.01*1e-4*E_lee_p(:,:,:,tau)    
  E_lee_m(:,:,:,tau) = E_lee_m(:,:,:,tau) -  maskW*dt_tracer*0.01*1e-4*E_lee_m(:,:,:,tau)  
    
end subroutine set_forcing

subroutine set_topography
 use main_module   
 use config_module   
 kbot=1
end subroutine set_topography


subroutine set_diagnostics
end subroutine set_diagnostics

subroutine set_particles
end subroutine set_particles





