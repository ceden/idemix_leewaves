
!=======================================================================
! 1 D test case 
!======================================================================= 

module config_module
 ! use this module only locally in this file
 implicit none
 real*8, parameter :: N_0 = 30e-4
 real*8, parameter :: h0=2000.
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
    
 !enable_idemix = .true.

 alpha_v = 0.
 mu0 = 0.
 
 dt_mom     = 360
 dt_tracer  = dt_mom
 
 runlen =  86400.*4*30*10
 enable_diag_ts_monitor = .true.; ts_monint =  6*3600! dt_mom
 enable_diag_snapshots  = .true.; snapint  =   6*3600! dt_mom
 

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
        
       Nsqr(i,j,k) =  N_0**2
       Nsqrt(i,j,k) = N_0**2
     enddo
   enddo
 enddo


 do k=1,nz
   u(:,:,k,tau)    =  0.1
   u(:,:,k,taum1)  =  0.1
 enddo
 
 h_rms = 10
 k_s = 1./10e3

 
end subroutine set_initial_conditions


subroutine set_forcing
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





