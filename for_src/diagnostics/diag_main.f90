



subroutine init_diagnostics
!=======================================================================
! initialize diagnostic routines
!=======================================================================
 use main_module
 use idemix_module
 use diagnostics_module   
 implicit none

 if (my_pe==0) print'(/,a)','Diagnostic setup:'

 if (enable_diag_ts_monitor) then
   if (my_pe==0) print'(a,e12.6,a,f10.2,a)',' time step monitor every ',ts_monint,' seconds/',ts_monint/dt_tracer,' time steps'
 endif

 if (enable_diag_snapshots) then
    if (my_pe==0) print'(a,e12.6,a,f10.2,a)',' writing snapshots every ',snapint,' seconds/',snapint/dt_tracer,' time steps'
    call init_snap_cdf_idemix
 endif

end subroutine init_diagnostics




subroutine diagnose
!=======================================================================
! call diagnostic routines
!=======================================================================
 use main_module   
 use diagnostics_module   
 use idemix_module   
 implicit none
 logical :: GM_strfct_diagnosed
 real*8 :: time
 real*8 :: fxa,fxb,fxc,fxd,fxu,fxl,fxcl,fxe
 integer :: i,j,k

 GM_strfct_diagnosed = .false.
 time = itt*dt_tracer

 if ( enable_diag_ts_monitor .and.  modulo(time,ts_monint) < dt_tracer ) then
    if (my_pe==0 )  print'(a,i10.10,a,e10.4,a)',' itt=',itt,' time=',time,'s'
       
   
       fxa = 0.
       fxb = 0.
       fxd = 0.
       fxu = 0.
       fxe = 0.
       do k=1,nz
        do i=is_pe,ie_pe
         do j=js_pe,je_pe
          fxa = fxa + mean_flow_to_E_lee(i,j,k)*dzw(k)
          fxb = fxb + udiss(i,j,k)*dzw(k)
          fxe = fxe + 2*iw_diss_lee(i,j,k)*dzw(k)
          fxd = fxd + dE_lee_p(i,j,k,tau)*dzw(k)
          fxd = fxd + dE_lee_m(i,j,k,tau)*dzw(k)
          fxu = fxu + u(i,j,k,tau)*du(i,j,k,tau)*dzt(k)
         enddo
        enddo
       enddo 
       call global_sum(fxa) 
       call global_sum(fxb) 
       call global_sum(fxd)
       call global_sum(fxu)  
       call global_sum(fxe)
       fxc = 0.
       fxcl = 0
       fxl = 0
        do i=is_pe,ie_pe
         do j=js_pe,je_pe
          fxc = fxc + flux_lee_tot(i,j)
          fxl = fxl - lee_stress_u(i,j)*u(i,j,1,tau)
         enddo
        enddo
        call global_sum(fxc) 
        call global_sum(fxcl) 
        call global_sum(fxl) 
        
        if (my_pe==0 )  then
          
          print*,' '
          print*,'diss of mean flow :',rho_0* fxb,' W/m^2 ' 
          print*,'lee wave flux     :',rho_0*fxl
          print*,'dU^2/2dt          :',rho_0*fxu
          print*,'error             :',rho_0*(fxu-fxb+fxl)
          print*,' '
          print*,'mean flow to lee  :',rho_0*fxa
          print*,'lee wave flux     :',rho_0*fxc
          print*,'dissipation       :',rho_0*fxe
          print*,'rate of change    :',rho_0*fxd
          print*,'error dE_lee_s/dt :',rho_0*(fxd-fxc-fxa+fxe)
          print*,' '
          print*,'error in exchange :',rho_0*(fxb+fxa),rho_0*(fxl-fxc)
          print*,' '
          
         endif
      
 
    call diag_snap_idemix
 endif

end subroutine diagnose




subroutine panic_snap
end subroutine panic_snap

