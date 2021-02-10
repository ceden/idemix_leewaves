




program main
!=======================================================================
!      Top level driver 
!=======================================================================
      use main_module   
      use idemix_module   
      use diagnostics_module   
      use timing_module   
      implicit none
      integer :: otaum1

      if (my_pe==0) print '(/,a,f5.3)' ,'here is idemix3 version ',version

      n_pes_i = 1; n_pes_j = 1
      call tic('setup')
!---------------------------------------------------------------------------------
!      Initialize model
!---------------------------------------------------------------------------------
      itt=0
      call setup
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
      enditt = itt+int(runlen/dt_tracer)
      if (my_pe == 0 ) then
        print'(/,a,e8.2,a)','Starting integration for ',runlen,' s'
        print'(a,i10,a,i10,/)',' from time step ',itt,' to ',enditt
      endif
      call toc('setup')
!---------------------------------------------------------------------------------
!      Begin main model loop
!---------------------------------------------------------------------------------
      do while (itt < enditt) 
        call tic('main loop')
        call set_forcing

        call set_idemix_parameter

       

        call tic('idemix')
        if (enable_idemix_lee_simple) then
         call leewaves_simple
        else
         call integrate_leewaves
        endif
        call toc('idemix')

        if (enable_momentum_equation) call momentum

        call tic('diag')
        call diagnose
        call toc('diag')
 
        ! shift time 
        otaum1=taum1; taum1= tau; tau  = taup1; taup1= otaum1; itt=itt+1
!---------------------------------------------------------------------------------
!       End main model loop
!---------------------------------------------------------------------------------
      enddo
      if (my_pe==0) print'(/a)','end of integration'
end program main


subroutine setup
!=======================================================================
!  setup everything 
!=======================================================================
  use main_module   
  use idemix_module   
  implicit none
  integer :: k

  if (my_pe==0) print'(/a/)','setting up everything'

!--------------------------------------------------------------
! allocate everything
!--------------------------------------------------------------
  call set_parameter
  call pe_decomposition
  call allocate_main_module
  call allocate_idemix_module

!--------------------------------------------------------------
!  Grid
!--------------------------------------------------------------
  call set_grid
  call calc_grid


!--------------------------------------------------------------
! Coriolis
!--------------------------------------------------------------
  call set_coriolis

!--------------------------------------------------------------
! topography
!--------------------------------------------------------------
  call set_topography
  call calc_topo

!--------------------------------------------------------------
! initial condition and forcing 
!--------------------------------------------------------------
  call set_initial_conditions
  call set_forcing

  !call wgrid_to_tgrid( is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,Nsqr,Nsqrt)  
  do k=2,nz
   dNsqrdz(:,:,k) =  (Nsqr(:,:,k)-Nsqr(:,:,k-1))/dzt(k)
  enddo
  dNsqrdz(:,:,1) = dNsqrdz(:,:,2)
!--------------------------------------------------------------
! initialize diagnostics
!--------------------------------------------------------------
  call init_diagnostics
  call set_diagnostics

end subroutine setup






subroutine momentum
 use main_module   
 use idemix_module   
 implicit none
 integer :: k
 
    do k=1,nz-1
     flux_top(:,:,k) = -(2/pi)*Lambda_0(:,:,k)*E_lee_d(:,:,k,tau)*sign(1d0,u_bot)
    enddo          
   flux_top(:,:,nz) = 0.0 
   
   
   k=1; du(:,:,k,tau) = -(flux_top(:,:,k) - lee_stress_u(:,:))/dzt(k)*maskU(:,:,k)
   do k=2,nz
    du(:,:,k,tau) = -(flux_top(:,:,k) - flux_top(:,:,k-1))/dzt(k)*maskU(:,:,k)
   enddo
   
   do k=1,nz-1
    udiss(:,:,k) = (u(:,:,k+1,tau)-u(:,:,k,tau))*flux_top(:,:,k)/dzw(k)  
   enddo
   udiss(:,:,nz)=0.0
   
   ! explicit vertical friction
   if (.false.) then 
    do k=1,nz-1
     flux_top(:,:,k) = 2d-3*(u(:,:,k+1,tau)-u(:,:,k,tau))/dzw(k)  
    enddo   
    k=1; du(:,:,k,tau) = du(:,:,k,tau) + (flux_top(:,:,k) )/dzt(k)*maskU(:,:,k)
    do k=2,nz
     du(:,:,k,tau) = du(:,:,k,tau) + (flux_top(:,:,k) - flux_top(:,:,k-1))/dzt(k)*maskU(:,:,k)
    enddo
    do k=1,nz-1
     udiss(:,:,k) = udiss(:,:,k) - (u(:,:,k+1,tau)-u(:,:,k,tau))*flux_top(:,:,k)/dzw(k)  
    enddo
   endif
   u(:,:,:,taup1)= u(:,:,:,tau) + dt_mom*du(:,:,:,tau)*maskU
                                            
   call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,taup1)) 
   call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,u(:,:,:,taup1))
                                        


end subroutine momentum



