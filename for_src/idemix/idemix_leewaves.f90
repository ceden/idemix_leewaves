


subroutine leewaves_simple
 use main_module
 use idemix_module
 implicit none
 integer :: i,j,k,ks
 real*8 :: fxa,L_u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: h(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: L_ww(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: L_s(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: L_tot(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 
  h = 0.
  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     h(i,j) = h(i,j) + dzw(k)*maskW(i,j,k)
    enddo
   enddo
  enddo
  
  L_u = 0.
  L_ww = 0.
  L_s = 0.
  do k=1,nz-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     !fxa = u_bot(i,j)*uz(i,j,k)*maskW(i,j,k)
     !L_u(i,j)  = L_u(i,j)  - dzw(k)*sign(1./max(1e-12,abs(fxa)),fxa )*u_bot(i,j)**2
     !L_ww(i,j) = L_ww(i,j) + dzw(k)*c_lee(i,j,k)/(2*alpha_c(i,j,k)*5e-5) 
     !L_s(i,j)  = L_s(i,j)  + dzw(k)*c_lee(i,j,k)/alpha_v 
     
     L_u(i,j)  = L_u(i,j)  - dzw(k)*u_bot(i,j)*uz(i,j,k)*maskW(i,j,k)/max(1e-12,u_bot(i,j)**2) 
     L_ww(i,j) = L_ww(i,j) + dzw(k)*(2*alpha_c(i,j,k)*E_gm) /max(1e-12,c_lee(i,j,k) )*maskW(i,j,k)
     L_s(i,j)  = L_s(i,j)  + dzw(k)*alpha_v/max(1e-12,c_lee(i,j,k) )*maskW(i,j,k)
     
    enddo
   enddo
  enddo
  !where( h>0 ) L_u =  sign( max(1e-12,abs(L_u)/h), L_u) 
  !where( h>0 ) L_ww = max(1e-12,L_ww/h ) 
  !where( h>0 ) L_s =  max(1e-12, L_s/h )
  where( h>0 ) L_u =  L_u/h
  where( h>0 ) L_ww = L_ww/h  
  where( h>0 ) L_s =  L_s/h 
  L_tot = sqrt(L_s*L_ww)
   
  do k=1,nz-1
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
      ks = max(1,kbot(i,j))
      fxa = flux_lee_tot(i,j)/c_lee(i,j,ks) 
      E_lee_d(i,j,k,tau) = maskW(i,j,k)*fxa*exp(-(zw(k)+h(i,j))*L_u(i,j) ) &
            *sinh(zw(k)*L_tot(i,j))/sinh(-h(i,j)*L_tot(i,j))
      !E_lee_d(i,j,k,tau) = maskW(i,j,k)*fxa*exp(-(zw(k)+h(i,j))/L_u(i,j) ) &
      !      *sinh(zw(k)/L_tot(i,j))/sinh(-h(i,j)/L_tot(i,j))   
    enddo
   enddo  
  enddo
  
  !print*, E_lee_d(is_pe,js_pe,nz,tau)
  !stop
  
  E_lee_p(:,:,:,:) = E_lee_d(:,:,:,:)
  E_lee_m(:,:,:,:) = 0   
  
  E_lee_s(:,:,:,:) = E_lee_p(:,:,:,:)   
  
  mean_flow_to_E_lee = Tc_lee*E_lee_d(:,:,:,tau)                                                                                                                                    
end subroutine 




subroutine superbee(flux,vel,var)
  use main_module
 use idemix_module
 implicit none
 integer :: k,i,j,kp2,km1
 real*8 :: flux(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,0:nz)
 real*8 :: var(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: vel(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: Rjp,Rj,Rjm,uCFL,Cr
 real*8 :: Limiter
 
 Limiter(Cr)=max(0.D0,max(min(1.D0,2.D0*Cr), min(2.D0,Cr))) 
 
 do k=1,nz-2
       kp2=min(nz-1,k+2); !if (kp2>np) kp2=3
       km1=max(1,k-1) !if (km1<1) km1=np-2
       do j=js_pe,je_pe
        do i=is_pe,ie_pe
         Rjp=(var(i,j,kp2)-var(i,j,k+1))!*maskW(i,j,k+1)
         Rj =(var(i,j,k+1)-var(i,j,k  ))!*maskW(i,j,k  )
         Rjm=(var(i,j,k  )-var(i,j,km1))!*maskW(i,j,km1)
         uCFL = ABS( vel(i,j,k)*dt_tracer/dzt(k) )
         IF (Rj.NE.0.) THEN
          IF (vel(i,j,k).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (vel(i,j,k).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         flux(i,j,k) = vel(i,j,k)*(var(i,j,k+1)+var(i,j,k))*0.5d0   &
                                -ABS(vel(i,j,k))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0
        enddo
       enddo
 enddo
end subroutine 





subroutine integrate_leewaves
 use main_module
 use idemix_module
 implicit none
 integer :: i,j,k,ks
 real*8 :: flux_p(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,0:nz)
 real*8 :: flux_m(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,0:nz)

  flux_m=0;flux_p=0
  call superbee(flux_p, c_lee,E_lee_p(:,:,:,tau))
  call superbee(flux_m,-c_lee,E_lee_m(:,:,:,tau))
  do i = is_pe,ie_pe
     do j = js_pe,je_pe
      ks = max(1,kbot(i,j))
      flux_m(i,j,ks-1) = - c_lee(i,j,ks)*E_lee_m(i,j,ks,tau)
      flux_p(i,j,ks-1) = flux_lee_tot(i,j) - flux_m(i,j,ks-1)
     enddo
  enddo
  k=nz-1; flux_p(:,:,k) =  c_lee(:,:,k)*E_lee_p(:,:,k,tau)
  k=nz-1; flux_m(:,:,k) = -flux_p(:,:,k)

  do k=1,nz-1
   dE_lee_p(:,:,k,tau) = -( flux_p(:,:,k)-flux_p(:,:,k-1) )/dzw(k) + Tc_lee(:,:,k)*E_lee_p(:,:,k,tau)
   dE_lee_m(:,:,k,tau) = -( flux_m(:,:,k)-flux_m(:,:,k-1) )/dzw(k) - Tc_lee(:,:,k)*E_lee_m(:,:,k,tau)
  enddo 
  E_lee_d(:,:,:,tau) = E_lee_p(:,:,:,tau) - E_lee_m(:,:,:,tau)
  E_lee_s(:,:,:,tau) = E_lee_p(:,:,:,tau) + E_lee_m(:,:,:,tau)
  
  mean_flow_to_E_lee = Tc_lee*E_lee_d(:,:,:,tau)
 
  iw_diss_lee = alpha_c/2*E_lee_s(:,:,:,tau)**2
  dE_lee_p(:,:,:,tau) =  dE_lee_p(:,:,:,tau) - alpha_v*E_lee_d(:,:,:,tau) - iw_diss_lee
  dE_lee_m(:,:,:,tau) =  dE_lee_m(:,:,:,tau) + alpha_v*E_lee_d(:,:,:,tau) - iw_diss_lee
                                        
  !iw_diss_lee = alpha_c*E_gm*(E_lee_p(:,:,:,tau) + E_lee_m(:,:,:,tau))
  !dE_lee_p(:,:,:,tau) =  dE_lee_p(:,:,:,tau) - alpha_v*E_lee_d(:,:,:,tau) - alpha_c*E_gm*E_lee_p(:,:,:,tau)
  !dE_lee_m(:,:,:,tau) =  dE_lee_m(:,:,:,tau) + alpha_v*E_lee_d(:,:,:,tau) - alpha_c*E_gm*E_lee_m(:,:,:,tau)
                                        
                                        
  E_lee_p(:,:,:,taup1) = E_lee_p(:,:,:,tau) +  maskW*dt_tracer*dE_lee_p(:,:,:,tau)  
  E_lee_m(:,:,:,taup1) = E_lee_m(:,:,:,tau) +  maskW*dt_tracer*dE_lee_m(:,:,:,tau)  
                                                                                                                                               
end subroutine 






subroutine leewave_flux
  use main_module
  use idemix_module
  implicit none
  integer   :: i,j,ks
  real*8 :: N_loc, f_loc   ! stability freq and abs of coriolis parameter
  real*8 :: fxa,mu
 
  nu = 0.8
  mu = 2*(nu+1)
  
  do j = js_pe,je_pe
   do i = is_pe,ie_pe
    ks = max(1,kbot(i,j))
    u_bot(i,j) = 0.5*( u(i,j,ks,tau)*maskU(i,j,ks) + u(i-1,j,ks,tau)*maskU(i-1,j,ks) )
   enddo
  enddo
 
  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,u_bot) 
  call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,u_bot)
  


 ! faster version with approximated isotropic spectrum
 
  do j = js_pe,je_pe
   do i  = is_pe,ie_pe
    ks = kbot(i,j)
    if (ks>0) then
     if ((Nsqr(i,j,ks) .gt. 0) .and. (Nsqr(i,j,ks) .gt. coriolis_t(i,j)**2)) then
      N_loc = sqrt(Nsqr(i,j,ks))
      f_loc = max(0.2e-4, abs(coriolis_t(i,j)) )
      fxa = (max(1.,N_loc/f_loc)-1.)**(nu-0.15)
      fxa = 8*1.3*pi*h_rms(i,j)**2*nu*k_s(i,j)**(mu-2)*N_loc**(-mu+4)*fxa
      flux_lee_tot(i,j) = fxa*abs(u_bot(i,j))**(mu-1)
      lee_stress_u(i,j) = -fxa*u_bot(i,j)*abs(u_bot(i,j))**(mu-3)
      
      ! account for critical Froud number
      inv_fr(i,j) = h_rms(i,j)*N_loc/max(abs(u_bot(i,j)),1d-16)
      if (inv_fr(i,j) > inv_fr_c ) then
           flux_lee_tot(i,j) = flux_lee_tot(i,j)*(inv_fr_c/(inv_fr(i,j)+1d-16))**2
           lee_stress_u(i,j) = lee_stress_u(i,j)*(inv_fr_c/(inv_fr(i,j)+1d-16))**2
           
      endif    
     endif  ! N>0 and N>f
    endif ! ks>0
   enddo
  enddo

  
end subroutine 











