


subroutine set_idemix_parameter
!=======================================================================
! set main IDEMIX parameter 
!=======================================================================
  use main_module   
  use idemix_module   
  implicit none
  !include "mass.include"  ! include this on AIX which does not know function acosh, also link with -lmass
  real*8 :: bN0,fxa
  real*8 :: gofx2,hofx1
  integer :: i,j,k
  real*8 :: floc,Nloc!,fxa,fxb
  
 !---------------------------------------------------------------------------------
 ! Idemix 1 parameter
 !---------------------------------------------------------------------------------
  do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
    bN0=0.0
    do k=1,nz-1
     bN0 = bN0 + max(0d0,Nsqr(i,j,k))**0.5*dzw(k)*maskW(i,j,k) 
    enddo
    bN0 = bN0 + max(0d0,Nsqr(i,j,nz))**0.5*0.5*dzw(nz)*maskW(i,j,nz) 
    cstar(i,j) = max(1d-2,bN0/(pi*jstar) )
    do k=1,nz
     fxa = max(0d0,Nsqr(i,j,k))**0.5/(1d-22 + abs(coriolis_t(i,j)) )
     alpha_c(i,j,k) = mu0*acosh(max(1d0,fxa))*abs(coriolis_t(i,j))/cstar(i,j)**2 *maskW(i,j,k) 
     c0(i,j,k)=max(0d0, gamma*cstar(i,j)*gofx2(fxa)*maskW(i,j,k) )
     v0(i,j,k)=max(0d0, gamma*cstar(i,j)*hofx1(fxa)*maskW(i,j,k) )
    enddo
   enddo
  enddo

 !---------------------------------------------------------------------------------
 !  Shear
 !---------------------------------------------------------------------------------
 do k=1,nz-1
    uz(:,:,k) = (u(:,:,k+1,tau)-u(:,:,k,tau))/dzw(k)*maskU(:,:,k)*maskU(:,:,k+1)
 enddo
 uz(:,:,nz)=uz(:,:,nz-1); 
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,uz) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,uz)

 do k=1,nz
   do j=js_pe-onx,je_pe+onx
     do i=is_pe-onx,ie_pe+onx
        floc = max(abs(coriolis_t(i,j)),1D-6)
        Nloc = max(floc*1.2, sqrt( max(0d0, Nsqr(i,j,k))  ) )
        Lambda_0(i,j,k) = 0.65*floc/Nloc/log(10d0)
        
        !fxa = 0.5*floc/Nloc*(Nloc**2+0.5*floc**2)/(Nloc**2-floc**2 ) 
        !fxa = fxa * log( (Nloc+sqrt(Nloc**2-floc**2) )/(Nloc-sqrt(Nloc**2-floc**2)) )
        !fxa = fxa - 3./2.*floc/sqrt(Nloc**2-floc**2 ) 
        !Lambda_0(i,j,k) = fxa*(2./pi)/(1-2./pi*asin(floc/Nloc) )*maskW(i,j,k)
    enddo
   enddo
  enddo


  call leewave_flux 
 

  do k=1,nz-1
   do j = js_pe,je_pe
    do i = is_pe,ie_pe
     c_lee(i,j,k)  = (2/pi)*0.5*(Lambda_0(i,j,k)+Lambda_0(i,j,k+1))*abs( u_bot(i,j) )
     Tc_lee(i,j,k) = (2/pi)*Lambda_0(i,j,k)*sign(1d0,u_bot(i,j) )*0.5*(uz(i,j,k)+uz(i-1,j,k))                                                     
    enddo
   enddo 
  enddo



end subroutine set_idemix_parameter










function gofx2(x)
!=======================================================================
! a function g(x) 
!=======================================================================
 implicit none
 real*8 :: gofx2,x,c
 real*8, parameter :: pi = 3.14159265358979323846264338327950588
 x=max(3d0,x)
 c= 1.-(2./pi)*asin(1./x)
 gofx2 = 2/pi/c*0.9*x**(-2./3.)*(1-exp(-x/4.3))
end function gofx2

function hofx1(x)
!=======================================================================
! a function h(x) 
!=======================================================================
 implicit none
 real*8 :: hofx1,x
 real*8, parameter :: pi = 3.14159265358979323846264338327950588
 hofx1 = (2./pi)/(1.-(2./pi)*asin(1./x)) * (x-1.)/(x+1.)
end function hofx1

