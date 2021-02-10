




subroutine u_centered_grid(dyt,dyu,yt,yu,n)
!---------------------------------------------------------------------------------
! setup u-centered grid based in Delta yt and the relations
! dyt_i = yu_i - yu_i-1 , yu_i = 0.5(yt_i+yt_(i+1)) , dyu_i = yt_(i+1)-yt_i
!---------------------------------------------------------------------------------
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: dyt(n)
  real*8, intent(out) :: yu(n),yt(n),dyu(n)
  integer :: i 
  yu(1)=0
  do i=2,n  
   yu(i)=yu(i-1)+dyt(i)
  enddo
  yt(1)=yu(1)-dyt(1)*0.5
  do i=2,n
   yt(i) = 2*yu(i-1) - yt(i-1)
  enddo
  do i=1,n-1 
   dyu(i)= yt(i+1)-yt(i)
  enddo
  dyu(n)=2*dyt(n)- dyu(max(1,n-1))
end subroutine u_centered_grid





subroutine calc_grid
!---------------------------------------------------------------------------------
!    setup grid based on dxt,dyt,dzt and x_origin, y_origin
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: i,j
 real*8 :: aloc(nx,ny)
 real*8, dimension(1-onx:nx+onx) :: dxt_gl,dxu_gl,xt_gl,xu_gl
 real*8, dimension(1-onx:ny+onx) :: dyt_gl,dyu_gl,yt_gl,yu_gl

  aloc=0.
!--------------------------------------------------------------
! transfer from locally defined variables to global ones
!--------------------------------------------------------------
  aloc(is_pe:ie_pe,1) = dxt(is_pe:ie_pe)
  call pe0_recv_2D(nx,ny,aloc)
  call pe0_bcast(aloc,nx*ny)
  dxt_gl(1:nx) = aloc(:,1)

  if (enable_cyclic_x) then
   do i=1,onx
      dxt_gl(nx+i)=dxt_gl(i); dxt_gl(1-i)=dxt_gl(nx-i+1) 
   enddo
  else
   do i=1,onx
      dxt_gl(nx+i)=dxt_gl(nx); dxt_gl(1-i)=dxt_gl(1) 
   enddo
  endif

  aloc(1,js_pe:je_pe) = dyt(js_pe:je_pe)
  call pe0_recv_2D(nx,ny,aloc)
  call pe0_bcast(aloc,nx*ny)
  dyt_gl(1:ny) = aloc(1,:)

  if (enable_cyclic_y) then
   do j=1,onx
      dyt_gl(ny+j)=dyt_gl(j); dyt_gl(1-j)=dyt_gl(ny-j+1) 
   enddo
  else
   do j=1,onx
      dyt_gl(ny+j)=dyt_gl(ny); dyt_gl(1-j)=dyt_gl(1) 
   enddo
  endif
!--------------------------------------------------------------
! grid in east/west direction
!--------------------------------------------------------------
  call u_centered_grid(dxt_gl,dxu_gl,xt_gl,xu_gl,nx+2*onx)
  xt_gl=xt_gl-xu_gl(1)+x_origin
  xu_gl=xu_gl-xu_gl(1)+x_origin

  if (enable_cyclic_x) then
   do i=1,onx
       xt_gl(nx+i)=xt_gl(i); xt_gl(1-i)=xt_gl(nx-i+1) 
       xu_gl(nx+i)=xt_gl(i); xu_gl(1-i)=xu_gl(nx-i+1) 
       dxu_gl(nx+i)=dxu_gl(i); dxu_gl(1-i)=dxu_gl(nx-i+1) 
   enddo
  endif

!--------------------------------------------------------------
! grid in north/south direction
!--------------------------------------------------------------
  call u_centered_grid(dyt_gl,dyu_gl,yt_gl,yu_gl,ny+2*onx)
  yt_gl=yt_gl-yu_gl(1)+y_origin
  yu_gl=yu_gl-yu_gl(1)+y_origin

  if (enable_cyclic_y) then
   do j=1,onx
       yt_gl(ny+j)=yt_gl(j); yt_gl(1-j)=yt_gl(ny-j+1) 
       yu_gl(ny+j)=yt_gl(j); yu_gl(1-j)=yu_gl(ny-j+1) 
       dyu_gl(ny+j)=dyu_gl(j); dyu_gl(1-j)=dyu_gl(ny-j+1) 
   enddo
  endif
  
  if (coord_degree) then
!--------------------------------------------------------------
! convert from degrees to pseudo cartesian grid
!--------------------------------------------------------------
    dxt_gl=dxt_gl*degtom; dxu_gl=dxu_gl*degtom;
    dyt_gl=dyt_gl*degtom; dyu_gl=dyu_gl*degtom;
  endif

!--------------------------------------------------------------
!  transfer to locally defined variables
!--------------------------------------------------------------
  xt(is_pe-onx:ie_pe+onx)  = xt_gl(is_pe-onx:ie_pe+onx)
  xu(is_pe-onx:ie_pe+onx)  = xu_gl(is_pe-onx:ie_pe+onx)
  dxu(is_pe-onx:ie_pe+onx) = dxu_gl(is_pe-onx:ie_pe+onx)
  dxt(is_pe-onx:ie_pe+onx) = dxt_gl(is_pe-onx:ie_pe+onx)

  yt(js_pe-onx:je_pe+onx)  = yt_gl(js_pe-onx:je_pe+onx)
  yu(js_pe-onx:je_pe+onx)  = yu_gl(js_pe-onx:je_pe+onx)
  dyu(js_pe-onx:je_pe+onx) = dyu_gl(js_pe-onx:je_pe+onx)
  dyt(js_pe-onx:je_pe+onx) = dyt_gl(js_pe-onx:je_pe+onx)

!--------------------------------------------------------------
! grid in vertical direction
!--------------------------------------------------------------
  call u_centered_grid(dzt,dzw,zt,zw,nz)
  !dzw(nz)=dzt(nz) !*0.5 ! this is account for in the model directly
  zt = zt - zw(nz); zw = zw - zw(nz)  ! zero at zw(nz) 

!--------------------------------------------------------------
! metric factors
!--------------------------------------------------------------
  if (coord_degree) then
   do j=js_pe-onx,je_pe+onx
    cost(j) = cos( yt(j)/180.*pi ) 
    cosu(j) = cos( yu(j)/180.*pi ) 
    !tantr(j) = tan( yt(j)/180.*pi ) /radius
   enddo
  else
   cost=1.0;cosu=1.0;!tantr=0.0
  endif

end subroutine calc_grid



subroutine calc_topo
!--------------------------------------------------------------
! calulate masks, total depth etc
!--------------------------------------------------------------
 use main_module
 implicit none
 integer :: i,j,k

!--------------------------------------------------------------
! close domain
!--------------------------------------------------------------
  if (.not. enable_cyclic_y) then
    if (my_blk_j == 1)         kbot(:,1-onx:0)=0
    if (my_blk_j == n_pes_j)   kbot(:,ny+1:ny+onx)=0
  endif
  if (.not. enable_cyclic_x) then
    if (my_blk_i == 1)         kbot(1-onx:0,:)=0
    if (my_blk_i == n_pes_i)   kbot(nx+1:nx+onx,:)=0  
  endif
  call border_exchg_xy_int(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,kbot) 
  call setcyclic_xy_int   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,kbot)

!--------------------------------------------------------------
! Land masks
!--------------------------------------------------------------
  maskT = 0.0
  do k=1,nz
   do j=js_pe-onx,je_pe+onx
      do i=is_pe-onx,ie_pe+onx
        if ( kbot(i,j)/=0 .and. kbot(i,j) <= k ) maskT(i,j,k)=1.0
      enddo
   enddo
  enddo
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskT) 
  call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskT)
  maskU=maskT
  do i=is_pe-onx,ie_pe+onx-1
     maskU(i,:,:)=min(maskT(i,:,:),maskT(i+1,:,:))
  enddo
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskU) 
  call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskU)
  maskV=maskT
  do j=js_pe-onx,je_pe+onx-1
      maskV(:,j,:)=min(maskT(:,j,:),maskT(:,j+1,:))
  enddo
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskV) 
  call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskV)
  maskZ=maskT
  do j=js_pe-onx,je_pe+onx-1
   do i=is_pe-onx,ie_pe+onx-1
     maskZ(i,j,:)=min(maskT(i,j,:),maskT(i,j+1,:),maskT(i+1,j,:))
   enddo
  enddo
  call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskZ) 
  call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,maskZ)
  maskW=maskT
  do k=1,nz-1
    maskW(:,:,k)=min(maskT(:,:,k),maskT(:,:,k+1))
  enddo
end subroutine calc_topo










subroutine solve_tridiag(a,b,c,d,x,n)
      implicit none
!---------------------------------------------------------------------------------
!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!        b - the main diagonal
!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!        d - right part
!        x - the answer
!        n - number of equations
!---------------------------------------------------------------------------------
        integer,intent(in) :: n
        real*8,dimension(n),intent(in) :: a,b,c,d
        real*8,dimension(n),intent(out) :: x
        real*8,dimension(n) :: cp,dp
        real*8 :: m,fxa
        integer i
 
! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           fxa = 1D0/m
           cp(i) = c(i)*fxa
           dp(i) = (d(i)-dp(i-1)*a(i))*fxa
         enddo
! initialize x
         x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do
end subroutine solve_tridiag



subroutine wgrid_to_tgrid( is_,ie_,js_,je_,nz_,A,B)  
!---------------------------------------------------------------------------------
! global integral conserving interpolation from W to T grid
! for U-centered boxes, A and B can be identical
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: is_,ie_,js_,je_,nz_,k
 real*8, dimension(is_:ie_,js_:je_,nz_) :: A,B
 real*8 :: aloc(is_:ie_,js_:je_)
 
 aloc=A(:,:,nz)
 k=nz; B(:,:,k)=(0.5*dzw(k)*A(:,:,k)+dzw(k-1)*A(:,:,k-1))/(2*dzt(k))   !!!
 do k=nz-1,2,-1
  B(:,:,k)=(dzw(k)*A(:,:,k)+dzw(k-1)*A(:,:,k-1))/(2*dzt(k))  
 enddo
 B(:,:,1)=(dzw(1)*A(:,:,1))/(2*dzt(1))  
 !B(:,:,nz)=B(:,:,nz) + (dzw(nz)*aloc)/(2*dzt(nz))   !! dzw(nz) is only half as large
 B(:,:,nz)=B(:,:,nz) + (0.5*dzw(nz)*aloc)/(2*dzt(nz))   !! dzw(nz) is only half as large
end subroutine wgrid_to_tgrid



subroutine tgrid_to_wgrid( is_,ie_,js_,je_,nz_,A,B) 
!---------------------------------------------------------------------------------
! global integral conserving interpolation from T to W grid
! for U-centered boxes, A and B can be identical
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: is_,ie_,js_,je_,nz_,k
 real*8, dimension(is_:ie_,js_:je_,nz_) :: A,B
 real*8 :: aloc(is_:ie_,js_:je_)
 
 aloc=A(:,:,nz)
 k=nz; B(:,:,k)=(dzt(k)*A(:,:,k)+dzt(k-1)*A(:,:,k-1))/(2*0.5*dzw(k))  
 do k=nz-1,2,-1
  B(:,:,k)=(dzt(k)*A(:,:,k)+dzt(k-1)*A(:,:,k-1))/(2*dzw(k))  
 enddo
 B(:,:,1)=(dzt(1)*A(:,:,1))/(2*dzw(1))  
 B(:,:,nz)=B(:,:,nz) + (dzt(nz)*aloc)/(2*0.5*dzw(nz))  
end subroutine tgrid_to_wgrid
