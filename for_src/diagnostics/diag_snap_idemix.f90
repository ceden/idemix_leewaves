



subroutine init_snap_cdf_idemix
!=======================================================================
!     initialize NetCDF snapshot file for idemix variables
!=======================================================================
 use main_module   
 use idemix_module   
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 integer :: lon_tdim,lon_udim,itimedim
 integer :: lat_tdim,lat_udim,id,z_tdim,z_udim
 integer :: dims(4)
 character :: name*32, unit*16
 real*8, parameter :: spval = -1.0d33

 if (my_pe==0) print'(a)',' preparing file pyOM_idemix.cdf'

 call def_grid_cdf('pyOM_idemix.cdf')

 if (my_pe==0) then
      iret=nf_open('pyOM_idemix.cdf',NF_WRITE, ncid)
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
      call ncredf(ncid, iret)
      iret=nf_inq_dimid(ncid,'xt',lon_tdim)
      iret=nf_inq_dimid(ncid,'xu',lon_udim)
      iret=nf_inq_dimid(ncid,'yt',lat_tdim)
      iret=nf_inq_dimid(ncid,'yu',lat_udim)
      iret=nf_inq_dimid(ncid,'zt',z_tdim)
      iret=nf_inq_dimid(ncid,'zu',z_udim)
      iret=nf_inq_dimid(ncid,'Time',itimedim)

      dims = (/Lon_tdim,lat_tdim,z_udim,itimedim/)
      id  = ncvdef (ncid,'Nsqr', NCFLOAT,4,dims,iret)
      name = 'Square of stability frequency'; unit = '1/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_udim,lat_tdim,z_tdim,itimedim/)
      id  = ncvdef (ncid,'U', NCFLOAT,4,dims,iret)
      name = 'Mean flow'; unit = 'm/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'c0', NCFLOAT,4,dims,iret)
      name = 'vertical IW group velocity'; unit = 'm/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)  ! mu E^2 = E_t [m^2/s^3]  ->  mu [ m^2/s^3 (s^4/m^4) =  s/m^2 ]
      id  = ncvdef (ncid,'mu', NCFLOAT,4,dims,iret)
      name = 'Parameter'; unit = 'm/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
      id  = ncvdef (ncid,'cstar', NCFLOAT,3,dims,iret)
      name = 'modal gravity wave speed'; unit = 'm/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'v0', NCFLOAT,4,dims,iret)
      name = 'hor. IW group velocity'; unit = 'm/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'iw_diss', NCFLOAT,4,dims,iret)
      name = 'IW dissipation'; unit = 'm^2/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

  
        dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
        id  = ncvdef (ncid,'lee_flux', NCFLOAT,3,dims,iret)
        name = 'IW bottom forcing'; unit = 'm^3/s^3'
        call dvcdf(ncid,id,name,32,unit,16,spval)
       
        dims = (/Lon_udim,lat_tdim,iTimedim,1/)
        id  = ncvdef (ncid,'lee_stress_u', NCFLOAT,3,dims,iret)
        name = 'IW bottom stress'; unit = 'm^2/s^2'
        call dvcdf(ncid,id,name,32,unit,16,spval)
        
        dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
        id  = ncvdef (ncid,'inv_Fr', NCFLOAT,3,dims,iret)
        name = 'inverse Froude number'; unit = ' '
        call dvcdf(ncid,id,name,32,unit,16,spval)
    
        dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
        id  = ncvdef (ncid,'E_lee_p', NCFLOAT,4,dims,iret)
        name = 'Internal wave energy'; unit = 'm^2/s^2'
        call dvcdf(ncid,id,name,32,unit,16,spval) 

        dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
        id  = ncvdef (ncid,'E_lee_m', NCFLOAT,4,dims,iret)
        name = 'Internal wave energy'; unit = 'm^2/s^2'
        call dvcdf(ncid,id,name,32,unit,16,spval) 
      
        dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
        id  = ncvdef (ncid,'Tc_lee', NCFLOAT,4,dims,iret)
        name = 'inverse time scale'; unit = '1/s'
        call dvcdf(ncid,id,name,32,unit,16,spval) 
        
        dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
        id  = ncvdef (ncid,'Lambda_0', NCFLOAT,4,dims,iret)
        name = 'Parameter'; unit = ' '
        call dvcdf(ncid,id,name,32,unit,16,spval) 
        
        dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
        id  = ncvdef (ncid,'iw_diss_lee', NCFLOAT,4,dims,iret)
        name = 'Dissipation of wave energy'; unit = 'm^2/s^3'
        call dvcdf(ncid,id,name,32,unit,16,spval) 
      
      
       
      call ncclos (ncid, iret)
 endif
   

 if (my_pe==0) print'(a)',' done '
 call fortran_barrier

end subroutine init_snap_cdf_idemix


subroutine diag_snap_idemix
!=======================================================================
!     write to NetCDF snapshot file
!=======================================================================
 use main_module   
 use idemix_module   
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 real*8 :: bloc(nx,ny),fxa
 integer :: itdimid,ilen,itimeid,id,k
 real*8, parameter :: spval = -1.0d33

 if (my_pe == 0 ) then
   iret=nf_open('pyOM_idemix.cdf',NF_WRITE,ncid)
   print'(a)',' writing to file pyOM_idemix.cdf'
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   iret=nf_inq_dimid(ncid,'Time',itdimid)
   iret=nf_inq_dimlen(ncid, itdimid,ilen)
   ilen=ilen+1
   fxa = itt*dt_tracer/86400.0
   iret=nf_inq_varid(ncid,'Time',itimeid)
   iret= nf_put_vara_double(ncid,itimeid,ilen,1,fxa)
   call ncclos (ncid, iret)
 endif
 call fortran_barrier

 if (my_pe == 0 ) then
    iret=nf_open('pyOM_idemix.cdf',NF_WRITE,ncid)
    iret=nf_inq_dimid(ncid,'Time',itdimid)
    iret=nf_inq_dimlen(ncid, itdimid,ilen)
 endif

 ! cstar
 bloc(is_pe:ie_pe,js_pe:je_pe) = cstar(is_pe:ie_pe,js_pe:je_pe)
 where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe==0) then
    iret=nf_inq_varid(ncid,'cstar',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
 endif

 
 
  bloc(is_pe:ie_pe,js_pe:je_pe) = flux_lee_tot(is_pe:ie_pe,js_pe:je_pe)
  where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'lee_flux',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif

  bloc(is_pe:ie_pe,js_pe:je_pe) = lee_stress_u(is_pe:ie_pe,js_pe:je_pe)
  where( maskU(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'lee_stress_u',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif
  
  bloc(is_pe:ie_pe,js_pe:je_pe) = inv_fr(is_pe:ie_pe,js_pe:je_pe)
  where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'inv_Fr',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif 
  
 

 ! 3D fields

 do k=1,nz

   ! stability frequency
   bloc(is_pe:ie_pe,js_pe:je_pe) = Nsqr(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'Nsqr',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   bloc(is_pe:ie_pe,js_pe:je_pe) = alpha_c(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'mu',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   bloc(is_pe:ie_pe,js_pe:je_pe) = u(is_pe:ie_pe,js_pe:je_pe,k,tau)
   where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'U',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

 
   ! Internal wave dissipation
   bloc(is_pe:ie_pe,js_pe:je_pe) = iw_diss(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'iw_diss',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! Internal wave vertical group velocity
   bloc(is_pe:ie_pe,js_pe:je_pe) = c0(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'c0',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! Internal wave horizontal group velocity
   bloc(is_pe:ie_pe,js_pe:je_pe) = v0(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'v0',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif



    
      bloc(is_pe:ie_pe,js_pe:je_pe) = E_lee_p(is_pe:ie_pe,js_pe:je_pe,k,tau)
      where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
      call pe0_recv_2D(nx,ny,bloc)
      if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'E_lee_p',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
      endif

      bloc(is_pe:ie_pe,js_pe:je_pe) = E_lee_m(is_pe:ie_pe,js_pe:je_pe,k,tau)
      where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
      call pe0_recv_2D(nx,ny,bloc)
      if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'E_lee_m',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
      endif
    
      bloc(is_pe:ie_pe,js_pe:je_pe) = Tc_lee(is_pe:ie_pe,js_pe:je_pe,k)
      where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
      call pe0_recv_2D(nx,ny,bloc)
      if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'Tc_lee',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
      endif      
 
      bloc(is_pe:ie_pe,js_pe:je_pe) = Lambda_0(is_pe:ie_pe,js_pe:je_pe,k)
      where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
      call pe0_recv_2D(nx,ny,bloc)
      if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'Lambda_0',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
      endif      
  
  
      bloc(is_pe:ie_pe,js_pe:je_pe) = iw_diss_lee(is_pe:ie_pe,js_pe:je_pe,k)
      where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
      call pe0_recv_2D(nx,ny,bloc)
      if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'iw_diss_lee',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
      endif      
  
           
  
    
 enddo

 if (my_pe==0)   iret = nf_close (ncid)

end subroutine diag_snap_idemix








subroutine def_grid_cdf(filename)
!=======================================================================
!      Define standard grid in netcdf file
!=======================================================================
 use main_module   
 implicit none
 include "netcdf.inc"
 character*(*) filename
 integer :: ncid,iret,n
 integer :: lon_tdim,lon_udim,itimedim
 integer :: lon_tid,lon_uid,itimeid
 integer :: lat_tdim,lat_udim,lat_uid,lat_tid
 integer :: z_tdim,z_udim,z_tid,z_uid
 character :: name*24, unit*16

 if (my_pe==0) then
      iret = nf_create (filename, nf_clobber, ncid)
      if (iret /= 0) call halt_stop('NETCDF:'//nf_strerror(iret))
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
!     dimensions
      lon_tdim  = ncddef(ncid, 'xt', nx , iret)
      Lon_udim  = ncddef(ncid, 'xu', nx , iret)
      lat_tdim  = ncddef(ncid, 'yt', ny , iret)
      Lat_udim  = ncddef(ncid, 'yu', ny , iret)
      iTimedim  = ncddef(ncid, 'Time', nf_unlimited, iret)
!     grid variables
      Lon_tid  = ncvdef (ncid,'xt',NCFLOAT,1,lon_tdim,iret)
      Lon_uid  = ncvdef (ncid,'xu',NCFLOAT,1,lon_udim,iret)
      Lat_tid  = ncvdef (ncid,'yt',NCFLOAT,1,lat_tdim,iret)
      Lat_uid  = ncvdef (ncid,'yu',NCFLOAT,1,lat_udim,iret)
      itimeid  = ncvdef (ncid,'Time', NCFLOAT,1,itimedim,iret)
!     attributes of the grid
      if (coord_degree) then
       name = 'Longitude on T grid     '; unit = 'degrees E'
       call ncaptc(ncid, Lon_tid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lon_tid, 'units',     NCCHAR, 16, unit, iret) 
       name = 'Longitude on U grid     '; unit = 'degrees E'
       call ncaptc(ncid, Lon_uid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lon_uid, 'units',     NCCHAR, 16, unit, iret) 
       name = 'Latitude on T grid     '; unit = 'degrees N'
       call ncaptc(ncid, Lat_tid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lat_tid, 'units',     NCCHAR, 16, unit, iret) 
       name = 'Latitude on U grid     '; unit = 'degrees N'
       call ncaptc(ncid, Lat_uid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lat_uid, 'units',     NCCHAR, 16, unit, iret) 
      else
       name = 'zonal axis T grid     '; unit = 'km'
       call ncaptc(ncid, Lon_tid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lon_tid, 'units',     NCCHAR, 16, unit, iret) 
       name = 'zonal axis U grid     '; unit = 'km'
       call ncaptc(ncid, Lon_uid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lon_uid, 'units',     NCCHAR, 16, unit, iret) 
       name = 'meridional axis T grid'; unit = 'km'
       call ncaptc(ncid, Lat_tid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lat_tid, 'units',     NCCHAR, 16, unit, iret) 
       name = 'meridional axis U grid'; unit = 'km'
       call ncaptc(ncid, Lat_uid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lat_uid, 'units',     NCCHAR, 16, unit, iret) 
      endif

      name = 'Time '; unit = 'days'
      call ncaptc(ncid, itimeid, 'long_name', NCCHAR, 24, name, iret) 
      call ncaptc(ncid, itimeid, 'units',     NCCHAR, 16, unit, iret) 
      call ncaptc(ncid, iTimeid,'time_origin',NCCHAR, 20,'01-JAN-1900 00:00:00', iret)

      z_tdim    = ncddef(ncid, 'zt',  nz, iret)
      z_udim    = ncddef(ncid, 'zu',  nz, iret)
      z_tid  = ncvdef (ncid,'zt', NCFLOAT,1,z_tdim,iret)
      z_uid  = ncvdef (ncid,'zu', NCFLOAT,1,z_udim,iret)
      name = 'Height on T grid     '; unit = 'm'
      call ncaptc(ncid, z_tid, 'long_name', NCCHAR, 24, name, iret) 
      call ncaptc(ncid, z_tid, 'units',     NCCHAR, 16, unit, iret) 
      name = 'Height on U grid     '; unit = 'm'
      call ncaptc(ncid, z_uid, 'long_name', NCCHAR, 24, name, iret) 
      call ncaptc(ncid, z_uid, 'units',     NCCHAR, 16, unit, iret) 

      call ncendf(ncid, iret)
      iret= nf_put_vara_double(ncid,z_tid,1,nz ,zt)
      iret= nf_put_vara_double(ncid,z_uid,1,nz ,zw)
      call ncclos (ncid, iret)
 endif


 do n=0,n_pes-1
   call fortran_barrier
   if (my_pe==n) then
      iret=nf_open(filename,NF_WRITE,ncid)
      iret=nf_inq_varid(ncid,'xt',lon_tid)
      iret=nf_inq_varid(ncid,'xu',lon_uid)
      iret=nf_inq_varid(ncid,'yt',lat_tid)
      iret=nf_inq_varid(ncid,'yu',lat_uid)
      if (coord_degree) then
       iret= nf_put_vara_double(ncid,lon_Tid,is_pe,ie_pe-is_pe+1 ,xt(is_pe:ie_pe))
       iret= nf_put_vara_double(ncid,lon_uid,is_pe,ie_pe-is_pe+1 ,xu(is_pe:ie_pe))
       iret= nf_put_vara_double(ncid,lat_Tid,js_pe,je_pe-js_pe+1 ,yt(js_pe:je_pe))
       iret= nf_put_vara_double(ncid,lat_uid,js_pe,je_pe-js_pe+1 ,yu(js_pe:je_pe))
      else
       iret= nf_put_vara_double(ncid,lon_Tid,is_pe,ie_pe-is_pe+1 ,xt(is_pe:ie_pe)/1e3)
       iret= nf_put_vara_double(ncid,lon_uid,is_pe,ie_pe-is_pe+1 ,xu(is_pe:ie_pe)/1e3)
       iret= nf_put_vara_double(ncid,lat_Tid,js_pe,je_pe-js_pe+1 ,yt(js_pe:je_pe)/1e3)
       iret= nf_put_vara_double(ncid,lat_uid,js_pe,je_pe-js_pe+1 ,yu(js_pe:je_pe)/1e3)
      endif
      call ncclos (ncid, iret)
   endif
   call fortran_barrier
 enddo
end subroutine def_grid_cdf




subroutine dvcdf(ncid,ivarid,name,iname,unit,iunit,spval)
!=======================================================================
!     define some standard attributes of variable ivarid in NetCDF file ncid 
!=======================================================================
 implicit none
 integer ncid,ivarid,iname,iunit,iret
 character (len=*) :: name, unit
 real*8 :: spval, vv
 include "netcdf.inc"
 vv=spval
 call ncaptc(ncid,ivarid, 'long_name', NCCHAR,iname , name, iret) 
 if (iret.ne.0) print*,nf_strerror(iret)
 call ncaptc(ncid,ivarid, 'units',     NCCHAR,iunit, unit, iret) 
 if (iret.ne.0) print*,nf_strerror(iret)
 call ncapt (ncid,ivarid, 'missing_value',NCDOUBLE,1,vv,iret)
 if (iret.ne.0) print*,nf_strerror(iret)
 call ncapt (ncid,ivarid, '_FillValue', NCDOUBLE, 1,vv, iret)
 if (iret.ne.0) print*,nf_strerror(iret)
end subroutine dvcdf




