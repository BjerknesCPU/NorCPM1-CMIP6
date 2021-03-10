module m_get_micom_grid
contains 
subroutine get_micom_grid(modlon, modlat, depths, mindx, meandx, nx, ny)
   use netcdf
   use nfw_mod

   implicit none
   integer, intent(in) :: nx,ny
   real, dimension(nx,ny), intent(out) :: modlon,modlat,depths
   real,intent(out)   :: mindx,meandx
   integer ncid, x_ID, y_ID, vLON_ID, vLAT_ID, vDEPTH_ID, vPDX_ID, vPDY_ID
   real, dimension(nx,ny):: pdx,pdy

   logical ex

   inquire(file='grid.nc',exist=ex)
   if (ex) then
     ! Reading the grid file
      call nfw_open('grid.nc', nf_nowrite, ncid)
     ! Get dimension id in netcdf file ...
      call nfw_inq_varid('grid.nc', ncid,'plon' ,vLON_ID)
      call nfw_inq_varid('grid.nc', ncid,'plat' ,vLAT_ID)
      call nfw_inq_varid('grid.nc', ncid,'pdepth' ,vDEPTH_ID)
      call nfw_inq_varid('grid.nc', ncid,'pdx' ,vPDX_ID)
      call nfw_inq_varid('grid.nc', ncid,'pdy' ,vPDY_ID)
      call nfw_get_var_double('grid.nc', ncid, vLON_ID, modlon)
      call nfw_get_var_double('grid.nc', ncid, vLAT_ID, modlat)
      call nfw_get_var_double('grid.nc', ncid, vDEPTH_ID, depths)
      call nfw_get_var_double('grid.nc', ncid, vPDX_ID, pdx)
      call nfw_get_var_double('grid.nc', ncid, vPDY_ID, pdy)
      call nfw_close('grid.nc', ncid)
      mindx = min(real(minval(pdx,mask=depths>1.)), real(minval(pdy,mask=depths>1.)))
      meandx=sum  (pdx,mask=depths>1. .and. depths < 1e25) / &
          count(depths>1. .and. depths < 1e25)
   else
      stop 'ERROR: file grid.nc is missing'
   endif
end subroutine  get_micom_grid
end module  m_get_micom_grid
