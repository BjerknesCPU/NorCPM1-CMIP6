module m_get_micom_dim
contains 
subroutine get_micom_dim(nx,ny,nz)
   use netcdf
   use nfw_mod

   implicit none
   integer, intent(out) :: nx,ny,nz
   integer ncid, x_ID, y_ID, z_ID

   logical ex

   inquire(file='forecast001.nc',exist=ex)

   if (ex) then
     ! Reading the grid file
      call nfw_open('forecast001.nc', nf_nowrite, ncid)
     ! Get dimension id in netcdf file ...
     call nfw_inq_dimid('forecast001.nc', ncid, 'x', x_ID)
     call nfw_inq_dimid('forecast001.nc', ncid, 'y', y_ID)
     call nfw_inq_dimid('forecast001.nc', ncid, 'kk', z_ID)
     !Get the dimension
     call nfw_inq_dimlen('forecast001.nc', ncid, x_ID, nx)
     call nfw_inq_dimlen('forecast001.nc', ncid, y_ID, ny)
     call nfw_inq_dimlen('forecast001.nc', ncid, z_ID, nz)
   else
      stop 'ERROR: file forecast001.nc is missing'
   endif
end subroutine  get_micom_dim

subroutine get_climato_dim(nx,ny,nz)
  use netcdf
  use nfw_mod

   implicit none
   integer, intent(out) :: nx,ny,nz
   integer ncid, x_ID, y_ID, z_ID
   
   logical ex
   
   inquire(file='mean_mod.nc',exist=ex)

   if (ex) then
      ! Reading the grid file
      call nfw_open('mean_mod.nc', nf_nowrite, ncid)
      ! Get dimension id in netcdf file ...
      call nfw_inq_dimid('mean_mod.nc', ncid, 'x', x_ID)
      call nfw_inq_dimid('mean_mod.nc', ncid, 'y', y_ID)
      call nfw_inq_dimid('mean_mod.nc', ncid, 'depth', z_ID)
      !Get the dimension
      call nfw_inq_dimlen('mean_mod.nc', ncid, x_ID, nx)
      call nfw_inq_dimlen('mean_mod.nc', ncid, y_ID, ny)
      call nfw_inq_dimlen('mean_mod.nc', ncid, z_ID, nz)
   else
      stop 'ERROR: file mean_mod.nc is missing'
   endif
 end subroutine  get_climato_dim
end module  m_get_micom_dim
