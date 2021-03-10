module m_get_cice_dim
contains 
subroutine get_cice_dim(ncat,ikdm,skdm)
   use netcdf
   use nfw_mod

   implicit none
   integer, intent(out) :: ncat,ikdm,skdm
   integer :: ncid, ncat_ID, ikdm_ID, skdm_ID

   logical ex

   inquire(file='forecast_ice001.nc',exist=ex)
   if (ex) then
     ! Reading the grid file
      call nfw_open('forecast_ice001.nc', nf_nowrite, ncid)
     ! Get dimension id in netcdf file ...
     call nfw_inq_dimid('forecast_ice001.nc', ncid, 'ncat', ncat_ID)
     call nfw_inq_dimid('forecast_ice001.nc', ncid, 'ntilyr', ikdm_ID)
     call nfw_inq_dimid('forecast_ice001.nc', ncid, 'ntslyr', skdm_ID)
     !Get the dimension
     call nfw_inq_dimlen('forecast_ice001.nc', ncid, ncat_ID, ncat)
     call nfw_inq_dimlen('forecast_ice001.nc', ncid, ikdm_ID, ikdm)
     call nfw_inq_dimlen('forecast_ice001.nc', ncid, skdm_ID, skdm)
   else
      stop 'ERROR: file forecast_ice001.nc is missing'
   endif
end subroutine  get_cice_dim
end module  m_get_cice_dim
