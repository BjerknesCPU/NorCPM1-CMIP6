module m_get_micom_fld
use netcdf
use nfw_mod
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! KAL -- This routine reads one of the fields from the model, specified  
! KAL -- by name, vertical level and time level
! KAL -- This routine is really only effective for the new restart files.
subroutine get_micom_fld_new(memfile,fld,cfld,vlevel,tlevel,nx,ny)
#if defined (QMPI)
   use qmpi, only : qmpi_proc_num
#else
   use qmpi_fake
#endif
   implicit none
   integer,      intent(in)            :: nx,ny  ! Grid dimension
   real, dimension(nx,ny), intent(out) :: fld    ! output fld
   character(len=*), intent(in)        :: memfile! base name of input files
   character(len=*), intent(in)        :: cfld   ! name of fld
   integer, intent(in)                 :: tlevel ! time level
   integer, intent(in)                 :: vlevel ! vertical level

   real, dimension(nx,ny) :: readfld
   integer :: ex, ncid, vFIELD_ID
   integer, allocatable :: ns(:), nc(:)
   inquire(file=trim(memfile)//'.nc',exist=ex)
   if (ex) then
     ! Reading the observation file of satellite
     call nfw_open(trim(memfile)//'.nc', nf_nowrite, ncid)
     call nfw_inq_varid(trim(memfile)//'.nc', ncid,trim(cfld),vFIELD_ID)
     if (vlevel==0) then
      allocate(ns(3))
      allocate(nc(3))
      ns(1)=1
      ns(2)=1
      ns(3)=1
      nc(1)=nx
      nc(2)=ny
      nc(3)=1
      call nfw_get_vara_double(trim(memfile)//'.nc', ncid, vFIELD_ID, ns, nc, readfld)
      fld=readfld(:,:)
     else
      allocate(ns(4))
      allocate(nc(4))
      ns(1)=1
      ns(2)=1
      ns(3)=vlevel
      ns(4)=1
      nc(1)=nx
      nc(2)=ny
      nc(3)=1
      nc(4)=1
      call nfw_get_vara_double(trim(memfile)//'.nc', ncid, vFIELD_ID, ns, nc, readfld)
      fld=readfld(:,:)
     endif
     call nfw_close(trim(memfile)//'.nc', ncid)
  else
     print *, 'ERROR: forecast file is missing '//trim(memfile)//'.nc'
     stop
  endif
end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! KAL -- This routine reads one of the fields from the model, specified  
! KAL -- by name and time level
! KAL -- This routine is really only effective for vertical interpolation
subroutine get_micom_fld(memfile,fld,cfld,tlevel,nx,ny,nz)
#if defined (QMPI)
   use qmpi, only : qmpi_proc_num
#else
   use qmpi_fake
#endif
   implicit none
   integer,      intent(in)            :: nx,ny,nz ! Grid dimension
   real, dimension(nx,ny,nz), intent(out) :: fld   ! output fld
   character(len=*), intent(in)        :: memfile  ! base name of input files
   character(len=*), intent(in)        :: cfld     ! name of fld
   integer, intent(in)                 :: tlevel   ! time level

   real, dimension(nx,ny,nz) :: readfld
   integer :: ex, ncid, vFIELD_ID
   integer, allocatable :: ns(:), nc(:)
   inquire(file=trim(memfile)//'.nc',exist=ex)
   if (ex) then
     ! Reading the observation file of satellite
     call nfw_open(trim(memfile)//'.nc', nf_nowrite, ncid)
     call nfw_inq_varid(trim(memfile)//'.nc', ncid,trim(cfld),vFIELD_ID)

     allocate(ns(4))
     allocate(nc(4))
     ns(1)=1
     ns(2)=1
     ns(3)=1
     ns(4)=1
     nc(1)=nx
     nc(2)=ny
     nc(3)=nz
     nc(4)=1
     call nfw_get_vara_double(trim(memfile)//'.nc', ncid, vFIELD_ID, ns, nc, readfld)
     fld=readfld(:,:,:)
     call nfw_close(trim(memfile)//'.nc', ncid)
  else
     print *, 'ERROR: forecast file is missing '//trim(memfile)//'.nc'
     stop
  endif
end subroutine get_micom_fld

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! KAL -- This routine reads one dimension of the fields from the model outputs, specified
! KAL -- by name and vertical level
! KAL -- This routine is really only effective for vertical interpolation  
subroutine get_micom_fld_1d(memfile,fld,cfld,nz)
#if defined (QMPI)
   use qmpi, only : qmpi_proc_num
#else
   use qmpi_fake
#endif
   implicit none
   integer,      intent(in)            :: nz ! Grid dimension
   real, dimension(nz), intent(out)    :: fld   ! output fld 
   character(len=*), intent(in)        :: memfile  ! base name of input files 
   character(len=*), intent(in)        :: cfld     ! name of fld

   real, dimension(nz) :: readfld
   integer :: ex, ncid, vFIELD_ID

   inquire(file=trim(memfile)//'.nc',exist=ex)
   if (ex) then
      ! Reading data
      call nfw_open(trim(memfile)//'.nc', nf_nowrite, ncid)
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,trim(cfld),vFIELD_ID)

      call nfw_get_var_double(trim(memfile)//'.nc', ncid, vFIELD_ID, readfld)
      fld = readfld 
      call nfw_close(trim(memfile)//'.nc', ncid)
   else
      print *, 'ERROR: file is missing '//trim(memfile)//'.nc'
      stop
   endif
end subroutine get_micom_fld_1d

end module m_get_micom_fld


