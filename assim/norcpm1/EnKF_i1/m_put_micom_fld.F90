module m_put_micom_fld
use qmpi, only : master
use netcdf
use nfw_mod
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! KAL - This is for the new file type
subroutine put_micom_fld(memfile,fld,iens,cfld,vlevel,tlevel,nx,ny)
   implicit none
   integer, intent(in) :: nx,ny
   integer,                intent(in)  :: iens   ! Ensemble member to read
   real, dimension(nx,ny), intent(in)  :: fld    ! output fld
   character(len=*),       intent(in)  :: memfile! base name of input files
   character(len=8),       intent(in)  :: cfld   ! name of fld
   integer,                intent(in)  :: tlevel ! time level
   integer,                intent(in)  :: vlevel ! vertical level

   integer :: ex, ncid, vFIELD_ID
   integer, allocatable :: ns(:), nc(:)


   !if (master .and. iens.eq.1) then
   !   print *,'Dumping  ',trim(cfld),vlevel
   !endif
   inquire(file=trim(memfile)//'.nc',exist=ex)
     ! Reading the observation file of satellite
     call nfw_open(trim(memfile)//'.nc', or(nf_write,nf_share), ncid)
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
      call nfw_put_vara_double(trim(memfile)//'.nc', ncid, vFIELD_ID, ns, nc, fld)
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
      call nfw_put_vara_double(trim(memfile)//'.nc', ncid, vFIELD_ID, ns, nc, fld)
     endif
     call nfw_close(trim(memfile)//'.nc', ncid)
end subroutine



end module m_put_micom_fld


