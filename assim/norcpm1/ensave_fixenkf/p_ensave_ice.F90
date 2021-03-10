! File:          p_ensmean.F90
!
! Created:       Francois counillon
!
! Last modified: 24/05/2014
!
! Purpose:       Create a light and smart average of the ensemble
!
! Description:  
!            This program recompute the kfpla from the model output
!                         ensure that the sum of dp = pb
!            Input is the mem
!
! Modifications:
!
program ensave_ice
use netcdf
use nfw_mod
   implicit none

   integer*4, external :: iargc
   integer imem                  ! ensemble member
   character(len=80) :: oldfile,cbasis, char80
   logical          :: ex
   character(len=8) :: cfld, ctmp
   character(len=3) :: cproc,cmem
   integer          :: tlevel, vlevel, nproc
   integer          :: idm,jdm,kdm
   real, allocatable, dimension(:,:,:)     :: aicen_ave, vicen_ave,Tsfcn_ave
   real, allocatable, dimension(:,:,:)     :: aicentmp, vicentmp,Tsfcntmp,nbtmp
   integer :: ios,ios2,nmem
   integer :: dimids(3)
   integer :: i,j,k
   integer, allocatable :: ns(:), nc(:)
   integer :: ncid, x_ID, y_ID, z_ID
   integer :: vAICEN_ID, vVICEN_ID, vTSFCN_ID,vNB_ID
   integer :: ncid2



   if (iargc()==2 ) then
      call getarg(1,cbasis)
      call getarg(2,ctmp)
      read(ctmp,*) nmem
   else
      print *,'"ensave_ice" -- A light and smart ensave for ice'
      print *
      print *,'usage: '
      print *,'   ensave_ice cbasis ensemble_size'
      print *,'ex:   ensave analysis 30'
      call exit(1)
   endif
   write(cmem,'(i3.3)') 1
   oldfile=trim(cbasis)//cmem//'.nc'
   print *, 'reading dim from file:',oldfile
   inquire(exist=ex,file=trim(oldfile))
   if (.not.ex) then
      write(*,*) 'Can not find '//trim(oldfile)
      stop '(ensave)'
   end if
   ! Reading the restart file
   call nfw_open(trim(oldfile), nf_nowrite, ncid)
   ! Get dimension id in netcdf file ...
   !nb total of data
   call nfw_inq_dimid(trim(oldfile), ncid, 'ni', x_ID)
   call nfw_inq_dimid(trim(oldfile), ncid, 'nj', y_ID)
   call nfw_inq_dimid(trim(oldfile), ncid, 'ncat', z_ID)
   !nb total of track
   call nfw_inq_dimlen(trim(oldfile), ncid, x_ID, idm)
   call nfw_inq_dimlen(trim(oldfile), ncid, y_ID, jdm)
   call nfw_inq_dimlen(trim(oldfile), ncid, z_ID, kdm)
   call nfw_close(trim(oldfile), ncid)
  ! print *, 'The model dimension is :',idm,jdm,kdm
   allocate(aicen_ave (idm,jdm,kdm))
   allocate(vicen_ave(idm,jdm,kdm))
   allocate(Tsfcn_ave(idm,jdm,kdm))
   allocate(aicentmp (idm,jdm,kdm))
   allocate(vicentmp(idm,jdm,kdm))
   allocate(Tsfcntmp(idm,jdm,kdm))
   allocate(nbtmp(idm,jdm,kdm))
   aicen_ave(:,:,:)=0.
   vicen_ave(:,:,:)=0.
   Tsfcn_ave(:,:,:)=0.
   nbtmp(:,:,:)=0.
   !Reading dp 
   allocate(ns(3))
   allocate(nc(3))
   do imem=1,nmem
      write(cmem,'(i3.3)') imem
      oldfile=trim(cbasis)//cmem//'.nc'
      inquire(exist=ex,file=trim(oldfile))
      if (.not.ex) then
         write(*,*) 'Can not find '//oldfile
         stop '(ensave)'
      end if
      call nfw_open(trim(oldfile), nf_nowrite, ncid)
      ns(1)=1
      ns(2)=1
      ns(3)=1
      nc(1)=idm
      nc(2)=jdm
      nc(3)=kdm
      call nfw_inq_varid(trim(oldfile), ncid,'aicen',vAICEN_ID)
      call nfw_get_vara_double(trim(oldfile), ncid, vAICEN_ID, ns, nc, aicentmp)
      call nfw_inq_varid(trim(oldfile), ncid,'vicen',vVICEN_ID)
      call nfw_get_vara_double(trim(oldfile), ncid, vVICEN_ID, ns, nc, vicentmp)
      call nfw_inq_varid(trim(oldfile), ncid,'Tsfcn',vTSFCN_ID)
      call nfw_get_vara_double(trim(oldfile), ncid, vTSFCN_ID, ns, nc, Tsfcntmp)
      ! DP correction
      do j=1,jdm
      do i=1,idm
        !only if not land mask
         do k = 1, kdm
           if (aicentmp(i,j,k)>0.) then  
              Tsfcn_ave(i,j,k)=Tsfcn_ave(i,j,k)+Tsfcntmp(i,j,k)
              nbtmp(i,j,k)=nbtmp(i,j,k)+1
              aicen_ave(i,j,k)=aicen_ave(i,j,k)+aicentmp(i,j,k)
              vicen_ave(i,j,k)=vicen_ave(i,j,k)+vicentmp(i,j,k)
           end if
         end do
      end do
      end do
   enddo !mem
   do j=1,jdm
   do i=1,idm
      do k = 1, kdm
        if (nbtmp(i,j,k).eq.0.0) then 
          Tsfcn_ave(i,j,k)=0.0!Tsfcn_ave(i,j,k)
          aicen_ave(i,j,k)=0.0
          vicen_ave(i,j,k)=0.0
        else 
          Tsfcn_ave(i,j,k)=Tsfcn_ave(i,j,k)/nbtmp(i,j,k)
          aicen_ave(i,j,k)=aicen_ave(i,j,k)/nbtmp(i,j,k)
          vicen_ave(i,j,k)=vicen_ave(i,j,k)/nbtmp(i,j,k)

        end if 
        !aicen_ave(i,j,k)=aicen_ave(i,j,k)/nmem
        !vicen_ave(i,j,k)=vicen_ave(i,j,k)/nmem
      end do
   end do
   end do

   oldfile=trim(cbasis)//'_avg.nc'
   call nfw_create(trim(oldfile), nf_clobber, ncid)
   call nfw_def_dim(trim(oldfile), ncid, 'ni', idm, dimids(1))
   call nfw_def_dim(trim(oldfile), ncid, 'nj', jdm, dimids(2))
   call nfw_def_dim(trim(oldfile), ncid, 'ncat', kdm, dimids(3))
   call nfw_def_var(trim(oldfile), ncid, 'aicen',nf_float, 3, dimids, vAICEN_ID)
   call nfw_def_var(trim(oldfile), ncid, 'vicen',nf_float, 3, dimids, vVICEN_ID)
   call nfw_def_var(trim(oldfile), ncid, 'nbmem',nf_float, 3, dimids, vNB_ID)
   call nfw_def_var(trim(oldfile), ncid, 'Tsfcn',nf_float, 3, dimids, vTSFCN_ID)
   call nfw_enddef(trim(oldfile), ncid)
   call nfw_put_var_double(trim(oldfile), ncid, vAICEN_ID, aicen_ave(:,:,:))
   call nfw_put_var_double(trim(oldfile), ncid, vVICEN_ID, vicen_ave(:,:,:))
   call nfw_put_var_double(trim(oldfile), ncid, vTSFCN_ID, Tsfcn_ave(:,:,:))
   call nfw_put_var_double(trim(oldfile), ncid, vNB_ID, nbtmp(:,:,:))   
   call nfw_close(trim(oldfile), ncid)
end program
