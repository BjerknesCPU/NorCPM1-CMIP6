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
program ensave
#if defined(QMPI)
  use qmpi
#else
  use qmpi_fake
#endif
  use netcdf
  use nfw_mod
  use distribute
  implicit none

  integer*4, external :: iargc

  integer imem ! ensemble member
  character(len=80):: oldfile, cbasis, ctmp
  logical          :: ex
  character(len=3) :: cmem
  integer          :: idm,jdm,kdm
  integer :: nmem
  integer :: dimids(3)
  integer :: i,j,k
  integer :: ncid, x_ID, y_ID, z_ID
  integer :: vTEM_ID,vSAL_ID,vDP_ID,vNB_ID
  integer, allocatable :: ns(:), nc(:)
  real, allocatable, dimension(:,:,:) :: dp_ave,saln_ave,temp_ave
  real, allocatable, dimension(:,:,:) :: dptmp,salntmp,temptmp,nbtmp
  real, allocatable, dimension(:,:,:) :: glb_dp_ave,glb_saln_ave,glb_temp_ave,nb

#if defined(QMPI)
  call start_mpi()
#endif

  if (iargc()==2 ) then
     call getarg(1,cbasis)
     call getarg(2,ctmp)
     read(ctmp,*) nmem
  else
     print *,'"ensave" -- A light and smart ensave'
     print *
     print *,'usage: '
     print *,'   ensave cbasis ensemble_size'
     print *,'ex:   ./ensave analysis 30'
     print *,'ex:   mpirun -n 30 ./ensave analysis 30'
     call exit(1)
  endif

  ! Get dimension in the first netcdf file ...
  write(cmem,'(i3.3)') 1
  oldfile=trim(cbasis)//cmem//'.nc'
  if (master) then
     print *, 'calculating ensemble mean...'
  end if
  inquire(exist=ex,file=trim(oldfile))
  if (.not.ex) then
     write(*,*) 'Can not find '//trim(oldfile)
     stop '(ensave)'
  end if
  call nfw_open(trim(oldfile), nf_nowrite, ncid)
  !nb total of data
  call nfw_inq_dimid(trim(oldfile), ncid, 'x', x_ID)
  call nfw_inq_dimid(trim(oldfile), ncid, 'y', y_ID)
  call nfw_inq_dimid(trim(oldfile), ncid, 'kk', z_ID)
  !nb total of track
  call nfw_inq_dimlen(trim(oldfile), ncid, x_ID, idm)
  call nfw_inq_dimlen(trim(oldfile), ncid, y_ID, jdm)
  call nfw_inq_dimlen(trim(oldfile), ncid, z_ID, kdm)
  call nfw_close(trim(oldfile), ncid)
!  if (master) then
!     print *, 'The model dimension is :',idm,jdm,kdm
!  end if

  call distribute_iterations(nmem)
#if defined(QMPI)
  call barrier()
#endif

  !Allocate tmp variables
  allocate(dp_ave(idm,jdm,kdm),dptmp(idm,jdm,kdm))
  allocate(saln_ave(idm,jdm,kdm),salntmp(idm,jdm,kdm))
  allocate(temp_ave(idm,jdm,kdm),temptmp(idm,jdm,kdm))
  allocate(nbtmp(idm,jdm,kdm),ns(4), nc(4))
  temp_ave=0
  saln_ave=0
  dp_ave=0
  nbtmp=0
  do imem = my_first_iteration, my_last_iteration
     write(cmem,'(i3.3)') imem
     oldfile=trim(cbasis)//cmem//'.nc'
     !print *, trim(oldfile)
     inquire(exist=ex,file=trim(oldfile))
     if (.not.ex) then
        write(*,*) 'Can not find '//oldfile
        stop '(ensave)'
     end if
     call nfw_open(trim(oldfile), nf_nowrite, ncid)
     ns(1)=1
     ns(2)=1
     ns(3)=1
     ns(4)=1
     nc(1)=idm
     nc(2)=jdm
     nc(3)=kdm
     nc(4)=1
     !Get variable
     call nfw_inq_varid(trim(oldfile), ncid,'dp',vDP_ID)
     call nfw_get_vara_double(trim(oldfile), ncid, vDP_ID, ns, nc, dptmp)
     call nfw_inq_varid(trim(oldfile), ncid,'temp',vTEM_ID)
     call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTEM_ID, ns, nc, temptmp)
     call nfw_inq_varid(trim(oldfile), ncid,'saln',vSAL_ID)
     call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vSAL_ID, ns, nc, salntmp)
     ! DP correction
     do j=1,jdm
        do i=1,idm
           !only if not land mask
           if (dptmp(i,j,1)<100000000000.) then
              do k = 1, kdm
                 if (dptmp(i,j,k)>0.) then
                    temp_ave(i,j,k)=temp_ave(i,j,k)+temptmp(i,j,k)
                    saln_ave(i,j,k)=saln_ave(i,j,k)+salntmp(i,j,k)
                    nbtmp(i,j,k)=nbtmp(i,j,k)+1
                 end if
                 dp_ave(i,j,k)=dp_ave(i,j,k)+dptmp(i,j,k)
              end do
           end if
        end do
     end do
     call nfw_close(trim(oldfile), ncid)
  end do
  deallocate(dptmp,temptmp,salntmp)
#if defined(QMPI)
  call barrier()
#endif

  ! Collect local variables to global variables
  allocate(nb(idm,jdm,kdm))
  allocate(glb_dp_ave(idm,jdm,kdm))
  allocate(glb_temp_ave(idm,jdm,kdm))
  allocate(glb_saln_ave(idm,jdm,kdm))
#if defined(QMPI)
  call mpi_reduce(nbtmp,nb,idm*jdm*kdm,mpi_double,mpi_sum,0,mpi_comm_world,ierr)
  call mpi_reduce(temp_ave,glb_temp_ave,idm*jdm*kdm,mpi_double,mpi_sum,0,mpi_comm_world,ierr)
  call mpi_reduce(saln_ave,glb_saln_ave,idm*jdm*kdm,mpi_double,mpi_sum,0,mpi_comm_world,ierr)
  call mpi_reduce(dp_ave,glb_dp_ave,idm*jdm*kdm,mpi_double,mpi_sum,0,mpi_comm_world,ierr)
#else
  nb=nbtmp
  glb_temp_ave=temp_ave
  glb_saln_ave=saln_ave
  glb_dp_ave=dp_ave
#endif
  if (master) then
     do j=1,jdm
        do i=1,idm
           do k = 1,kdm
              glb_temp_ave(i,j,k)=glb_temp_ave(i,j,k)/nb(i,j,k)
              glb_saln_ave(i,j,k)=glb_saln_ave(i,j,k)/nb(i,j,k)
           end do
        end do
     end do
     glb_dp_ave=glb_dp_ave/nmem

     !Write variables to a netcdf file...
     oldfile=trim(cbasis)//'_avg.nc'
     call nfw_create(trim(oldfile), nf_clobber, ncid)

     call nfw_def_dim(trim(oldfile), ncid, 'x', idm, dimids(1))
     call nfw_def_dim(trim(oldfile), ncid, 'y', jdm, dimids(2))
     call nfw_def_dim(trim(oldfile), ncid, 'z', kdm, dimids(3))
     call nfw_def_var(trim(oldfile), ncid, 'temp',nf_float, 3, dimids, vTEM_ID)
     call nfw_def_var(trim(oldfile), ncid, 'saln',nf_float, 3, dimids, vSAL_ID)
     call nfw_def_var(trim(oldfile), ncid, 'nbmem',nf_float, 3, dimids, vNB_ID)
     call nfw_def_var(trim(oldfile), ncid, 'dp',nf_float, 3, dimids, vDP_ID)
     call nfw_enddef(trim(oldfile), ncid)
     
     call nfw_put_var_double(trim(oldfile), ncid, vTEM_ID, glb_temp_ave)
     call nfw_put_var_double(trim(oldfile), ncid, vSAL_ID, glb_saln_ave)
     call nfw_put_var_double(trim(oldfile), ncid, vDP_ID, glb_dp_ave)
     call nfw_put_var_double(trim(oldfile), ncid, vNB_ID, nb)
     call nfw_close(trim(oldfile), ncid)
  end if
#if defined(QMPI)
  call barrier()
  call stop_mpi()
#endif
end program ensave
