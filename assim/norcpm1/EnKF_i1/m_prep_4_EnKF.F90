! File:          m_prep_4_EnKF.F90
!
! Created:       ???
!
! Last modified: 29/06/2010
!
! Purpose:       Calculation of HA ("S")
!
! Description:   Calculates HA by going sequentially through each data type.
!
! Modifications:
!                05/08/2016 Yiguo WANG: add SSH data type
!                22/07/2016 Yiguo WANG: add ANOMALY option for profile data
!                09/03/2016 Yiguo WANG: replace get_S() by get_S_spline()
!                29/07/2010 PS:
!                  - merged insitu_QC() with generic obs_QC(). Moved
!                    insitu_writeforecast() to the point after QC.
!                29/06/2010 PS:
!                  - added generic observation QC: increase the observation
!                    error when observation and ensemble mean are much too far
!                    away than expected
!                Prior history:
!                  Not documented.

module m_prep_4_EnKF

  integer, parameter, private :: STRLEN = 512

  private read_mean_ssh

contains

  ! This subroutine uses the observation and ensembles from the model
  ! to prepare the input to the EnKF analysis scheme.
  ! The output from this routine is used directly in the global analysis
  ! while the output has to be run through a "filter" to be used in the
  ! local analysis scheme.

  ! S = HA     (ensemble observation anomalies)
  ! d = d - Hx (innovations) 
  !
  ! S is calculated in two steps:
  ! 1. S = HE
  ! 2. S = S - repmat(s, 1, m), 
  !    where s = mean(S')';
  ! Note that in reality (with HYCOM/MICOM) H is different for each member... 
  ! So that HX must be read "HX" rather than "H * X".
  !
  subroutine prep_4_EnKF(nrens, d, S, depths, meandx, nx, ny, nz)
#if defined (QMPI)
    use qmpi, only : master, stop_mpi
#else
    use qmpi_fake, only : master, stop_mpi
#endif
    use mod_measurement
    use m_obs
    use m_Generate_element_Si
    use m_get_micom_fld
    use m_parameters
    implicit none

    integer, intent(in) :: nx, ny, nz ! Model size
    integer, intent(in) :: nrens ! Size of ensemble
    real, intent(in) :: depths(nx, ny)
    real, intent(in) :: meandx ! mean grid size
    real, intent(inout) :: d(nobs)
    real, intent(inout) :: S(nobs, nrens)

    real :: x(nobs)

    integer :: i, j, m, iens
    real*4, dimension(nx,ny) :: fldr4
    real :: readfld(nx, ny)
    real :: C(nobs)
    
    ! hard-coded for now
    integer, parameter :: drnx = 152, drny = 132
    real*4, dimension(drnx, drny) :: modzon, modmer

    integer :: reclSLA, ios, reclDRIFT
    character*3 :: cmem
    character*2 :: day

    logical :: ex

    character(STRLEN) :: fname
    integer :: iuobs

    ! FANF: For track assim we launch m_Generate_Si for each day (t=1:Wd)
    ! which fills in S at the approriate indices.
    ! Wd is the assimilation cycle (i.e. 7 days)
    !
    integer :: Wd, t
    integer :: tlevel
    real :: field2(nx, ny), field3(nx, ny) ! auxiliary fields (e.g. mean SSH, 
                                           ! field bias estimate etc.)
    integer :: nthisobs, thisobs(nobs)

    if (any(obs(:) % id == 'TSLA ')) then
       Wd = 6
    else
       Wd = 0
    endif

    ! security check
    !
    if (any(obs(:) % id == 'SSH  ') .or. any(obs(:) % id == 'SLA  ')) then
       if (any(obs(:) % id == 'SLA  ')) then
          inquire(exist = ex, file = 'model_SLA.uf')
          if (.not.ex) then
             if (master) print *,'model_SLA.uf does not exist'
             call stop_mpi()
          end if
       end if
       if (any(obs(:) % id == 'SSH  ')) then
          inquire(exist = ex, file = 'model_SSH.uf')
          if (.not.ex) then
             if (master) print *,'model_SSH.uf does not exist'
!             call stop_mpi()
          end if
       end if
    end if

    ! construct S=HA
    !
    do iuobs = 1, nuobs

       if (master) then
          print *, 'prep_4_EnKF: now preparing "', trim(unique_obs(iuobs)), '" observations'
       end if

       if (trim(unique_obs(iuobs)) == 'ICEC') then
          do iens = 1, nrens
             write(cmem,'(i3.3)') iens
             call get_micom_fld_new(trim('forecast'//cmem), readfld,&
                  'icec', 0, 1, nx, ny)
             call Generate_element_Si(S(:, iens), unique_obs(iuobs),&
                  readfld, depths, nx, ny, nz, 0) 
          end do

       elseif (trim(unique_obs(iuobs)) == 'SST') then
          do iens = 1, nrens
             write(cmem,'(i3.3)') iens
             call get_micom_fld_new(trim('forecast'//cmem), readfld, &
                  'temp', 1, 1, nx, ny)
#ifdef ANOMALY
!Read the monthly mean 
                call get_micom_fld_new('mean_mod', field2, &
                  'sst', 0, 1, nx, ny)
                readfld = readfld - field2
#endif
             call Generate_element_Si(S(:, iens), unique_obs(iuobs),&
                  readfld, depths, nx, ny, nz, 0) 
          end do

       elseif (trim(unique_obs(iuobs)) == 'SSH') then
          do iens = 1, nrens
             write(cmem,'(i3.3)') iens
             call get_micom_fld_new(trim('forecast'//cmem), readfld, &
                  'sealv', 1, 1, nx, ny)
             ! Convert to meter
             readfld = readfld/100.
             
             !Read the monthly mean
             call get_micom_fld_new('mean_mod', field2, &
                  'sealv', 0, 1, nx, ny)
             readfld = readfld - field2
             
             call Generate_element_Si(S(:, iens), unique_obs(iuobs),&
                  readfld, depths, nx, ny, nz, 0)
          end do

       elseif (trim(unique_obs(iuobs)) == 'SLA' .or.&
            trim(unique_obs(iuobs)) == 'TSLA') then

          if (trim(unique_obs(iuobs)) == 'TSLA') then
             call read_mean_ssh(field2, nx, ny)
          end if
          
          do iens = 1, nrens
                write(cmem,'(i3.3)') iens
             call get_micom_fld_new(trim('forecast'//cmem), readfld, &
              'ssh', 1, 1, nx, ny)
             if (prm_prmestexists('msshb')) then
                write(cmem,'(i3.3)') iens
                call get_micom_fld_new(trim('forecast'//cmem), field3,&
                     'msshb', 0, 1, nx, ny)
                readfld = readfld - field3 ! mean SSH bias for this member
             end if
             if (trim(unique_obs(iuobs)) == 'TSLA') then
                readfld = readfld - field2 ! mean SSH
             end if
             
             call Generate_element_Si(S(:, iens), unique_obs(iuobs),&
                  readfld, depths, nx, ny, nz, t)
          end do
          if (master) then
             print *, 'forming S, day', t
             print *, '  # of non-zero ens observations = ', count(S /= 0.0)
             print *, '  # of zero ens observations = ', count(S == 0.0)
             print *, '  # of non-zero observations for member 1 = ', count(S(:, 1) /= 0.0)
          end if

       elseif (trim(unique_obs(iuobs)) == 'SAL' .or.&
            trim(unique_obs(iuobs)) == 'TEM' .or.&
            trim(unique_obs(iuobs)) == 'GSAL' .or.&
            trim(unique_obs(iuobs)) == 'GTEM') then

          if (master) then
             print *, '  Interpolating ensemble vectors to the locations of "',&
                  trim(unique_obs(iuobs)), '" observations'
          end if

#ifdef ANOMALY
          !Read the monthly mean
          call get_climato_spline(C, trim(unique_obs(iuobs)), nobs, obs)
#endif
          !
          ! for each ensemble member process all profiles "in parallel",
          ! reading the fields layer by layer
          !
          do iens = 1, nrens
             call get_S_spline(S(:, iens), trim(unique_obs(iuobs)), nobs, obs, iens)
#ifdef ANOMALY
             S(uobs_begin(iuobs) : uobs_end(iuobs), iens) = S(uobs_begin(iuobs) : uobs_end(iuobs), iens) &
                  - C(uobs_begin(iuobs) : uobs_end(iuobs))
#endif
          end do
          if (master) then
             print *, '  Interpolation completed'
          end if
          
       else
          if (master) then 
             print *,'ERROR: unknown obs type ' // trim(unique_obs(iuobs))
          end if
          call stop_mpi()
       end if
    end do ! iuobs

    ! some generic QC - relax fitting if the model and obs are too far apart
    !
    call obs_QC(nrens, S)

    ! add calculated HA to to observations-<type>.nc files for each data type
    !
    do iuobs = 1, nuobs
       if (master) then
          nthisobs = 0
          do m = 1, nobs
             if (trim(unique_obs(iuobs)) == trim(obs(m) % id)) then
                nthisobs = nthisobs + 1
                thisobs(nthisobs) = m
             end if
          end do

          ! add forecast values to the observation-<TYPE>.nc files
          !
          call add_forecast(unique_obs(iuobs), S(thisobs(1 : nthisobs), :), obs(thisobs(1 : nthisobs)))

          ! append the superobed values (and modified observation error
          ! variances) to the file with pre-processed observations (SAL.nc,
          ! TEM.nc, GSAL.nc or GTEM.nc)
          !
          if (trim(unique_obs(iuobs)) == 'SAL' .or.&
               trim(unique_obs(iuobs)) == 'TEM' .or.&
               trim(unique_obs(iuobs)) == 'GSAL' .or.&
               trim(unique_obs(iuobs)) == 'GTEM') then
          
             call insitu_writeforecast(unique_obs(iuobs), nobs, nrens, S, obs)
          end if
       end if
    end do

    if (master) then
       print *, 'm_prep_4_EnKF: end calculating S = HA'
    end if

    x = sum(S, DIM = 2) / real(nrens)
    if (master) print*,'m_prep_4_EnKF: end calculating Hx'
    if (master) then
       print *, 'Hx range = ', minval(x), '-', maxval(x)
       print *, 'mean(Hx) = ', sum(x) / real(nobs)
    end if
    if (master) then
       print *, 'observation range = ', minval(obs % d), '-', maxval(obs % d)
       print *, 'mean(observation) = ', sum(obs % d) / nobs
    end if
    ! Compute HA = HE - mean(HE)
    !
    if (master) print*,'prep_4_EnKF(): calculating S = S - x'
    do j = 1, nrens
       S(:, j) = S(:, j) - x
    enddo
    ! Compute innovation
    !
    d = obs % d - x
    if (master) then
       print *, '  innovation range = ', minval(d), '-', maxval(d)
       if (minval(d) < -1000.0d0) then
          print *, 'm_prep_4_EnKF: error: innovation too small detected'
          call stop_mpi()
       end if
       if (maxval(d) > 1000.0d0) then
          print *, 'm_prep_4_EnKF: error: innovation too big detected'
          call stop_mpi()
       end if
    end if

  end subroutine prep_4_EnKF


  subroutine read_mean_ssh(mean_ssh, nx, ny)
#if defined (QMPI)
    use qmpi
#else
    use qmpi_fake
#endif
    use m_parameters

    integer, intent(in) :: nx, ny
    real, intent(out):: mean_ssh(nx, ny)
    logical :: exists

    inquire(file = trim(MEANSSHFNAME), exist = exists)
    if (.not. exists) then
       if (master) then
          print *,'ERROR: read_mean_ssh(): file "', trim(MEANSSHFNAME), '" not found'
       end if
       stop
    end if
       
    open (10, file = trim(MEANSSHFNAME), status = 'unknown',form = 'unformatted')
    read (10) mean_ssh
    close (10)
  end subroutine read_mean_ssh


  ! This subroutine adds forecast observations (i.e Hx) to the NetCDF
  ! observation files for each data type.
  !
  subroutine add_forecast(obstag, S, obs)
    use mod_measurement
    use nfw_mod
    implicit none
    
    character(OBSTYPESTRLEN), intent(in) :: obstag
    real, dimension(:, :), intent(in) :: S
    type(measurement), dimension(:) :: obs

    character(STRLEN) :: fname
    logical :: exists
    integer :: ncid
    integer :: dids(2), dimlen
    logical :: addsobs
    integer :: for_id, inn_id, forvar_id, slon_id, slat_id,&
         sdepth_id, sipiv_id, sjpiv_id, sd_id, svar_id,&
         newvar_id
    
    real, allocatable, dimension(:) :: x, Svar, innovation
  
    integer :: m, p, o

    write(fname, '(a, a, a)') 'observations-', trim(obstag), '.nc'
    inquire(file = trim(fname), exist = exists)
    if (.not. exists) then
       print *, 'file "', trim(fname), 'not found, skip adding forecast'
       return
    else
       print *, 'dumping forecast to "', trim(fname), '"'
    end if

    p = size(S, DIM = 1)
    m = size(S, DIM = 2)

    allocate(x(p), Svar(p), innovation(p))

    x = sum(S, DIM = 2) / real(m);
    Svar = 0.0
    do o = 1, p
       Svar(o) = sum((S(o, :) - x(o))** 2)
    end do
    Svar = Svar / real(m - 1)
    innovation = obs % d - x
  
    addsobs = .false.
    call nfw_open(fname, nf_write, ncid)
    call nfw_inq_dimid(fname, ncid, 'nobs', dids(1))
    call nfw_inq_dimlen(fname, ncid, dids(1), dimlen)

    call nfw_redef(fname, ncid)
    if (dimlen == p) then
       dids(2) = dids(1)
    elseif (.not. nfw_dim_exists(ncid, 'nsobs')) then
       addsobs = .true.
       call nfw_def_dim(fname, ncid, 'nsobs', p, dids(2))
       call nfw_def_var(fname, ncid, 'slon', nf_float, 1, dids(2), slon_id)
       call nfw_def_var(fname, ncid, 'slat', nf_float, 1, dids(2), slat_id)
       call nfw_def_var(fname, ncid, 'sdepth', nf_float, 1, dids(2), sdepth_id)
       call nfw_def_var(fname, ncid, 'sipiv', nf_int, 1, dids(2), sipiv_id)
       call nfw_def_var(fname, ncid, 'sjpiv', nf_int, 1, dids(2), sjpiv_id)
       call nfw_def_var(fname, ncid, 'sd', nf_float, 1, dids(2), sd_id)
       call nfw_def_var(fname, ncid, 'svar', nf_float, 1, dids(2), svar_id)
    end if
    if (.not. nfw_var_exists(ncid, 'innovation')) then
       call nfw_def_var(fname, ncid, 'innovation', nf_double, 1, dids(2), inn_id)
    else
       call nfw_inq_varid(fname, ncid, 'innovation', inn_id)
    end if
    if (.not. nfw_var_exists(ncid, 'forecast')) then
       call nfw_def_var(fname, ncid, 'forecast', nf_double, 1, dids(2), for_id)
    else
       call nfw_inq_varid(fname, ncid, 'forecast', for_id)
    end if
    if (.not. nfw_var_exists(ncid, 'forecast_variance')) then
       call nfw_def_var(fname, ncid, 'forecast_variance', nf_double, 1, dids(2), forvar_id)
    else
       call nfw_inq_varid(fname, ncid, 'forecast_variance', forvar_id)
    end if
    if (.not. addsobs) then
       if (dimlen == p) then
          if (.not. nfw_var_exists(ncid, 'new_var')) then
             call nfw_def_var(fname, ncid, 'new_var', nf_double, 1, dids(2), newvar_id)
          else
             call nfw_inq_varid(fname, ncid, 'new_var', newvar_id)
          end if
       else
          if (.not. nfw_var_exists(ncid, 'new_svar')) then
             call nfw_inq_dimid(fname, ncid, 'nsobs', dids(2))
             call nfw_def_var(fname, ncid, 'new_svar', nf_double, 1, dids(2), newvar_id)
          else
             call nfw_inq_varid(fname, ncid, 'new_svar', newvar_id)
          end if
       end if
    end if
    call nfw_enddef(fname, ncid)

    call nfw_put_var_double(fname, ncid, forvar_id, Svar)
    call nfw_put_var_double(fname, ncid, for_id, x)
    call nfw_put_var_double(fname, ncid, inn_id, innovation)
    if (addsobs) then
       call nfw_put_var_double(fname, ncid, slon_id, obs % lon)
       call nfw_put_var_double(fname, ncid, slat_id, obs % lat)
       call nfw_put_var_double(fname, ncid, sdepth_id, obs % depth)
       call nfw_put_var_int(fname, ncid, sipiv_id, obs % ipiv)
       call nfw_put_var_int(fname, ncid, sjpiv_id, obs % jpiv)
       call nfw_put_var_double(fname, ncid, sd_id, obs % d)
       call nfw_put_var_double(fname, ncid, svar_id, obs % var)
    else
       call nfw_put_var_double(fname, ncid, newvar_id, obs % var)
    end if

    call nfw_close(fname, ncid)

    deallocate(x)
    deallocate(Svar)
    deallocate(innovation)
  end subroutine add_forecast

end module m_prep_4_EnKF
