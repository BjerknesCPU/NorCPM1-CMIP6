module genaero_intr

  !---------------------------------------------------------------------------------
  ! Module to interface the aerosol parameterizations with CAM
  ! written by PJR (extensively modified by Ø. Seland for use in CAm-Oslo) 
  !---------------------------------------------------------------------------------
! OS Based on and to a large extent copied from dust_intr

  use shr_kind_mod,only: r8 => shr_kind_r8
  use spmd_utils,  only: masterproc
  use ppgrid,      only: pcols, pver,pverp
  use physconst,   only: mwdry, mwh2o,gravit,rair
  use constituents,only: pcnst, cnst_add, cnst_name, cnst_get_ind
  use rgrid,       only: fullgrid
  use abortutils,  only: endrun
  use cam_logfile, only: iulog
  use aerosoldef,     only:ixac,ncui,ixae
  implicit none

  private          ! Make default type private to the module

  save


  integer :: ncyear


  !
  ! Public interfaces
  !
  public genaero_wet_intr                             ! interface to wet deposition
  public genaero_emis_intr
  
  public genaero_drydep_intr
contains

  subroutine genaero_wet_intr (state, ptend, nstep, dt, lat, clat, cme, &
prain, evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr, &
lo3,lh2o2,wetdepflx) 
    !----------------------------------------------------------------------- 
    ! 
    ! Author O. Seland
    ! Copied from general wet deposition routine described below
    !
    ! Purpose: 
    ! Interface to wet processing of aerosols (source and sinks).
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: B.A. Boville
    ! 
    !-----------------------------------------------------------------------
    use cam_history,   only: outfld
    use physics_types, only: physics_update,physics_state, physics_ptend,physics_tend
    use aerosoldef
    use mass,         only: cmidry,cmidpdivg     

    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments:
    !
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    integer, intent(in) :: nstep
    integer, intent(in) :: lat(pcols)                  ! latitude index for S->N storage
    real(r8), intent(in) :: clat(pcols)                    ! latitude 
    real(r8), intent(in) :: cme(pcols,pver)            ! local condensation of cloud water
    real(r8), intent(in) :: prain(pcols,pver)            ! production of rain
    real(r8), intent(in) :: evapr(pcols,pver)            ! evaporation of rain
    real(r8), intent(in) :: cldn(pcols,pver)            ! cloud fraction
    real(r8), intent(in) :: cldc(pcols,pver)            ! convective cloud fraction
    real(r8), intent(in) :: cldv(pcols,pver)            ! cloudy volume undergoing scavenging
    real(r8), intent(in) :: conicw(pcols, pver)
    real(r8), intent(in) :: cmfdqr(pcols, pver)
    real(r8), intent(in) :: rainmr(pcols, pver) ! rain mixing ratio

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies

    real(r8), intent(inout) :: fracis(pcols,pver,pcnst)         ! fraction of transported species that are insoluble
  real(r8), intent(in) :: lo3(pcols,pver) ! O3 concentration
  real(r8), intent(in) :: lh2o2(pcols,pver) ! H2O2 concentration
  real(r8), intent(out) :: wetdepflx(pcols,pcnst) ! constituent wetdep fluxes (kg/m2/s)
!
! Local variables
!
   type(physics_tend ) :: tend  
! Physics tendencies (empty, needed for physics_update call)
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
    integer :: ix
    real(r8) :: sflx(pcols)            ! deposition flux

    real(r8) :: obuf(1)
    real(r8), intent(in) :: calday        ! current calendar day
    real(r8) :: iscavt(pcols, pver)
   real(r8) :: scavt(pcols,pver,pcnst)
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    integer :: mm,i,k
    integer :: m                                  ! tracer index
    integer :: ixcldliq
    integer :: ixcldice
    real(r8) totcond(pcols, pver) ! total condensate
    real(r8) :: sol_fact
    real(r8) :: fice(pcols,pver)     ! Fraction of cwat that is ice. Taken from Iversen and Seland 2003
    real(r8) :: cmi2d(pcols)


   real(r8) :: aqprod(pcols)     ! aqueous production of sulphate(kg/m2/s)
   real(r8) :: qstate(pcols,pver,pcnst)

    !-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    sflx(:)=0._r8

    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    totcond(:ncol,:) = state%q(:ncol,:,ixcldliq) + state%q(:ncol,:,ixcldice)


  do k=1,pver
     do i=1,ncol
! If warmer than -10 degrees C then water phase
!
         if (state%t(i,k) > 273.16_r8) fice(i,k) = 0.0_r8
!
! If colder than 0 degrees C but warmer than -25 C mixed phase
!
         if (state%t(i,k) <= 273.16_r8 .and. state%t(i,k) >= 248.16_r8) then
            fice(i,k) =(273.16_r8-state%t(i,k)) / 25.0_r8
         end if
!
! If colder than -25 degrees C then ice phase
!
         if (state%t(i,k) < 248.16_r8) fice(i,k) = 1.0_r8
      end do
   end do	
      qstate(:,:,:)=state%q(:,:,:)
      ptend%lq(ixac:pcnst)=.true.   

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Toggle the chemistry and wet deposition routines !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if( mod( nstep/2, 2 ) .eq. 0 ) then

!   ptend%name  = 'aqgas'
!   call t_startf ('aqgas')
      call aqgas(lchnk,   ncol,    dt,              &
                   qstate, state%t, state%pmid,    cldv,        &
                   state%pdel,   conicw,       &
                   cmfdqr,  prain,   cme,   evapr,             &
                   wetdepflx,                  &
                   aqprod,  fracis,  &
                    lo3,lh2o2,totcond,ixcldliq) 
!   call t_stopf ('aqgas')

!   call physics_update (state, tend, ptend, dt)

   scavt(:,:,:)=0._r8

   do m=ixac,pcnst
     if(species_class(m)==spec_class_aerosol) then

            call wetdepaer (lchnk,ncol,dt, state%pdel, cldv, &
          cmfdqr, prain, cme,                                &
          evapr, fice, state%q(:,:,ixcldliq), qstate(:,:,m),  &
          conicw, scavconv(m),                                 &
          scavt(:,:,m), fracis(:,:,m),scavin(m),scavbe(m),scavbecon(m))   

          do k=1,pver
             do i=1,ncol
                qstate(i,k,m)=qstate(i,k,m)+scavt(i,k,m)*dt
!                ptend%q(i,k,m)=scavt(i,k,m)
             end do
          end do     
       call cmidry( lchnk,ncol,state%pdel, state%q(:,:,1), scavt(:,:,m), cmi2d)
       do i=1,ncol
         wetdepflx(i,m)=cmi2d(i)
       end do          
     end if       
   end do
!   call physics_update (state, tend, ptend, dt)


               !###############################
      else     !# toggle wetdep and chemistry #
               !###############################

   scavt(:,:,:)=0._r8

   do m=ixac,pcnst
     if(species_class(m)==spec_class_aerosol) then
            call wetdepaer (lchnk,ncol,dt, state%pdel, cldv, &
          cmfdqr, prain, cme,                                &
          evapr, fice, state%q(:,:,ixcldliq), qstate(:,:,m),   &
          conicw, scavconv(m),                                 &
          scavt(:,:,m), fracis(:,:,m),scavin(m),scavbe(m),scavbecon(m))   

  
          do k=1,pver
             do i=1,ncol
                qstate(i,k,m)=qstate(i,k,m)+scavt(i,k,m)*dt
!                ptend%q(i,k,m)=scavt(i,k,m)
             end do
          end do     

       call cmidpdivg( lchnk,ncol,state%pdel, scavt(:,:,m), cmi2d)
!       call cmidry( lchnk,ncol,state%pdel, state%q(:,:,1), scavt(:,:,m), cmi2d)
       do i=1,ncol
         wetdepflx(i,m)=cmi2d(i)
       end do          
     end if       
   end do
!   call physics_update (state, tend, ptend, dt)



!   ptend%name  = 'aqgas'
!   call t_startf ('aqgas')
      call aqgas(lchnk,   ncol,    dt,              &
                    qstate ,state%t, state%pmid,    cldv,        &
                   state%pdel,   conicw,       &
                   cmfdqr,  prain,   cme,   evapr,             &
                   wetdepflx,                   &
                   aqprod,  fracis,     &
                    lo3,lh2o2,totcond,ixcldliq) 
                 
!   call t_stopf ('aqgas')

!   call physics_update (state, tend, ptend, dt)


   end if


     ptend%name = 'wetdepchem'
!      ptend%q(:ncol,:,ixac:pcnst)=0._r8
     ptend%q(:ncol,:,ixac:pcnst)=(qstate(:ncol,:,ixac:pcnst)- &
        state%q(:ncol,:,ixac:pcnst))/dt

     call outfld('S4AQ ',aqprod,pcols   ,lchnk )



    return

  end subroutine genaero_wet_intr

  subroutine genaero_emis_intr (state, ptend, dt)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to emission of all carbons
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
    use cam_history,    only: outfld
    use physics_types,  only: physics_state, physics_ptend
    use phys_grid,      only: get_lon_all_p, get_lat_all_p, get_rlat_all_p
    use time_manager,   only: get_curr_date, get_perp_date, get_curr_calday, &
                              is_perpetual
    use error_messages, only: alloc_err, handle_err

!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies

    integer  latidx(pcols)                  ! latitude index 
    integer  lonidx(pcols)                  ! longitude index
    integer lchnk
    integer ncol
    integer i
    integer m, mm, ix
    integer istat
!
    real(r8) :: calday        ! current calendar day
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    real(r8), allocatable :: cflx(:,:)  !surface flux , (pcols,nphob)

    calday = get_curr_calday()
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if

    lchnk = state%lchnk
    ncol = state%ncol
    
    call get_lat_all_p(lchnk, ncol, latidx)
    call get_lon_all_p(lchnk, ncol, lonidx)


!    allocate( cflx(pcols,nphob), stat=istat )
!    call alloc_err( istat, 'carbon_emis_intr', 'cflx', pcols*nphob )
!    call caersf (ncol, nphob, latidx, lonidx, cflx) 


!    ix = caer_idx1()
    ptend%name  = ptend%name//'+emis'
    do m = 1,ncui
       mm = ixac + m - 1
!!       ptend%lq(mm) = .true. ! tendencies for all carbon on
       ptend%q(:ncol,pver,mm) = 0._r8 ! zero all carbon tends
!       ptend%q(:ncol,pver,mm) = cflx(:ncol,m)*gravit/state%pdel(:ncol,plev) ! source
!       call outfld(trim(cnst_name(mm))//'SF', cflx(:,m), pcols, lchnk)
    end do ! m

!    deallocate( cflx, stat=istat )
!    call handle_err( istat, &
!         'ERROR deallocating memory for cflx in routine carbon_emis_intr')


    return
  end subroutine genaero_emis_intr 



  subroutine genaero_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
       fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, month, landfrac, &
       icefrac, ocnfrac,fvin,ram1in,cam_out)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to parameterized greenhouse gas chemisty (source/sink).
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: P. J. Rasch
! 
!-----------------------------------------------------------------------
    use cam_history,       only: outfld
    use physics_types,     only: physics_state, physics_ptend
    use phys_grid,         only: get_lat_all_p,get_lon_all_p
    use constituents,      only: cnst_name
    use drydep_mod,        only: setdvel, ddflux,calcram,d3ddflux
    use dust_sediment_mod, only: dust_sediment_tend
    use physconst,         only: rair
    use ppgrid,            only: pverp
    use aerosoldef
    use camsrfexch_types,  only: cam_out_t
    use time_manager,      only : get_ref_date
    use aero_to_srf, only: set_srf_drydep
    use mo_mass_xforms,    only : h2o_to_vmr
    use shr_const_mod,only : SHR_CONST_MWDAIR
    use time_manager,      only : get_curr_calday
    
    use wv_saturation,     only : aqsat
    use drydep_so2,        only : dvel_so2,set_soilwso2

!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in) :: state          ! Physics state variables

    integer, intent(in) :: lat(pcols)                  ! latitude index for S->N storage
    real(r8), intent(in) :: clat(pcols)                 ! latitude 
    real(r8), intent(in) :: fsds(pcols)                 ! longwave down at sfc
    real(r8), intent(in) :: obklen(pcols)                 ! longwave down at sfc
    real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel
    real(r8), intent(in) :: ts(pcols)                     ! sfc temp
    real(r8), intent(in) :: landfrac(pcols)               ! land fraction
    real(r8), intent(in) :: icefrac(pcols)                ! ice fraction
    real(r8), intent(in) :: ocnfrac(pcols)                ! ocean fraction
    real(r8), intent(in) :: hflx(pcols)                  ! sensible heat flux
    real(r8), intent(in) :: prect(pcols)                     ! prect
    real(r8), intent(in) :: snowh(pcols)                     ! snow depth
    real(r8), intent(in) :: pblh(pcols)                     ! pbl height
    integer, intent(in)  :: month
    real(r8), intent(in) :: wvflx(pcols)       ! water vapor flux
    real(r8), intent(in) :: fvin(pcols)        ! for dry dep velocities from land model
    real(r8), intent(in) :: ram1in(pcols)       ! for dry dep velocities from land model
    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state
!
! Local variables
!
    integer :: i,k,m                                  ! tracer index
    integer :: mm                                  ! tracer index
    integer :: ioff                               ! offset for ghg indices
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
    integer :: ix
    real(r8) :: tvs(pcols,pver)
    real(r8) :: tvsgas(pcols,pver)
    real(r8) :: obuf(1)
    real(r8)  sflx(pcols,pcnst)            ! deposition flux
    real(r8) :: sflxgas(pcols)            ! deposition flux
    real(r8), parameter :: mil  = .001_r8
    logical, parameter :: dyn_soilw = .false.
    real(r8) :: mill,milo,mili  ! Vd land/ocean/ice
    real(r8) :: calday
    real(r8) wet_radius(pcols,pver,ncaer)   ! Effective radius including water
    real(r8) wet_rho(pcols,pver,ncaer)   ! Aerosol density including water
    real(r8), parameter :: vlc_sed_carb  = 0.5e-3_r8 !sedimentation velocity carb cm/s
    real(r8):: pvcarb(pcols,pverp)   !sedimentation velocity by column/level
    real(r8):: rho(pcols,pver)      !density

    real(r8)::  dep_trb(pcols)       !kg/m2/s
    real(r8)::  dep_dry(pcols)       !kg/m2/s (total of grav and trb)
    real(r8)::  dep_grv(pcols)       !kg/m2/s (total of grav and trb)
    real(r8) :: pvmzaer(pcols,pverp)    ! sedimentation velocity in Pa
    real(r8) :: fv(pcols)         ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: ram1(pcols)       ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: dqdt_tmp(pcols,pver)   ! temporary array to hold tendency for 1 species
    real(r8) :: vlc_dry(pcols,pver)       ! dep velocity
    real(r8) :: vlc_grv(pcols,pver)       ! dep velocity
    real(r8)::  vlc_trb(pcols)            ! dep velocity
    real(r8) :: sigpart(pcols,pver)
    integer  :: ncdate,yr,mon,day,sec    
    real(r8) :: wind_speed(pcols)        ! surface wind speed (m/s)

    real(r8)     ::  depvel(pcols)       ! dry deposition velocity (cm/s) for SO2
    real(r8)     ::  dvel(pcols)                ! dry deposition velocity (cm/s)
   real(r8) :: dryso2(pcols),dryso4(pcols),drybc(pcols)
   real(r8) :: drypom(pcols),drysalt(pcols),drydust(pcols)
   integer  :: latndx(pcols)                     ! chunk lat indices
   integer  :: lonndx(pcols)                     ! chunk lon indices

    real(r8) :: relhum(pcols,pver)                          ! relative humidity
    real(r8) :: h2ovmr(pcols,pver)                          ! water vapor volume
    logical  :: table_soilw
    real(r8) :: soilw(pcols)
    real(r8) :: satv(pcols,pver)                            ! wrk array for relative humidity
    real(r8) :: satq(pcols,pver)                            ! wrk array for relative humidity
    real(r8) :: mbar(pcols,pver)
    

!-----------------------------------------------------------------------
!    ix = caer_idx1()
!    ioff  = ix - 1
    lchnk = state%lchnk
    ncol  = state%ncol

    calday = get_curr_calday( )

    call get_lat_all_p(lchnk, ncol, latndx)
    call get_lon_all_p(lchnk, ncol, lonndx)
     call calcaersize(ncol, &
                 state%t, state%q(1,1,1), state%pmid,   &
                 state%pdel, wet_radius,wet_rho,relhum)

    call calcram(ncol,landfrac,icefrac,ocnfrac,obklen,&
         ustar,ram1in,ram1,state%t(:,pver),state%pmid(:,pver),&
         state%pdel(:,pver),fvin,fv)

!    tvs(:ncol) = state%t(:ncol,pver)!*(1+state%q(:ncol,pver)
!    Added 3D fluxes
    tvs(:ncol,:) = state%t(:ncol,:)!*(1+state%q(:ncol,k)
    tvsgas(:ncol,:) = state%t(:ncol,:)*(1+state%q(:ncol,:,1))

!   write(iulog,*) ' carbon drydep invoked '

    do m=1,ncaer
      mm=m+ixae-1
      ptend%lq(mm)=.true.
      sigpart(:,:)=sgpart(mm)
      call aero_depvel_part( ncol, state%t(:,:), state%pmid(:,:), ram1, fv,  & 
            vlc_dry(:,:), vlc_trb(:), vlc_grv(:,:),  &
            wet_radius(:,:,m), wet_rho(:,:,m), sigpart(:,:), 3, lchnk)
!      Write(6,*) vlc_dry(:,pver-1)
    
      call d3ddflux( ncol, vlc_dry(:,:), state%q(:,:,mm), state%pmid, &
                               state%pdel, tvs, sflx(:,mm), ptend%q(:,:,mm), dt )
#ifdef SHORTRUN
       call outfld( trim(cnst_name(mm))//'DRY', sflx(:,mm), pcols, lchnk)
       call outfld( trim(cnst_name(mm))//'DV', vlc_dry(:,pver), pcols, lchnk)
#endif
!       call outfld( trim(cnst_name(mm))//'DV', vlc_grv(:,pver), pcols, lchnk)
      end do
         mbar(:,:) = SHR_CONST_MWDAIR! dry mass

       call get_ref_date(yr, mon, day, sec)
       ncdate = yr*10000 + mon*100 + day
       wind_speed(:ncol) = sqrt( state%u(:ncol,pver)*state%u(:ncol,pver) + state%v(:ncol,pver)*state%v(:ncol,pver) )

    call aqsat( state%t, state%pmid, satv, satq, pcols, &
         ncol, pver, 1, pver )


!       if( .not. dyn_soilw .and. table_soilw ) then
!          call set_soilwso2( soilw, ncol , lonndx, latndx, calday )
!       end if
        soilw(:)=0._r8

!       call drydep( ncdate, state2d%ts, state2d%ps,  &
!            wind_speed, state%q(:,pver), state%t(:,pver), state%pmid(:,pver), &
!            prect, snowh, fsds, depvel, sflxgas, mmr, &
!            tvsgas, soilw, relhum(:,pver:pver), ncol, lonndx, latndx, lchnk )


       call dvel_so2( ncdate, ts, state%ps,  &
            wind_speed, state%q(:,1,pver), state%t(:,pver), state%pmid(:,pver), &
             prect,snowh, fsds, depvel, &
            tvsgas(:,pver), soilw, relhum(:,pver), ncol, lonndx, latndx, lchnk ) 

!    call setdvel( ncol, landfrac, icefrac, ocnfrac, mill, milo, mili, dvel )

        do m = 1, ncgas
        mm = ixac + m - 1

        if (mm.eq.l_so2) then
          ptend%lq(mm) =.TRUE.
          do i=1,ncol
             dvel(i)=depvel(i)
          end do
        else
          dvel(:)=0._r8
        end if  
#ifdef SHORTRUN
       call outfld( trim(cnst_name(mm))//'DRY', sflx(:,mm), pcols, lchnk)

       call outfld( trim(cnst_name(mm))//'DV', dvel(:), pcols, lchnk)
#endif
!     if(species_class(m)==spec_class_aerosol) then 
       call ddflux( ncol, dvel, state%q(:,pver,mm), state%pmid(:,pver), & 
     tvs(:,pver), sflx(:,mm ))
!#endif
       ptend%q(:ncol,pver,mm) = sflx(:ncol,mm)*gravit*state%rpdel(:ncol,pver)
!#ifdef MATCH
!       call outfld( trim(cnst_name(mm))//'DRY', sflx, pcols, lat, obuf )
!#else


!#endif


    end do


! 


! set flags for tracer tendencies 
    ptend%name  = ptend%name//'+drydep'

! Combine total dry-deposition species into SO4,BC,POM,DUST and SS

   do i=1,ncol  

!     drymsa(i) = sflx(i,l_msa)
     dryso2(i) = sflx(i,l_so2)
     dryso4(i) = sflx(i,l_so4_n)+sflx(i,l_so4_na)   &
+ sflx(i,l_so4_a1)+sflx(i,l_so4_a2)+sflx(i,l_so4_ac)+sflx(i,l_so4_pr)
     drybc(i) = sflx(i,l_bc_n)+sflx(i,l_bc_ax) &
+ sflx(i,l_bc_ni)+sflx(i,l_bc_a)+sflx(i,l_bc_ai)+sflx(i,l_bc_ac)
     drypom(i) = sflx(i,l_om_ni)&
+sflx(i,l_om_ai)+sflx(i,l_om_ac)
     drydust(i)= sflx(i,l_dst_a2)+sflx(i,l_dst_a3)
     drysalt(i)= sflx(i,l_ss_a1)+sflx(i,l_ss_a2)+sflx(i,l_ss_a3)

   end do
    call set_srf_drydep(sflx, cam_out)
!      call outfld('DRY_MSA ',drymsa   ,pcols   ,lchnk   )
      call outfld('DRY_SO2 ',dryso2   ,pcols   ,lchnk   )
      call outfld('DRY_SO4 ',dryso4   ,pcols   ,lchnk   )
      call outfld('DRY_BC ',drybc   ,pcols   ,lchnk   )
      call outfld('DRY_POM ',drypom   ,pcols   ,lchnk   )
      call outfld('DRY_DUST ',drydust   ,pcols   ,lchnk   )
      call outfld('DRY_SS ',drysalt   ,pcols   ,lchnk   )


!
! record tendencies on history files
!    do m = 1, 4
!       call outfld ('DRYDEP'//srcnam(m),ptend%q(:,:,ioff+m),pcols,lchnk)
!    end do

    return
  end subroutine genaero_drydep_intr



end module genaero_intr
