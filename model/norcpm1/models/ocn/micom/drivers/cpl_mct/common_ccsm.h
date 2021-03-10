c --- ------------------------------------------------------------------
c --- common blocks related to the application of CCSM forcing fields
c --- ------------------------------------------------------------------
c
c --- daily average forcing fields
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     .  swa_da,        ! solar heat flux
     .  nsf_da,        ! non-solar heat flux
!+tht
     .  nsf_da_atm,    ! non-solar heat flux
     .  sst_da_atm,
!-tht
     .  hmlt_da,       ! heat flux due to melting
     .  lip_da,        ! liquid water flux
     .  sop_da,        ! solid precipitation
     .  eva_da,        ! evaporation
     .  rnf_da,        ! runoff, liquid
     .  rfi_da,        ! runoff, frozen
     .  fmltfz_da,     ! fresh water flux due to melting/freezing
     .  sfl_da,        ! salt flux
     .  ztx_da,        ! u component of wind stress
     .  mty_da,        ! v component of wind stress
c+sy 2016/10/13
     .  ztx_da_atm,        ! u component of wind stress to atmos
     .  mty_da_atm,        ! v component of wind stress to atmos
c-sy 2016/10/13
     .  ustarw_da,     ! friction velocity for open water
     .  slp_da,        ! sea level pressure
     .  abswnd_da,     ! wind speed at measurement height -zu-
     .  atmco2_da,     ! atmospheric co2 concentration
     .  ficem_da       ! ice concentration
c
c --- time level indices for CCSM fields
      integer ll1,ll2
c
      common /ccsm/ swa_da,nsf_da,
!+tht
     .              nsf_da_atm,sst_da_atm,
!-tht
     .              hmlt_da,lip_da,sop_da,eva_da,rnf_da,
     .              rfi_da,fmltfz_da,sfl_da,ztx_da,mty_da,
c+tht
     .              ztx_da_atm,mty_da_atm,
c-tht
     .              ustarw_da,
     .              slp_da,abswnd_da,atmco2_da,ficem_da,ll1,ll2
