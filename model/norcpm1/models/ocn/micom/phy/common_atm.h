c --- ------------------------------------------------------------------
c --- common blocks related to the reading and interpolation of
c --- atmospheric data
c --- ------------------------------------------------------------------
c
      real, dimension(atm_idm,atm_jdm,12) ::
     .  atm_icec,      ! climatological ice cover
     .  atm_sktclm     ! climatological surface temperature
c
      real, dimension(atm_idm,atm_jdm) ::
     .  atm_lon,       ! longitudes
     .  atm_lat,       ! latitudes
     .  atm_topo       ! topography
c
      integer, dimension(atm_idm,atm_jdm) ::
     .  atm_msk        ! mask
c
      real*4, dimension(atm_nwgt,itdm,jtdm) ::
     .  atm_wgt        ! interpolation weights
c
      integer*2, dimension(atm_nwgt,itdm,jtdm) ::
     .  atm_iwgt,      ! interpolation adresses
     .  atm_jwgt       ! interpolation adresses
c
      integer, dimension(itdm,jtdm) ::
     .  itp            ! mask for the full domain
c
      real
     .  atm_ice_swgt,  ! smoothing weight for fields over ice
     .  atm_rnf_swgt   ! smoothing weight for fields over ice
c
      integer
     .  atm_ice_nsmt,  ! number of smoothing iterations over ice
     .  atm_rnf_nsmt   ! number of smoothing iterations over ice
c
      common /atm1/ atm_icec,atm_sktclm,atm_lon,atm_lat,atm_topo,
     .              atm_msk,atm_wgt,atm_iwgt,atm_jwgt,itp,atm_ice_swgt,
     .              atm_rnf_swgt,atm_ice_nsmt,atm_rnf_nsmt
c
c --- constants set in 'atmdat'
      real
     . atm_mval        ! value of a point that should not have a
                       ! physical value
     .,atm_fval        ! value of a point that do not have a physical
                       ! value but should have one
     .,atm_ice_csmt    ! constant determining how much the atm. fields
                       ! are smoothed over ice covered regions
     .,atm_rnf_csmt    ! constant determining how much the runoff is
                       ! smoothed at the coastal discharge points
     .,atm_crnf        ! runoff adjustment factor
     .,atm_cswa        ! short-wave radiation adjustment factor
c
      common /atm2/ atm_mval,atm_fval,atm_ice_csmt,atm_rnf_csmt,
     .              atm_crnf,atm_cswa
