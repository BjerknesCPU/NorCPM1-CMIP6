c --- ------------------------------------------------------------------
c --- common blocks related to the application of monthly climatological
c --- forcing fields
c --- ------------------------------------------------------------------
c
c --- monthly forcing fields
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12) ::
     .  taud,          ! wind stress
     .  tauxd,tauyd,   ! wind stress components
     .  dswrfl,        ! short-wave surface radiation
     .  nlwrfs,        ! net long-wave surface radiation
     .  shtflx,        ! sensible heat flux
     .  lhtflx,        ! latent heat flux
     .  precip,        ! precipitation
     .  clouds,        ! cloud cover
     .  slpres,        ! sea level pressure
     .  runoff         ! runoff
c
c --- interpolation parameters monthly climatological fields
      real xx
      integer ll1,ll2,ll3,ll4,ll5
c
      common /clim/ taud,tauxd,tauyd,dswrfl,nlwrfs,shtflx,lhtflx,precip,
     .              clouds,slpres,runoff,xx,ll1,ll2,ll3,ll4,ll5
