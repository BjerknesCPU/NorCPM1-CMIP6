c --- ------------------------------------------------------------------
c --- common blocks related to the computation of air-sea fluxes
c --- ------------------------------------------------------------------
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
c
c --- computed transfer coeffisients, gustiness squared, and air density
c --- for dataset
     .  cd_d,ch_d,ce_d,wg2_d,rhoa,
c
c --- computed transfer coeffisients and gustiness squared for model
     .  cd_m,ch_m,ce_m,wg2_m
c
      common /asf/ cd_d,ch_d,ce_d,wg2_d,rhoa,cd_m,ch_m,ce_m,wg2_m
c
c --- constants set in 'asfdat'
      real zu,zt,zq,emiss,cpair,stefanb
      integer tciter
      common /asfpar/ zu,zt,zq,emiss,cpair,stefanb,tciter
