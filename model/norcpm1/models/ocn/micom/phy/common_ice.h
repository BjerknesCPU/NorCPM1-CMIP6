c --- ------------------------------------------------------------------
c --- common blocks for the ice and snow part of the model
c --- ------------------------------------------------------------------
c   
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  ficem,hicem,tsrfm,hsnwm,ticem,ustari,tauxice,tauyice,uicem,
     .  vicem,iagem
c 
      common /icesnw/ ficem,hicem,tsrfm,hsnwm,ticem,ustari,tauxice,
     .                tauyice,uicem,vicem,iagem
c
c --- constants set in 'icedat'
      real
     .  albi_f,        ! --                 max albedo over ice
     .  albi_m,        ! --                 max albedo over melting ice
     .  albs_f,        ! --                 albedo over snow
     .  albs_m,        ! --                 albedo over melting snow
     .  rhoice,        ! kg / m^3           density of ice
     .  rhosnw,        ! kg / m^3           density of snow
     .  rkice,         ! w / (m k)          ice conductivity
     .  fusi,          ! j / kg             heat of fusion of ice
     .  fuss,          ! j / kg             heat of fusion of snow
     .  fice_max,      ! --                 maximum fractional ice cover
     .  tice_m,        ! k                  melting point of ice
     .  tsnw_m,        ! k                  melting point of snow
     .  hice_nhmn,     ! m                  min. ice thickness northern hemi.
     .  hice_shmn,     ! m                  min. ice thickness southern hemi.
     .  gamma,         ! 1 / s              snow aging timescale
     .  sice,          ! per mil            salinity of seaice
     .  sref,          ! per mil            global ref. surface salinity
     .  rksnw,         ! w / (m k)          snow conductivity
     .  cwi,           ! --                 ice-ocean heat transfer coeff.
     .  cuc,           ! w / (m^2 k)        const. for heat flux associated
                       !                    with under-cooled water, resulting
                       !                    in a temp. adjustment of a 20 m
                       !                    mixed layer towards freezing point
                       !                    with an e-folding timescale of
                       !                    approx. one day.
     .  cdiff,         ! m^2 / s            horizontal diffusivity
     .  cdfac          !                    cdiff is multiplied by cdfac in
                       !                    channels of one grid point width
c
      common /iceprm/ albi_f,albi_m,albs_f,albs_m,rhoice,rhosnw,rkice,
     .                fusi,fuss,fice_max,tice_m,tsnw_m,hice_nhmn,
     .                hice_shmn,gamma,sice,sref,rksnw,cwi,cuc,cdiff,
     .                cdfac
