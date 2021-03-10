c --- ------------------------------------------------------------------
c --- common blocks related to the application of param pert
c --- ------------------------------------------------------------------
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  rrm0,          ! wind stress
     .  rbdmc2,        ! Vertical diffusion
     .  regc           ! c in Eden and Greatbatch
c
c
      common /assim/ rrm0,rbdmc2,regc
