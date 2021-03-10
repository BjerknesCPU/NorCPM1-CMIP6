c --- ------------------------------------------------------------------
c --- common blocks related to the reading of synoptic atmospheric
c --- forcing fields
c --- ------------------------------------------------------------------
c
c --- weights and indexes for distributing runoff to coastal points
      real*4, dimension(atm_abdm,atm_idm,atm_jdm) :: rnf_weight
      integer*2, dimension(atm_abdm,atm_idm,atm_jdm) ::
     .  rnf_ocdpi,rnf_ocdpj
c
c --- path to atmospheric synoptic forcing fields
      character*120 atm_path
      integer atm_path_len
c
      common /atmsyn/ rnf_weight,rnf_ocdpi,rnf_ocdpj,atm_path,
     .                 atm_path_len
