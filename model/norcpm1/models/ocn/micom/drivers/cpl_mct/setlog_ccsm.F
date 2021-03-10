      subroutine setlog_ccsm
c
c --- ------------------------------------------------------------------
c --- Direct stdout to CCSM defined log file
c --- ------------------------------------------------------------------
c
      use mod_xc
      use shr_file_mod, only: shr_file_getUnit, shr_file_setIO
c
      implicit none
c
      logical exists
c
      if (mnproc.eq.1) then
        inquire(file='ocn_modelio.nml',exist=exists)
        if (exists) then
          lp = shr_file_getUnit()
          call shr_file_setIO('ocn_modelio.nml',lp)
        endif
      endif
c
      return
      end
