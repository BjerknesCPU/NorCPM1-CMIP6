      subroutine rdlim(nday1,nday2)
c
c --- ------------------------------------------------------------------
c --- Read limits file
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_dia
#ifdef CCSMCOUPLED
      use data_mct, only: runid_mct, runtype_mct
      use mod_nctools
#endif
#ifdef NUDGE_TS
      use mod_nudge_ts
#endif
c
      implicit none
c
      integer :: nday1,nday2
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_forc.h"
#ifdef SYN
#  include "param_atm.h"
#  include "common_atmsyn.h"
#else 
      character*120 :: atm_path
#endif
c
      integer :: i,j,m,n
#ifdef CCSMCOUPLED
      character*120 :: rstfnm
      logical :: fexist
#endif
      integer :: idate,idate0
      namelist /limits/ nday1,nday2,idate,idate0,runid,baclin,batrop,
     .  mdv2hi,mdv2lo,mdv4hi,mdv4lo,mdc2hi,mdc2lo,vsc2hi,vsc2lo,vsc4hi,
     .  vsc4lo,cbar,cb,cwbdts,cwbdls,egc,eggam,egmndf,egmxdf,egidfq,
     .  trxday,srxday,trxdpt,srxdpt,trxlim,srxlim,aptflx,apsflx,ditflx,
     .  disflx,srxbal,srxsrc,smtfrc,sprfac,path,path1,path2,atm_path,
     .  rstfrq,rstfmt,rstcmp
c
c --- read namelists from limits file
      GLB_AVEPERIO=0
#ifdef CCSMCOUPLED
      open (unit=nfu,file='ocn_in',status='old',action='read',recl=80)
#else
      open (unit=nfu,file='limits',status='old',action='read',recl=80)
#endif
      read (unit=nfu,nml=LIMITS)
#ifdef NUDGE_TS
      read (unit=nfu,nml=NUDGETS)
#endif
      read (unit=nfu,nml=DIAPHY)
      close (unit=nfu)
c
c --- print namelist on stdout 
      if (mnproc.eq.1) then
        write (lp,'(a)') 'MICOM LIMITS-NAMELIST:' 
        write (lp,*) 'NDAY1',NDAY1
        write (lp,*) 'NDAY2',NDAY2
        write (lp,*) 'IDATE',IDATE
        write (lp,*) 'IDATE0',IDATE0
        write (lp,*) 'RUNID ',trim(RUNID)
        write (lp,*) 'BACLIN',BACLIN
        write (lp,*) 'BATROP',BATROP
        write (lp,*) 'MDV2HI',MDV2HI
        write (lp,*) 'MDV2LO',MDV2LO
        write (lp,*) 'MDV4HI',MDV4HI
        write (lp,*) 'MDV4LO',MDV4LO
        write (lp,*) 'MDC2HI',MDC2HI
        write (lp,*) 'MDC2LO',MDC2LO
        write (lp,*) 'VSC2HI',VSC2HI
        write (lp,*) 'VSC2LO',VSC2LO
        write (lp,*) 'VSC4HI',VSC4HI
        write (lp,*) 'VSC4LO',VSC4LO
        write (lp,*) 'CBAR',CBAR
        write (lp,*) 'CB',CB
        write (lp,*) 'CWBDTS',CWBDTS
        write (lp,*) 'CWBDLS',CWBDLS
        write (lp,*) 'EGC',egc
        write (lp,*) 'EGGAM',eggam
        write (lp,*) 'EGMNDF',egmndf
        write (lp,*) 'EGMXDF',egmxdf
        write (lp,*) 'EGIDFQ',egidfq
        write (lp,*) 'TRXDAY',TRXDAY
        write (lp,*) 'SRXDAY',SRXDAY
        write (lp,*) 'TRXDPT',TRXDPT
        write (lp,*) 'SRXDPT',SRXDPT
        write (lp,*) 'TRXLIM',TRXLIM
        write (lp,*) 'SRXLIM',SRXLIM
        write (lp,*) 'APTFLX',APTFLX
        write (lp,*) 'APSFLX',APSFLX
        write (lp,*) 'DITFLX',DITFLX
        write (lp,*) 'DISFLX',DISFLX
        write (lp,*) 'SRXBAL',SRXBAL
        write (lp,*) 'SRXSRC ',trim(SRXSRC)
        write (lp,*) 'SMTFRC',SMTFRC
        write (lp,*) 'SPRFAC',SPRFAC
        write (lp,*) 'PATH ',trim(PATH)
        write (lp,*) 'PATH1 ',trim(PATH1)
        write (lp,*) 'PATH2 ',trim(PATH2)
        write (lp,*) 'ATM_PATH ',trim(ATM_PATH)
        write (lp,*) 'RSTFRQ',RSTFRQ
        write (lp,*) 'RSTFMT',RSTFMT
        write (lp,*) 'RSTCMP',RSTCMP
      endif 
c
c --- determine number of io groups and print namelist
      nphy=0  
      do n=1,nphymax 
        if (GLB_AVEPERIO(n).gt.0) nphy=nphy+1 
      end do     
      if (trxday.eq.0.) then
        SRF_SURRLX(1:nphy)=0
      endif
      if (srxday.eq.0.) then
        SRF_SALRLX(1:nphy)=0
      endif
      if (mnproc.eq.1) then
        write(lp,'(a)') ' ' 
        write(lp,'(a)') 'MICOM IO_NAMELIST ' 
        write(lp,*) 'GLB_FNAMETAG',GLB_FNAMETAG(1:nphy)
        write(lp,*) 'GLB_AVEPERIO',GLB_AVEPERIO(1:nphy)
        write(lp,*) 'GLB_FILEFREQ',GLB_FILEFREQ(1:nphy)
        write(lp,*) 'GLB_COMPFLAG',GLB_COMPFLAG(1:nphy)
        write(lp,*) 'GLB_NCFORMAT',GLB_NCFORMAT(1:nphy)
        write(lp,*) 'SRF_ABSWND  ',SRF_ABSWND(1:nphy)
        write(lp,*) 'SRF_ALB     ',SRF_ALB(1:nphy)
        write(lp,*) 'SRF_BRNFLX  ',SRF_BRNFLX(1:nphy)
        write(lp,*) 'SRF_BRNPD   ',SRF_BRNPD(1:nphy)
        write(lp,*) 'SRF_DFL     ',SRF_DFL(1:nphy)
        write(lp,*) 'SRF_EVA     ',SRF_EVA(1:nphy)
        write(lp,*) 'SRF_FMLTFZ  ',SRF_FMLTFZ(1:nphy)
        write(lp,*) 'SRF_FICE    ',SRF_FICE(1:nphy)
        write(lp,*) 'SRF_HICE    ',SRF_HICE(1:nphy)
        write(lp,*) 'SRF_HMLTFZ  ',SRF_HMLTFZ(1:nphy)
        write(lp,*) 'SRF_HSNW    ',SRF_HSNW(1:nphy)
        write(lp,*) 'SRF_IAGE    ',SRF_IAGE(1:nphy)
        write(lp,*) 'SRF_FICE    ',SRF_FICE(1:nphy)
        write(lp,*) 'SRF_HICE    ',SRF_HICE(1:nphy)
        write(lp,*) 'SRF_HMLTFZ  ',SRF_HMLTFZ(1:nphy)
        write(lp,*) 'SRF_HSNW    ',SRF_HSNW(1:nphy)
        write(lp,*) 'SRF_IAGE    ',SRF_IAGE(1:nphy)
        write(lp,*) 'SRF_LIP     ',SRF_LIP(1:nphy)
        write(lp,*) 'SRF_MAXMLD  ',SRF_MAXMLD(1:nphy)
        write(lp,*) 'SRF_MLD     ',SRF_MLD(1:nphy)
        write(lp,*) 'SRF_MLDU    ',SRF_MLDU(1:nphy)
        write(lp,*) 'SRF_MLDV    ',SRF_MLDV(1:nphy)
        write(lp,*) 'SRF_MTY     ',SRF_MTY(1:nphy)
        write(lp,*) 'SRF_MXLU    ',SRF_MXLU(1:nphy)
        write(lp,*) 'SRF_MXLV    ',SRF_MXLV(1:nphy)
        write(lp,*) 'SRF_NSF     ',SRF_NSF(1:nphy)
        write(lp,*) 'SRF_RFIFLX  ',SRF_RFIFLX(1:nphy)
        write(lp,*) 'SRF_RNFFLX  ',SRF_RNFFLX(1:nphy)
        write(lp,*) 'SRF_SALFLX  ',SRF_SALFLX(1:nphy)
        write(lp,*) 'SRF_SALRLX  ',SRF_SALRLX(1:nphy)
        write(lp,*) 'SRF_SEALV   ',SRF_SEALV(1:nphy)
        write(lp,*) 'SRF_SFL     ',SRF_SFL(1:nphy)
        write(lp,*) 'SRF_SIGMX   ',SRF_SIGMX(1:nphy)
        write(lp,*) 'SRF_SOP     ',SRF_SOP(1:nphy)
        write(lp,*) 'SRF_SSS     ',SRF_SSS(1:nphy)
        write(lp,*) 'SRF_SST     ',SRF_SST(1:nphy)
        write(lp,*) 'SRF_SURFLX  ',SRF_SURFLX(1:nphy)
        write(lp,*) 'SRF_SURRLX  ',SRF_SURRLX(1:nphy)
        write(lp,*) 'SRF_SWA     ',SRF_SWA(1:nphy)
        write(lp,*) 'SRF_TAUX    ',SRF_TAUX(1:nphy)
        write(lp,*) 'SRF_TAUY    ',SRF_TAUY(1:nphy)
        write(lp,*) 'SRF_TICE    ',SRF_TICE(1:nphy)
        write(lp,*) 'SRF_TSRF    ',SRF_TSRF(1:nphy)
        write(lp,*) 'SRF_UB      ',SRF_UB(1:nphy)
        write(lp,*) 'SRF_UICE    ',SRF_UICE(1:nphy)
        write(lp,*) 'SRF_USTAR   ',SRF_USTAR(1:nphy)
        write(lp,*) 'SRF_VB      ',SRF_VB(1:nphy)
        write(lp,*) 'SRF_VICE    ',SRF_VICE(1:nphy)
        write(lp,*) 'SRF_ZTX     ',SRF_ZTX(1:nphy)
        write(lp,*) 'LYR_DIAFLX  ',LYR_DIAFLX(1:nphy)
        write(lp,*) 'LYR_DIFDIA  ',LYR_DIFDIA(1:nphy)
        write(lp,*) 'LYR_DIFINT  ',LYR_DIFINT(1:nphy)
        write(lp,*) 'LYR_DIFISO  ',LYR_DIFISO(1:nphy)
        write(lp,*) 'LYR_DP      ',LYR_DP(1:nphy)
        write(lp,*) 'LYR_DZ      ',LYR_DZ(1:nphy)
        write(lp,*) 'LYR_SALN    ',LYR_SALN(1:nphy)
        write(lp,*) 'LYR_TEMP    ',LYR_TEMP(1:nphy)
        write(lp,*) 'LYR_TRC     ',LYR_TRC(1:nphy)
        write(lp,*) 'LYR_UFLX    ',LYR_UFLX(1:nphy)
        write(lp,*) 'LYR_UTFLX   ',LYR_UTFLX(1:nphy)
        write(lp,*) 'LYR_USFLX   ',LYR_USFLX(1:nphy)
        write(lp,*) 'LYR_UMFLTD  ',LYR_UMFLTD(1:nphy)
        write(lp,*) 'LYR_UTFLTD  ',LYR_UTFLTD(1:nphy)
        write(lp,*) 'LYR_UTFLLD  ',LYR_UTFLLD(1:nphy)
        write(lp,*) 'LYR_USFLTD  ',LYR_USFLTD(1:nphy)
        write(lp,*) 'LYR_USFLLD  ',LYR_USFLLD(1:nphy)
        write(lp,*) 'LYR_UVEL    ',LYR_UVEL(1:nphy)
        write(lp,*) 'LYR_VFLX    ',LYR_VFLX(1:nphy)
        write(lp,*) 'LYR_VTFLX   ',LYR_VTFLX(1:nphy)
        write(lp,*) 'LYR_VSFLX   ',LYR_VSFLX(1:nphy)
        write(lp,*) 'LYR_VMFLTD  ',LYR_VMFLTD(1:nphy)
        write(lp,*) 'LYR_VTFLTD  ',LYR_VTFLTD(1:nphy)
        write(lp,*) 'LYR_VTFLLD  ',LYR_VTFLLD(1:nphy)
        write(lp,*) 'LYR_VSFLTD  ',LYR_VSFLTD(1:nphy)
        write(lp,*) 'LYR_VSFLLD  ',LYR_VSFLLD(1:nphy)
        write(lp,*) 'LYR_VVEL    ',LYR_VVEL(1:nphy)
        write(lp,*) 'LYR_WFLX    ',LYR_WFLX(1:nphy)
        write(lp,*) 'LYR_WFLX2   ',LYR_WFLX2(1:nphy)
        write(lp,*) 'LYR_TKE     ',LYR_TKE(1:nphy)
        write(lp,*) 'LYR_GLS_PSI ',LYR_GLS_PSI(1:nphy)
        write(lp,*) 'LYR_IDLAGE  ',LYR_IDLAGE(1:nphy)
        write(lp,*) 'LVL_SALN    ',LVL_SALN(1:nphy)
        write(lp,*) 'LVL_TEMP    ',LVL_TEMP(1:nphy)
        write(lp,*) 'LVL_TRC     ',LVL_TRC(1:nphy)
        write(lp,*) 'LVL_UFLX    ',LVL_UFLX(1:nphy)
        write(lp,*) 'LVL_UTFLX   ',LVL_UTFLX(1:nphy)
        write(lp,*) 'LVL_USFLX   ',LVL_USFLX(1:nphy)
        write(lp,*) 'LVL_UMFLTD  ',LVL_UMFLTD(1:nphy)
        write(lp,*) 'LVL_UTFLTD  ',LVL_UTFLTD(1:nphy)
        write(lp,*) 'LVL_UTFLLD  ',LVL_UTFLLD(1:nphy)
        write(lp,*) 'LVL_USFLTD  ',LVL_USFLTD(1:nphy)
        write(lp,*) 'LVL_USFLLD  ',LVL_USFLLD(1:nphy)
        write(lp,*) 'LVL_UVEL    ',LVL_UVEL(1:nphy)
        write(lp,*) 'LVL_VFLX    ',LVL_VFLX(1:nphy)
        write(lp,*) 'LVL_VTFLX   ',LVL_VTFLX(1:nphy)
        write(lp,*) 'LVL_VSFLX   ',LVL_VSFLX(1:nphy)
        write(lp,*) 'LVL_VMFLTD  ',LVL_VMFLTD(1:nphy)
        write(lp,*) 'LVL_VTFLTD  ',LVL_VTFLTD(1:nphy)
        write(lp,*) 'LVL_VTFLLD  ',LVL_VTFLLD(1:nphy)
        write(lp,*) 'LVL_VSFLTD  ',LVL_VSFLTD(1:nphy)
        write(lp,*) 'LVL_VSFLLD  ',LVL_VSFLLD(1:nphy)
        write(lp,*) 'LVL_VVEL    ',LVL_VVEL(1:nphy)
        write(lp,*) 'LVL_WFLX    ',LVL_WFLX(1:nphy)
        write(lp,*) 'LVL_WFLX2   ',LVL_WFLX2(1:nphy)
        write(lp,*) 'LVL_TKE     ',LVL_TKE(1:nphy)
        write(lp,*) 'LVL_GLS_PSI ',LVL_GLS_PSI(1:nphy)
        write(lp,*) 'LVL_IDLAGE  ',LVL_IDLAGE(1:nphy)
        write(lp,*) 'MSC_MMFLXL  ',MSC_MMFLXL(1:nphy)
        write(lp,*) 'MSC_MMFLXD  ',MSC_MMFLXD(1:nphy)
        write(lp,*) 'MSC_MMFTDL  ',MSC_MMFTDL(1:nphy)
        write(lp,*) 'MSC_MMFTDD  ',MSC_MMFTDD(1:nphy)
        write(lp,*) 'MSC_MHFLX   ',MSC_MHFLX(1:nphy)
        write(lp,*) 'MSC_MHFTD   ',MSC_MHFTD(1:nphy)
        write(lp,*) 'MSC_MHFLD   ',MSC_MHFLD(1:nphy)
        write(lp,*) 'MSC_MSFLX   ',MSC_MSFLX(1:nphy)
        write(lp,*) 'MSC_MSFTD   ',MSC_MSFTD(1:nphy)
        write(lp,*) 'MSC_MSFLD   ',MSC_MSFLD(1:nphy)
        write(lp,*) 'MSC_VOLTR   ',MSC_VOLTR(1:nphy)
        write(lp,*) ' '
      endif 
c
c --- convert integer dates 
      nyear=sign(abs(idate)/10000,idate)
      nmonth=abs(idate)/100-abs(nyear)*100
      nday=abs(idate)-abs(nyear)*10000-nmonth*100
      nyear0=sign(abs(idate0)/10000,idate0)
      nmonth0=abs(idate0)/100-abs(nyear0)*100
      nday0=abs(idate0)-abs(nyear0)*10000-nmonth0*100
c
c --- set path lengths
      path_len=1
      do while (path_len.lt.80.and.path(path_len:path_len).ne.' ')
        path_len=path_len+1
      enddo
      path_len=path_len-1
      path1_len=1
      do while (path1_len.lt.80.and.path1(path1_len:path1_len).ne.' ')
        path1_len=path1_len+1
      enddo
      path1_len=path1_len-1
      path2_len=1
      do while (path2_len.lt.80.and.path2(path2_len:path2_len).ne.' ')
        path2_len=path2_len+1
      enddo
      path2_len=path2_len-1
#ifdef SYN
      atm_path_len=1
      do while (atm_path_len.lt.80.and.
     .          atm_path(atm_path_len:atm_path_len).ne.' ')
        atm_path_len=atm_path_len+1
      enddo
      atm_path_len=atm_path_len-1
#endif
c
#ifdef CCSMCOUPLED
c --- override some information read from the micom namelist with
c --- information received from the coupler and from the restart file
c --- pointer if it exists
      runid(1:80)=runid_mct(1:80)
      inquire(file='rpointer.ocn',exist=fexist)
      if (fexist) then
        open (unit=nfu,file='rpointer.ocn')
        read (nfu,'(a)') rstfnm
        close (unit=nfu)
        time=0.
        if (mnproc.eq.1) then
          inquire(file=path2(1:path2_len)//rstfnm,exist=fexist)
          if (fexist) then
            call ncfopn(path2(1:path2_len)//rstfnm,'r',' ',1)
            call ncgetr('time',time)
          else
            write (lp,*) 'Could not find proper restart file!'
            call xchalt('(rdlim)')
            stop '(rdlim)' 
          endif
        endif
        call xcmaxr(time)
        nday1=nint(time)
        if (runtype_mct.ne.'initial') then
          n=1
          do while (n.lt.120.and.rstfnm(n:n).ne.' ')
            n=n+1
          enddo
          n=n-1
          read (rstfnm(n-18:n-15),'(i4)') nyear
          read (rstfnm(n-13:n-12),'(i2)') nmonth
          read (rstfnm(n-10:n-9 ),'(i2)') nday
        endif
      endif
#endif
c
c --- set experiment id length
      runid_len=1
      do while (runid_len.lt.80.and.runid(runid_len:runid_len).ne.' ')
        runid_len=runid_len+1
      enddo
      runid_len=runid_len-1
c
      if (nday1.lt.0.or.nday2.lt.0) then
        if (mnproc.eq.1) then
          write (lp,*) 'Integrations days must be positive!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (nday2.lt.nday1) then
        if (mnproc.eq.1) then
          write (lp,'(2a)') ' First day of integration must be less ',
     .      'than or equal to last day!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (nday1.eq.0.and.
     .    (nyear.ne.nyear0.or.nmonth.ne.nmonth0.or.nday.ne.nday0)) then
        if (mnproc.eq.1) then
          write (lp,'(2a)') ' When first integration day is zero, ',
     .      'model date and initial experiment date '
          write (lp,'(a)') ' must be equal!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (nyear.lt.nyear0.or.(nyear.eq.nyear0.and.nmonth.lt.nmonth0).or.
     .  (nyear.eq.nyear0.and.nmonth.eq.nmonth0.and.nday.lt.nday0)) then
        if (mnproc.eq.1) then
          write (lp,'(2a)') ' Model date must be greater or equal to ',
     .      'initial experiment date!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (mod(86400./baclin+epsil,1.).gt.2.*epsil) then
        if (mnproc.eq.1) then
          write (lp,'(2a)') ' Must have an integer number of ',
     .     'baroclinic time steps pr. day!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (trxday.eq.0.and.ditflx) then
        if (mnproc.eq.1) then
          write (lp,*) 'trxday=0. and ditflx=.true.. Inconsistent!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (srxday.eq.0.and.disflx) then
        if (mnproc.eq.1) then
          write (lp,*) 'srxday=0. and disflx=.true.. Inconsistent!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (srxsrc.ne.'PHC3.0'.and.srxsrc.ne.'CORE') then
        if (mnproc.eq.1) then
          write (lp,*) 'srxsrc must be either ''PHC3.0'' or ''CORE''!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
c --- get time steps in a day
      nstep_in_day=nint(86400./baclin)
c
c --- represent time between restarts in time steps
      rstmon=.false.
      rstann=.false.
      if (nint(rstfrq).eq.30) then
        rstmon=.true.
      elseif (nint(rstfrq).eq.365) then
        rstann=.true.
      endif
      rstfrq=nstep_in_day*max(1.,rstfrq)
c
c --- represent time between diagnostics in time steps
      do n=1,nphy 
        GLB_FILEFREQ(n)=max(GLB_AVEPERIO(n),GLB_FILEFREQ(n))
c
        if (GLB_AVEPERIO(n).lt.0) then
          diagfq_phy(n)=-real(nstep_in_day)/GLB_AVEPERIO(n)
        else
          diagfq_phy(n)=nstep_in_day*max(1,GLB_AVEPERIO(n))
        endif
        diagmon_phy(n)=.false.
        diagann_phy(n)=.false.
        if (GLB_AVEPERIO(n).eq.30) then
          diagmon_phy(n)=.true.
        elseif (GLB_AVEPERIO(n).eq.365) then
          diagann_phy(n)=.true.
        endif
c
        if (GLB_FILEFREQ(n).lt.0) then
          filefq_phy(n)=-real(nstep_in_day)/GLB_FILEFREQ(n)
        else
          filefq_phy(n)=nstep_in_day*max(1,GLB_FILEFREQ(n))
        endif
        filemon_phy(n)=.false.
        fileann_phy(n)=.false.
        if (GLB_FILEFREQ(n).eq.30) then
          filemon_phy(n)=.true.
        elseif (GLB_FILEFREQ(n).eq.365) then
          fileann_phy(n)=.true.
        endif
      enddo
c
      if (mnproc.eq.1) then
        write (lp,'(a,i10)')   'nday1:             ',nday1
        write (lp,'(a,i10)')   'nday2:             ',nday2
        write (lp,'(a,i4.4,a,i2.2,a,i2.2)') 'nyear ,nmonth ,nday : ',
     .                             nyear ,'.',nmonth ,'.',nday
        write (lp,'(a,i4.4,a,i2.2,a,i2.2)') 'nyear0,nmonth0,nday0: ',
     .                             nyear0,'.',nmonth0,'.',nday0
        write (lp,'(a,a60)')   'runid:             ',runid
        write (lp,'(a,f10.4)') 'baclin:            ',baclin
        write (lp,'(a,f10.4)') 'batrop:            ',batrop
        write (lp,'(a,f10.4)') 'mdv2hi:            ',mdv2hi
        write (lp,'(a,f10.4)') 'mdv2lo:            ',mdv2lo
        write (lp,'(a,f10.4)') 'mdv4hi:            ',mdv4hi
        write (lp,'(a,f10.4)') 'mdv4lo:            ',mdv4lo
        write (lp,'(a,f10.4)') 'mdc2hi:            ',mdc2hi
        write (lp,'(a,f10.4)') 'mdc2lo:            ',mdc2lo
        write (lp,'(a,f10.4)') 'vsc2hi:            ',vsc2hi
        write (lp,'(a,f10.4)') 'vsc2lo:            ',vsc2lo
        write (lp,'(a,f10.4)') 'vsc4hi:            ',vsc4hi
        write (lp,'(a,f10.4)') 'vsc4lo:            ',vsc4lo
        write (lp,'(a,f10.4)') 'cbar:              ',cbar
        write (lp,'(a,f10.4)') 'cb:                ',cb
        write (lp,'(a,f10.4)') 'cwbdts:            ',cwbdts
        write (lp,'(a,f10.4)') 'cwbdls:            ',cwbdls
        write (lp,'(a,f10.4)') 'trxday:            ',trxday
        write (lp,'(a,f10.4)') 'srxday:            ',srxday
        write (lp,'(a,f10.4)') 'trxdpt:            ',trxdpt
        write (lp,'(a,f10.4)') 'srxdpt:            ',srxdpt
        write (lp,'(a,f10.4)') 'trxlim:            ',trxlim
        write (lp,'(a,f10.4)') 'srxlim:            ',srxlim
        write (lp,'(a,l10)')   'aptflx:            ',aptflx
        write (lp,'(a,l10)')   'apsflx:            ',apsflx
        write (lp,'(a,l10)')   'ditflx:            ',ditflx
        write (lp,'(a,l10)')   'disflx:            ',disflx
        write (lp,'(a,a50)')   'srxsrc:            ',srxsrc
        write (lp,'(a,l10)')   'smtfrc:            ',smtfrc
        write (lp,'(a,l10)')   'sprfac:            ',sprfac
        write (lp,'(a,a60)')   'path:      ',path
        write (lp,'(a,a60)')   'path1:     ',path1
        write (lp,'(a,a60)')   'path2:     ',path2
#ifdef SYN
        write (lp,'(a,a60)')   'atm_path:  ',atm_path
#endif
c
        write (lp,101) mdv2hi,mdv2lo,mdc2hi,mdc2lo,vsc2hi,vsc2lo
 101    format (' turb. flux parameters:',1p/
     .    ' mdv2hi =',e9.2,' mdv2lo =',e9.2/
     .    ' mdc2hi =',e9.2,' mdc2lo =',e9.2/
     .    ' vsc2hi =',e9.2,' vsc2lo =',e9.2)
        write (lp,'(a,f10.4)') 'rstfrq:            ',rstfrq
        write (lp,'(a,l10)')   'rstmon:            ',rstmon
        write (lp,'(a,l10)')   'rstann:            ',rstann
        call flush(lp)
      endif
c
c --- 'lstep' = number of barotropic time steps per baroclinic time step.
c --- lstep   m u s t   be even.
      lstep=2*ceiling(.5*baclin/batrop)
      dlt=baclin/lstep
      if (mnproc.eq.1) then
        write (lp,'(i4,'' barotropic steps per baroclinic time step'')')
     .    lstep
        call flush(lp)
      endif
c
c --- model is to be integrated from time step 'nstep1' to 'nstep2'
      nstep1=nday1*nstep_in_day
      nstep2=nday2*nstep_in_day
c
      if (csdiag) then
        nstep2=nstep1+2
      endif
c
      return
      end
