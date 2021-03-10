      subroutine ininst
c
c --- ------------------------------------------------------------------
c --- Read limits file for coarse outer model
c --- ------------------------------------------------------------------
c
      implicit none
c
      integer nday1,nday2
c
#include "dimensions.h"
#include "common_blocks.h"
#include "common_clndr.h"
#include "param_nest.h"
#include "common_nest.h"
c
      integer ifield_o(idm_o,jdm_o),i,j,m,n,ibz,k,ifac
      real rfac
c
      open (unit=11,file='limits_outer',status='old',action='read')
      read (11,'(i4,x,i2,x,i2)') nyrotr,nmnotr,ndyotr ! Initial experiment date
      read (11,'(a3)') rngotr      ! Experiment name
      read (11,'(i3)' ) i          ! time between 2D diagnostics
      read (11,'(i3)' ) j          ! time between 3D diagnostics
      read (11,'(a)') pthotr       ! Path to diagnostic files
      close (unit=11)
c
      if (nyear.lt.nyrotr.or.(nyear.eq.nyrotr.and.nmonth.lt.nmnotr).or.
     .  (nyear.eq.nyrotr.and.nmonth.eq.nmnotr.and.nday.lt.ndyotr)) then
        write (lp,*)
     .'Model date must be greater or equal to initial outer model date!'
        stop
      endif
c
c --- trim the path
      n=1
      do while (pthotr(n:n).ne.':')
        n=n+1
        if (n.gt.80) then
          write (lp,*)
     .    'Place a colon after the end of useable characters in pthotr!'
          stop
        endif
      enddo
      do m=n,80
        pthotr(m:m)=' '
      enddo
c
c --- represent time between boundary data updates in time steps
      if (i.lt.0) then
        otr2dfq=-real(nstep_in_day)/i
      else
        otr2dfq=nstep_in_day*max(1,i)
      endif
      if (i.eq.30) then
        otr2dmon=.true.
      else
        otr2dmon=.false.
      endif
c
      if (j.lt.0) then
        otr3dfq=-real(nstep_in_day)/j
      else
        otr3dfq=nstep_in_day*max(1,j)
      endif
      if (j.eq.30) then
        otr3dmon=.true.
      else
        otr3dmon=.false.
      endif
c
      if (.not.otr3dmon.or..not.otr2dmon) then
        write (lp,*)
     .    'rdotr: only monthly boundary conditions are implemented!'
        stop
      endif
c
      write (lp,'(a,i4.4,a,i2.2,a,i2.2)') 'nyrotr,nmnotr,ndyotr: ',
     .                           nyrotr,'.',nmnotr,'.',ndyotr
      write (lp,'(a,a10)')   'rngotr:            ',rngotr
      if (otr2dmon) then
      write (lp,'(a,f10.4)') '2D diag. interval:    monthly'
      else
      write (lp,'(a,f10.4)') '2D diag. interval: ',otr2dfq/nstep_in_day
      endif
      if (otr3dmon) then
      write (lp,'(a,f10.4)') '3D diag. interval:    monthly'
      else
      write (lp,'(a,f10.4)') '3D diag. interval: ',otr3dfq/nstep_in_day
      endif
      write (lp,'(a,a60)')   'pthotr:     ',pthotr
c
c --- set path length
      lptotr=1
      do while (lptotr.lt.80.and.pthotr(lptotr:lptotr).ne.' ')
        lptotr=lptotr+1
      enddo
      lptotr=lptotr-1
c
c --- define weigths for nesting boundary zone
      do ibz=1,nbz
        wnst(ibz)=1.-tanh(cnst*real(ibz-1)/(nbz-1))
        write (lp,'(a,i2,a,f12.8)') 'wnst(',ibz,'): ',wnst(ibz)
      enddo
c
c --- read bathymetry for outer model
      write (lp,'(2a)') ' reading bathymetry file from ',
     .                  path(1:path_len)//'depth_o.forc'
      open (unit=9,file=path(1:path_len)//'depth_o.forc',status='old',
     .      action='read')
      read (9,'(4i8)') i,j,k,ifac
      if (i.ne.idm_o.or.j.ne.jdm_o.or.k.ne.1)
     .  stop '(wrong dimensions in depth.forc)'
      read (9,'(15i8)') ifield_o
      close (unit=9)
      rfac=1./real(ifac)
c$OMP PARALLEL DO
      do 6 j=1,jdm_o
      do 6 i=1,idm_o
      depths_o(i,j)=max(0.,real(ifield_o(i,j))*rfac)
      if (depths_o(i,j).gt.0.) depths_o(i,j)=max(50.,depths_o(i,j))
 6    continue
c$OMP END PARALLEL DO
c
c --- read boundary data from outer model
c
c$OMP PARALLEL DO
      do k=1,kk
        do j=1,jj
          do ibz=1,nbz
            temp_ow1(ibz,j,k)=huge
            temp_oe1(ibz,j,k)=huge
            saln_ow1(ibz,j,k)=huge
            saln_oe1(ibz,j,k)=huge
            utot_ow1(ibz,j,k)=huge
            utot_oe1(ibz,j,k)=huge
            vtot_ow1(ibz,j,k)=huge
            vtot_oe1(ibz,j,k)=huge
            dp_ow1(ibz,j,k)=huge
            dp_oe1(ibz,j,k)=huge
          enddo
        enddo
        do ibz=1,nbz
          do i=1,ii
            temp_os1(i,ibz,k)=huge
            temp_on1(i,ibz,k)=huge
            saln_os1(i,ibz,k)=huge
            saln_on1(i,ibz,k)=huge
            utot_os1(i,ibz,k)=huge
            utot_on1(i,ibz,k)=huge
            vtot_os1(i,ibz,k)=huge
            vtot_on1(i,ibz,k)=huge
            dp_os1(i,ibz,k)=huge
            dp_on1(i,ibz,k)=huge
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do j=1,jj
        ub_ow1(j)=huge
        ub_oe1(j)=huge
        vb_ow1(j)=huge
        vb_oe1(j)=huge
        pb_ow1(1,j)=huge
        pb_oe1(1,j)=huge
        pb_ow1(2,j)=huge
        pb_oe1(2,j)=huge
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do i=1,ii
        ub_os1(i)=huge
        ub_on1(i)=huge
        vb_os1(i)=huge
        vb_on1(i)=huge
        pb_os1(i,1)=huge
        pb_on1(i,1)=huge
        pb_os1(i,2)=huge
        pb_on1(i,2)=huge
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do j=1,jj
        do ibz=1,nbz
          ficem_ow1(ibz,j)=huge
          ficem_oe1(ibz,j)=huge
          hicem_ow1(ibz,j)=huge
          hicem_oe1(ibz,j)=huge
          hsnwm_ow1(ibz,j)=huge
          hsnwm_oe1(ibz,j)=huge
          iagem_ow1(ibz,j)=huge
          iagem_oe1(ibz,j)=huge
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do ibz=1,nbz
        do i=1,ii
          ficem_os1(i,ibz)=huge
          ficem_on1(i,ibz)=huge
          hicem_os1(i,ibz)=huge
          hicem_on1(i,ibz)=huge
          hsnwm_os1(i,ibz)=huge
          hsnwm_on1(i,ibz)=huge
          iagem_os1(i,ibz)=huge
          iagem_on1(i,ibz)=huge
        enddo
      enddo
c$OMP END PARALLEL DO
c
      write (*,*) 'nyear,nmonth,nday',nyear,nmonth,nday,nday_of_year
      if (otr3dmon) then
        i=2*31
      else
        i=max(1,nint(2.*otr3dfq/nstep_in_day))
      endif
      call pstdat(nyear,nmonth,nday,i)
      call inicln
      write (*,*) 'nyear,nmonth,nday',nyear,nmonth,nday,nday_of_year
      do nstep=nstep1-i*nstep_in_day+1,nstep1
        time=nstep*baclin/86400.
        call updcln(nstep)
        call rdotr
      enddo
      write (*,*) 'nyear,nmonth,nday',nyear,nmonth,nday,nday_of_year
c
      return
      end
