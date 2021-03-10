      program micom_ensemble_init
c
#undef DIAG
#define DIAG
c
      use mod_xc

#ifdef DIAG
      use mod_nctools
#endif
c
      implicit none
c
#include "common_blocks.h"
#include "common_assim.h"
#include "netcdf.inc"
      integer*4, external :: iargc
c
      real, dimension(itdm,jtdm) :: tmp2d,tmp2dsum  
      integer, dimension(4) :: start,count
      integer i,j,k,l,status,dimid,varid,ncid
      real q,dps,dpsum
      real mass1,saln1,dfreshwater1,dsalt1
      logical compensate_ice_freshwater
      character(len=80) :: oldfile
      character(len=3) :: cmem
      character(len=8) :: ctmp
      integer imem,imemmax,imemmax2,imemfirst,iarg
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: iuu,ivv,iqq

c
      real swtfrz,sig,sofsig,getlpi
      external swtfrz,sig,sofsig,getlpi

c
c --- ------------------------------------------------------------------
c --- Initialize SPMD processing
c --- ------------------------------------------------------------------
c
      call xcspmd(imem,imemmax)
      iarg=iargc()
      if (iarg==1 .or. iarg==2 ) then
         call getarg(iarg,ctmp)
         read(ctmp,*) imemmax2
         imemfirst=1
         if (iarg==2) then 
           call getarg(1,ctmp)
           read(ctmp,*) imemfirst
           imemmax2=imemmax2-imemfirst+1 
         end if 
         write(cmem,'(i3.3)') imem+imemfirst-1
         if (imemmax.ne.imemmax2) then 
           print *
           print *,'Wrong usage: '
           print *,'Number of tasks inconsistent with number of members'
           print *,'Expected number of members=',imemmax 
           print *,'From input argument=',imemmax2
           call exit(1)
         !else 
         !  print *,'Member=',imem,' Task=',mnproc
         endif 
      else
         print *
         print *,'Wrong usage: '
         print *,'   micom_ensemble_init ensemble_size'
         print *,'   or'
         print *,'   micom_ensemble_init first_member last_member'
         call exit(1)
      endif

c
c --- ------------------------------------------------------------------
c --- Initialize check sum algorithm
c --- ------------------------------------------------------------------
c
      call crcinit
c
c --- set variables usually set in namelist
      baclin=4800.
      batrop=96.
      path='./'
      path_len=2
c --- ------------------------------------------------------------------
c --- Initialize the geographic environment
c --- ------------------------------------------------------------------
c
      call inigeo
c
c --- ------------------------------------------------------------------
c --- Initialize various arrays
c --- ------------------------------------------------------------------
c
      call inivar
c
c --- ------------------------------------------------------------------
c --- Define coefficients for equation of state functions
c --- ------------------------------------------------------------------
c
      call eosini
c
c --- ------------------------------------------------------------------
c --- mass field initialization
c --- ------------------------------------------------------------------
c
      oldfile='forecast'//cmem//'.nc'
#ifdef DEnOI
      oldfile='dyna_MICOM_'//cmem//'.nc'
#endif
#ifdef DOSA
      oldfile='previous'//cmem//'.nc'
#endif
      if (mnproc.eq.1) then
c        write (lp,'(2a)') ' reading initial condition from ',
c     .                    path(1:path_len)//trim(oldfile)
c        call flush(lp)
c
c --- - open netcdf file
        status=nf_open(path(1:path_len)//trim(oldfile),nf_nowrite,
     .                 ncid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_open:'//oldfile,
     .                     nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
c
c --- - check dimensions
        status=nf_inq_dimid(ncid,'x',dimid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimid: x: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
        status=nf_inq_dimlen(ncid,dimid,i)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimlen: x: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
        status=nf_inq_dimid(ncid,'y',dimid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimid: y: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
        status=nf_inq_dimlen(ncid,dimid,j)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimlen: y: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
        if (i.ne.itdm.or.j.ne.jtdm) then
          write (lp,*) 'wrong dimensions in '//trim(oldfile)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
c
      endif
c
      start(1)=1
      start(2)=1
      start(4)=1
      count(1)=itdm
      count(2)=jtdm
      count(3)=1
      count(4)=1
c
c --- read reference potential density
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'sigmar',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: sigmar: ',nf_strerror(status)
          call xchalt('(micom_serial_init)')
                 stop '(micom_serial_init)'
        endif
      endif
      do k=1,kk
         if (mnproc.eq.1) then
            start(3)=k
            status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
            if (status.ne.nf_noerr) then
              write(lp,'(2a)') 'nf_get_vara_double: sigmar: ',
     .                      nf_strerror(status)
              call xchalt('(micom_serial)')
              stop '(micom_serial)'
            endif
          endif
       call xcaput(tmp2d,sigmar(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- read potential temperature
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'temp',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: temp: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: temp: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,temp(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- read salinity
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'saln',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: saln: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: saln: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,saln(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- read layer thickness
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'dp',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: dp: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: dp: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,dp(1-nbdy,1-nbdy,k),1)
      enddo
c --- read pbot
      if (mnproc.eq.1) then
        k=1
        status=nf_inq_varid(ncid,'pb',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: dp: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      k=1
      if (mnproc.eq.1) then
        start(3)=1
        status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_vara_double: pb: ',
     .                       nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      call xcaput(tmp2d,pb(1-nbdy,1-nbdy,k),1)


c --- read ub
      k=1
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'ub',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: ub: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      if (mnproc.eq.1) then
        start(3)=1
        status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_vara_double: ub: ',
     .                       nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      call xcaput(tmp2d,ub(1-nbdy,1-nbdy,k),1)
c --- read vb
      k=1
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vb',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vb: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      if (mnproc.eq.1) then
        start(3)=1
        status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_vara_double: ub: ',
     .                       nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      call xcaput(tmp2d,vb(1-nbdy,1-nbdy,k),1)
c --- read ubflx
      k=1
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'ubflx',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: ubflx: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      if (mnproc.eq.1) then
        start(3)=1
        status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_vara_double: ubflx: ',
     .                       nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      call xcaput(tmp2d,ubflx(1-nbdy,1-nbdy,k),1)
c --- read vbflx
      k=1
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vbflx',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vbflx: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      if (mnproc.eq.1) then
        start(3)=1
        status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_vara_double: ub: ',
     .                       nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      call xcaput(tmp2d,vbflx(1-nbdy,1-nbdy,k),1)
c --- read ubflxs
      k=1
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'ubflxs',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: ubflxs: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      if (mnproc.eq.1) then
        start(3)=1
        status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_vara_double: ubflxs: ',
     .                       nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      call xcaput(tmp2d,ubflxs(1-nbdy,1-nbdy,k),1)
c --- read vbflxs
      k=1
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vbflxs',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vbflxs: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      if (mnproc.eq.1) then
        start(3)=1
        status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_vara_double: vbflxs: ',
     .                       nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      call xcaput(tmp2d,vbflxs(1-nbdy,1-nbdy,k),1)


c --- read u
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'u',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: u: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: u: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,u(1-nbdy,1-nbdy,k),1)
      enddo
c --- read v
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'v',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: v: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: v: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,v(1-nbdy,1-nbdy,k),1)
      enddo
c --- read uflx
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'uflx',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: uflx: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: uflx: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,uflx(1-nbdy,1-nbdy,k),1)
      enddo

c --- read vflx
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vflx',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vflx: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: vflx: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,vflx(1-nbdy,1-nbdy,k),1)
      enddo
c --- read utflx
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'utflx',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: utflx: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: utflx: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,utflx(1-nbdy,1-nbdy,k),1)
      enddo

c --- read vtflx
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vtflx',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vtflx: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: vtflx: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,vtflx(1-nbdy,1-nbdy,k),1)
      enddo
c --- read usflx
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'usflx',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: usflx: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: usflx: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,usflx(1-nbdy,1-nbdy,k),1)
      enddo

c --- read vsflx
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vsflx',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vsflx: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: vsflx: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,vsflx(1-nbdy,1-nbdy,k),1)
      enddo

c ----- read umfltd
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'umfltd',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: umfltd: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: umfltd: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,umfltd(1-nbdy,1-nbdy,k),1)
      enddo

c ----- read vmfltd
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vmfltd',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vmfltd: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: vmfltd: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,vmfltd(1-nbdy,1-nbdy,k),1)
      enddo

c ----- read utfltd
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'utfltd',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: utfltd: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: utfltd: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,utfltd(1-nbdy,1-nbdy,k),1)
      enddo

c ----- read vtfltd
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vtfltd',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vtfltd: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: vtfltd: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,vtfltd(1-nbdy,1-nbdy,k),1)
      enddo
c ----- read utflld
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'utflld',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: utflld: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: utflld: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,utflld(1-nbdy,1-nbdy,k),1)
      enddo

c ----- read vtflld
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vtflld',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vtflld: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: vtflld: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,vtflld(1-nbdy,1-nbdy,k),1)
      enddo
c ----- read usfltd
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'usfltd',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: usfltd: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: usfltd: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,usfltd(1-nbdy,1-nbdy,k),1)
      enddo

c ----- read vsfltd
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vsfltd',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vsfltd: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: vsfltd: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,vsfltd(1-nbdy,1-nbdy,k),1)
      enddo
c ----- read usflld
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'usflld',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: usflld: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: usflld: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,usflld(1-nbdy,1-nbdy,k),1)
      enddo

c ----- read vsflld
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vsflld',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vsflld: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: vsflld: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
        endif
        call xcaput(tmp2d,vsflld(1-nbdy,1-nbdy,k),1)
      enddo
#ifdef DOSA
c --- read regc
      k=1
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'regc',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: regc: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      if (mnproc.eq.1) then
        start(3)=1
        status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_vara_double: regc: ',
     .                       nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      call xcaput(tmp2d,regc(1-nbdy,1-nbdy),1)

c --- read rrm0
      k=1
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'rrm0',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: rrm0: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      if (mnproc.eq.1) then
        start(3)=1
        status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_vara_double: rrm0: ',
     .                       nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      call xcaput(tmp2d,rrm0(1-nbdy,1-nbdy),1)

c --- read rbdmc2
      k=1
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'rbdmc2',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: rbdmc2: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      if (mnproc.eq.1) then
        start(3)=1
        status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_vara_double: rrm0: ',
     .                       nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
      endif
      call xcaput(tmp2d,rbdmc2(1-nbdy,1-nbdy),1)
#endif
      if (mnproc.eq.1) then
        status=nf_close(ncid)
      endif       
c
c --- read ice volume before and after assimilation 
      inquire(file='compensate_ice_freshwater',
     .  exist=compensate_ice_freshwater)
      if (compensate_ice_freshwater) then 
      if (mnproc.eq.1) then
        status=nf_open(path(1:path_len)//'viceold'//cmem//'.nc',
     .    nf_nowrite,ncid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_open:'//oldfile,
     .                     nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
        status=nf_inq_varid(ncid,'vicen',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vicen: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
        tmp2dsum=0. 
        do k=1,5
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: vicen: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
          tmp2dsum=tmp2dsum+tmp2d 
        enddo
        status=nf_close(ncid)
      endif
      call xcaput(tmp2dsum,viceold,1)
c
      if (mnproc.eq.1) then
        status=nf_open(path(1:path_len)//'forecast_ice'//cmem//'.nc',
     .    nf_nowrite,ncid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_open:'//oldfile,
     .                     nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
        status=nf_inq_varid(ncid,'vicen',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vicen: ',nf_strerror(status)
          call xchalt('(micom_serial)')
                 stop '(micom_serial)'
        endif
        tmp2dsum=0.
        do k=1,5
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: vicen: ',
     .                       nf_strerror(status)
            call xchalt('(micom_serial)')
                   stop '(micom_serial)'
          endif
          tmp2dsum=tmp2dsum+tmp2d
        enddo
        status=nf_close(ncid)
      endif
      call xcaput(tmp2dsum,vicenew,1)
      endif 
c
CCCCCC Start correction
      call xctilr(dp, 1,kk, 1,1, halo_ps)
      call xctilr(pb, 1,2, nbdy,nbdy, halo_ps)
     
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
           dpsum=0.
C When not updating the full water colum, sum(DP)!=PB
           do k=1,kk
             dpsum=dpsum+ dp(i,j,k) 
           enddo
           do k=1,kk
             dp(i,j,k)= dp(i,j,k)*pb(i,j,1)/dpsum
           enddo
C          if less that 20 m we split half-half
C          Otherwise 1st=10m second the rest
C          Cannot have MLD < 5m take from ater below
           dpsum=dp(i,j,1)+dp(i,j,2)
           if (dpsum .le. 20.*98060.) then
              if (dpsum.le.(5.*98060.)) then
                 dp(i,j,1)=2.5*98060.
                 dp(i,j,2)=2.5*98060.
                 dp(i,j,3) = dp(i,j,3) + (dpsum-(5.*98060.))
                 dpsum=5.*98060.
              else
                 dp(i,j,1)=dpsum/2.
                 dp(i,j,2)=dpsum/2.
              endif
           else
              dp(i,j,1)=10.*98060.
              dp(i,j,2)=dpsum-10.*98060.
           end if
           do k = 3, kk-1
               dp(i,j,k+1) = dp(i,j,k+1) + min(0.0,dp(i,j,k))
               dp(i,j,k) = max(dp(i,j,k),0.0)
           end do
CC               !!! Go backwards to fix lowermost layer.
           do k = kk, 2, -1
            dp(i,j,k-1) = dp(i,j,k-1) + min(0.0,dp(i,j,k))
            dp(i,j,k)   =   max(dp(i,j,k),0.0)
           end do 
           do k=1,kk
             dp(i,j,k+kk)= dp(i,j,k)
          enddo
        enddo
        enddo
      enddo
c
      call xctilr(temp, 1,kk, 1,1, halo_ps)
      call xctilr(saln, 1,kk, 1,1, halo_ps)
      call xctilr(sigma, 1,kk, 1,1, halo_ps)

c --- Sea ice update freshwater compensation - begin 
      if (compensate_ice_freshwater) then 
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c           mass in kg m-2 of first layer 
c           mass1=dp(i,j,1)*98060.*1000. ! bug 
            mass1=dp(i,j,1)/98060.*1000. 
c           salinity in g kg-1 
            saln1=saln(i,j,1)
c           freshwater to be added in kg m-2 (ice units=1000 kg m-2)
            dfreshwater1=-(vicenew(i,j)-viceold(i,j))*1000.
c           resulting salinity change in g kg-1 
            dsalt1=-dfreshwater1/mass1*saln1
c           limit change to below 1 psu 
            dsalt1=sign(min(abs(dsalt1),1.),dsalt1)              
            saln(i,j,1)=saln1+dsalt1
          enddo
          enddo
        enddo
      endif 
c --- Sea ice update freshwater compensation - end 
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          do k=1,2
C            temp(i,j,k)=max(swtfrz(saln(i,j,k)),temp(i,j,k))
            temp(i,j,k+kk)=temp(i,j,k)
            saln(i,j,k+kk)=saln(i,j,k)
            sigma(i,j,k)=sig(temp(i,j,k),saln(i,j,k))
            sigma(i,j,k+kk)=sigma(i,j,k)
          enddo
          do k=3,kk
C            temp(i,j,k)=max(temmin(i,j,k),temp(i,j,k))
C            saln(i,j,k   )=sofsig(sigmar(i,j,k),temp(i,j,k))
            temp(i,j,k+kk)=temp(i,j,k)
            saln(i,j,k+kk)=saln(i,j,k)
            sigma(i,j,k)=sig(temp(i,j,k),saln(i,j,k))
            sigma(i,j,k+kk)=sigma(i,j,k)
          enddo
          ustar(i,j)=0.
          buoyfl(i,j)=0.
        enddo
        enddo
      enddo
c
c$OMP PARALLEL DO PRIVATE(kn)
      do j=0,jj+1
        do k=1,kk
          !kn=k+nn
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
          enddo
          enddo
        enddo
      enddo

      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          phi(i,j,kk+1)=-g*depths(i,j)*1.e2
        enddo
        enddo
      enddo


c$OMP END PARALLEL DO
      call xctilr(dp, 1,kk, 1,1, halo_ps)
      call xctilr(p, 1,kk+1, 2,2, halo_ps)
      call xctilr(phi(1-nbdy,1-nbdy,kk+1), 1,1, 1,1, halo_ps)
      call xctilr(pb, 1,2, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- set layer thickness and bottom pressure
c --- ------------------------------------------------------------------
c
      do j=0,jj+1
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          pb(i,j,2)=pb(i,j,1)
          pb_p(i,j)=pb(i,j,1)
        enddo
        enddo
      enddo
c
      do  j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          pbu(i,j,1)=min(pb(i,j,1),pb(i-1,j,1))
          pbu(i,j,2)=pbu(i,j,1)
          pbu_p(i,j)=pbu(i,j,1)
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          pbv(i,j,1)=min(pb(i,j,1),pb(i,j-1,1))
          pbv(i,j,2)=pbv(i,j,1)
          pbv_p(i,j)=pbv(i,j,1)
        enddo
        enddo
      enddo
c
c$OMP PARALLEL DO PRIVATE(q)
      do j=-1,jj+2
        do k=1,kk
          do l=1,isu(j)
          do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            q=min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,k)=
     .        max(0.,min(q,.5*(p(i,j,k+1)+p(i-1,j,k+1))) 
     .              -min(q,.5*(p(i,j,k  )+p(i-1,j,k  ))))
          enddo
          enddo
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            q=min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,k)=
     .        max(0.,min(q,.5*(p(i,j,k+1)+p(i,j-1,k+1)))
     .              -min(q,.5*(p(i,j,k  )+p(i,j-1,k  ))))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            dpu(i,j,k+kk)=dpu(i,j,k)
            pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,k)
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            dpv(i,j,k+kk)=dpv(i,j,k)
            pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- initialize potential vorticity of barotropic flow
c --- ------------------------------------------------------------------
c
      call xctilr(pvtrop, 1,2, 1,3, halo_qs)
      call xctilr(corioq, 1,1, nbdy,nbdy, halo_qs)
c$OMP PARALLEL DO PRIVATE(q)
      do j=0,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          q=2./(pb_p(i,j)+pb_p(i-1,j))
          pvtrop(i,j  ,1)=corioq(i,j  )*q
          pvtrop(i,j+1,1)=corioq(i,j+1)*q
          pvtrop(i,j  ,2)=pvtrop(i,j  ,1)
          pvtrop(i,j+1,2)=pvtrop(i,j+1,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(q)
      do j=1,jj
        do l=1,isv(j)
        do i=max(0,ifv(j,l)),min(ii,ilv(j,l))
          q=2./(pb_p(i,j)+pb_p(i,j-1))
          pvtrop(i  ,j,1)=corioq(i  ,j)*q
          pvtrop(i+1,j,1)=corioq(i+1,j)*q
          pvtrop(i  ,j,2)=pvtrop(i  ,j,1)
          pvtrop(i+1,j,2)=pvtrop(i+1,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isq(j)
        do i=max(1,ifq(j,l)),min(ii,ilq(j,l))
          pvtrop(i,j,1)=corioq(i,j)*4./(pb_p(i,j  )+pb_p(i-1,j  )
     .                                 +pb_p(i,j-1)+pb_p(i-1,j-1))
          pvtrop(i,j,2)=pvtrop(i,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
#ifdef DOSA
c$OMP PARALLEL DO
      call xctilr(regc, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(rbdmc2, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(rrm0, 1,1, nbdy,nbdy, halo_ps)
      do j=0,jj+1
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          regc(i,j)=min(max(regc(i,j),0.40),1.60)
          rbdmc2(i,j)=min(max(rbdmc2(i,j),0.09),0.19)
          rrm0(i,j)=min(max(rrm0(i,j),0.6),2.4)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
#endif
c
c --- ------------------------------------------------------------------
c --- initialize fields related to the pressure gradient force
c --- ------------------------------------------------------------------
c
C      call xctilr(pgfxm, 1,2, 1,2, halo_uv)
      call xctilr(pgfxm(1-nbdy,1-nbdy,1), 1,1, 1,2, halo_uv)
      call xctilr(pgfym(1-nbdy,1-nbdy,1), 1,1, 1,2, halo_vv)
C      call xctilr(pgfym, 1,2, 1,2, halo_uv)
      call pgforc(2,1,kk,0,kk+1,1)
      call xctilr(ubflxs(1-nbdy,1-nbdy,1), 1,1, 1,1, halo_uv)
      call xctilr(vbflxs(1-nbdy,1-nbdy,1), 1,1, 1,1, halo_vv)
      call xctilr(ubflxs_p, 1,2, 2,2, halo_uv)
      call xctilr(vbflxs_p, 1,2, 2,2, halo_vv)
C      call xctilr(ubflxs_p, 1, 2, nbdy,nbdy, halo_us)
C      call xctilr(vbflxs, 1,1, 1,1, halo_us)
C      call xctilr(vbflxs_p, 1, 2, nbdy,nbdy, halo_us)
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            pgfx (i,j,k+kk)=pgfx (i,j,k)
            u    (i,j,k+kk)=u    (i,j,k)
            uflx (i,j,k+kk)=uflx (i,j,k)
            utflx(i,j,k+kk)=utflx(i,j,k)
            usflx(i,j,k+kk)=usflx(i,j,k)
            umfltd(i,j,k+kk)=umfltd(i,j,k)
            utfltd(i,j,k+kk)=utfltd(i,j,k)
            utflld(i,j,k+kk)=utflld(i,j,k)
            usfltd(i,j,k+kk)=usfltd(i,j,k)
            usflld(i,j,k+kk)=usflld(i,j,k)
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            pgfy (i,j,k+kk)=pgfy (i,j,k)
            v    (i,j,k+kk)=v    (i,j,k)
            vflx (i,j,k+kk)=vflx (i,j,k)
            vtflx(i,j,k+kk)=vtflx(i,j,k)
            vsflx(i,j,k+kk)=vsflx(i,j,k)
            vmfltd(i,j,k+kk)=vmfltd(i,j,k)
            vtfltd(i,j,k+kk)=vtfltd(i,j,k)
            vtflld(i,j,k+kk)=vtflld(i,j,k)
            vsfltd(i,j,k+kk)=vsfltd(i,j,k)
            vsflld(i,j,k+kk)=vsflld(i,j,k)
          enddo
          enddo
        enddo
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          pgfxm(i,j,2)   =pgfxm(i,j,1)
          xixp(i,j,2)    =xixp(i,j,1)
          xixm(i,j,2)    =xixm(i,j,1)
          ub(i,j,2)      =ub(i,j,1)
          ubflx(i,j,2)   =ubflx(i,j,1)
          ubflx_mn(i,j,1)   =ubflx(i,j,1)
          ubflx_mn(i,j,2)   =ubflx(i,j,1)
          ubflxs(i,j,2)  =ubflxs(i,j,1)
          ubflxs(i,j,3)  =ubflxs(i,j,1)
          ubflxs_p(i,j,1)=ubflxs(i,j,1)
          ubflxs_p(i,j,2)=ubflxs(i,j,1)
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          pgfym(i,j,2) =pgfym(i,j,1)
          xiyp(i,j,2)  =xiyp(i,j,1)
          xiym(i,j,2)  =xiym(i,j,1)
          vb(i,j,2)    =vb(i,j,1)
          vbflx(i,j,2) =vbflx(i,j,1)
          vbflx_mn(i,j,1) =vbflx(i,j,1)
          vbflx_mn(i,j,2) =vbflx(i,j,1)
          vbflxs(i,j,2)=vbflxs(i,j,1)
          vbflxs(i,j,3)=vbflxs(i,j,1)
          vbflxs_p(i,j,1)=vbflxs(i,j,1)
          vbflxs_p(i,j,2)=vbflxs(i,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO


c
c --- ------------------------------------------------------------------
c --- define first physical interior layer
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(k,dps)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          k=3
          dps=0.
          do while (dp(i,j,k).lt.epsil)
            dps=dps+dp(i,j,k)
            dp(i,j,k)=0.
            dp(i,j,k+kk)=0.
            k=k+1
            if (k.gt.kk) exit
          enddo
          if (k.gt.kk) then
            dp(i,j,2)=dp(i,j,2)+dps
            dp(i,j,2+kk)=dp(i,j,2)
          else
            dp(i,j,k)=dp(i,j,k)+dps
            dp(i,j,k+kk)=dp(i,j,k)
          endif
          kfpla(i,j,1)=k
          kfpla(i,j,2)=k
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c

#ifdef DIAG
       do j=1,jj
          do i=1,ii
            if ((ip(i,j)+ip(i-1,j)).ge.1) then
              iuu(i,j)=1
            else
              iuu(i,j)=0
            endif
            if ((ip(i,j)+ip(i,j-1)).ge.1) then
              ivv(i,j)=1
            else
              ivv(i,j)=0
            endif
            if ((iu(i,j)+iv(i,j)+iu(i,j-1)+iv(i-1,j)).ge.1) then
              iqq(i,j)=1
            else
              iqq(i,j)=0
            endif
          enddo
        enddo
c      print *,'dumping on file ', trim(oldfile)
      call ncfopn(trim(oldfile),'wold','c',1)
      call ncwrtr('u'    ,'x y kk2 time',u     (1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('v'    ,'x y kk2 time',v     (1-nbdy,1-nbdy,1),ivv
     .  ,1,1.,0.,8)
      call ncwrtr('dp'   ,'x y kk2 time',dp    (1-nbdy,1-nbdy,1),ip
     . ,1,1.,0.,8)
      call ncwrtr('dpold','x y kk2 time',dp    (1-nbdy,1-nbdy,1),ip
     . ,1,1.,0.,8)
      call ncwrtr('temp' ,'x y kk2 time',temp  (1-nbdy,1-nbdy,1),ip
     . ,1,1.,0.,8)
      call ncwrtr('saln' ,'x y kk2 time',saln  (1-nbdy,1-nbdy,1),ip
     . ,1,1.,0.,8)
      call ncwrtr('sigma','x y kk2 time',sigma (1-nbdy,1-nbdy,1),ip
     . ,1,1.,0.,8)
      call ncwrtr('pgfx' ,'x y kk2 time',pgfx  (1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('pgfy' ,'x y kk2 time',pgfy  (1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('uflx' ,'x y kk2 time',uflx  (1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('vflx' ,'x y kk2 time',vflx  (1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('utflx','x y kk2 time',utflx (1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('vtflx','x y kk2 time',vtflx (1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('usflx','x y kk2 time',usflx (1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('vsflx','x y kk2 time',vsflx (1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('umfltd','x y kk2 time',umfltd(1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('vmfltd','x y kk2 time',vmfltd(1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('utfltd','x y kk2 time',utfltd(1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('vtfltd','x y kk2 time',vtfltd(1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('utflld','x y kk2 time',utflld(1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('vtflld','x y kk2 time',vtflld(1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('usfltd','x y kk2 time',usfltd(1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('vsfltd','x y kk2 time',vsfltd(1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('usflld','x y kk2 time',usflld(1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('vsflld','x y kk2 time',vsflld(1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
CC    Dump k2 record
      call ncwrtr('pb'      ,'x y k2 time',pb      (1-nbdy,1-nbdy,1),ip
     . ,1,1.,0.,8)
      call ncwrtr('pb_mn'   ,'x y k2 time',pb      (1-nbdy,1-nbdy,1),ip
     . ,1,1.,0.,8)
      call ncwrtr('pbu'     ,'x y k2 time',pbu     (1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('pbv'     ,'x y k2 time',pbv     (1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('ub'      ,'x y k2 time',ub      (1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('vb'      ,'x y k2 time',vb      (1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('ubflx'   ,'x y k2 time',ubflx   (1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('vbflx'   ,'x y k2 time',vbflx   (1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('ubflx_mn','x y k2 time',ubflx_mn(1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('ubflxs_p','x y k2 time',ubflxs_p(1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('vbflx'   ,'x y k2 time',vbflx   (1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('vbflx_mn','x y k2 time',vbflx_mn(1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('vbflxs_p','x y k2 time',vbflxs_p(1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('pvtrop'  ,'x y k2 time',pvtrop  (1-nbdy,1-nbdy,1),iqq
     . ,1,1.,0.,8)
      call ncwrtr('pgfxm'   ,'x y k2 time',pgfxm   (1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('pgfym'   ,'x y k2 time',pgfym   (1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('xixp'    ,'x y k2 time',xixp    (1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('xixm'    ,'x y k2 time',xixm    (1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('xiyp'    ,'x y k2 time',xiyp    (1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
      call ncwrtr('xiym'    ,'x y k2 time',xiym    (1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
CC    Dump single 4D 3 time level
      call ncwrtr('ubflxs' ,'x y k3 time',ubflxs(1-nbdy,1-nbdy,1),iuu
     . ,1,1.,0.,8)
      call ncwrtr('vbflxs' ,'x y k3 time',vbflxs(1-nbdy,1-nbdy,1),ivv
     . ,1,1.,0.,8)
CC    Dump integer 4D
      call ncwrti('kfpla' ,'x y k2 time',kfpla (1-nbdy,1-nbdy,1),ip,1)


CC    Dump single 3D
      call ncwrtr('pb_p'   ,'x y time',pb_p (1-nbdy,1-nbdy),ip
     . ,1,1.,0.,8)
      call ncwrtr('pbu_p'  ,'x y time',pbu_p(1-nbdy,1-nbdy),iuu
     . ,1,1.,0.,8)
      call ncwrtr('pbv_p'  ,'x y time',pbv_p(1-nbdy,1-nbdy),ivv
     . ,1,1.,0.,8)
#ifdef DOSA
      call ncwrtr('rrm0'  ,'x y time',rrm0(1-nbdy,1-nbdy),ip
     . ,1,1.,0.,8)
      call ncwrtr('regc'  ,'x y time',regc(1-nbdy,1-nbdy),ip
     . ,1,1.,0.,8)
      call ncwrtr('rbdmc2'  ,'x y time',rbdmc2(1-nbdy,1-nbdy),ip
     . ,1,1.,0.,8)
#endif
      call ncfcls
#endif
c
c      if (mnproc.eq.1) print *, 'end of micom_ensemble_init'
      call xcstop(' ')
      end
