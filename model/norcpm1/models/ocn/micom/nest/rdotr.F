      subroutine rdotr
c
c --- ------------------------------------------------------------------
c --- Read boundary data from the coarse outer model
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "dimensions.h"
#include "common_blocks.h"
#include "common_clndr.h"
#include "param_nest.h"
#include "common_nest.h"
#include "netcdf.inc"
c
      integer ny,nm,nd,status,ncid,i,j,k,ibz
      real z_ow(nbz,jdm,kdm+1),z_oe(nbz,jdm,kdm+1),
     .     z_os(idm,nbz,kdm+1),z_on(idm,nbz,kdm+1),
     .     sealv_ow(nbz,jdm),sealv_oe(nbz,jdm),
     .     sealv_os(idm,nbz),sealv_on(idm,nbz)
      real phiup,philo,pup,plo,q
      character fname*25
c
      real getlpi
      external getlpi
c
      if (otr3dmon.and.
     .    abs(nday+mod(time,1.)-.5*nd_in_m(nmonth)).lt.epsil) then
c
      write (*,*) nday,nmonth,nyear,nstep,time
c --- - read 3D boundary data
c
c --- - copy boundary data to old time level
c$OMP PARALLEL DO
        do k=1,kk
          do j=1,jj
            do ibz=1,nbz
              temp_ow1(ibz,j,k)=temp_ow2(ibz,j,k)
              temp_oe1(ibz,j,k)=temp_oe2(ibz,j,k)
              saln_ow1(ibz,j,k)=saln_ow2(ibz,j,k)
              saln_oe1(ibz,j,k)=saln_oe2(ibz,j,k)
              utot_ow1(ibz,j,k)=utot_ow2(ibz,j,k)
              utot_oe1(ibz,j,k)=utot_oe2(ibz,j,k)
              vtot_ow1(ibz,j,k)=vtot_ow2(ibz,j,k)
              vtot_oe1(ibz,j,k)=vtot_oe2(ibz,j,k)
            enddo
          enddo
          do ibz=1,nbz
            do i=1,ii
              temp_os1(i,ibz,k)=temp_os2(i,ibz,k)
              temp_on1(i,ibz,k)=temp_on2(i,ibz,k)
              saln_os1(i,ibz,k)=saln_os2(i,ibz,k)
              saln_on1(i,ibz,k)=saln_on2(i,ibz,k)
              utot_os1(i,ibz,k)=utot_os2(i,ibz,k)
              utot_on1(i,ibz,k)=utot_on2(i,ibz,k)
              vtot_os1(i,ibz,k)=vtot_os2(i,ibz,k)
              vtot_on1(i,ibz,k)=vtot_on2(i,ibz,k)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c --- - get file name for 3D netCDF file
        ny=nyear
        nm=nmonth
        nd=nday
        if (otr3dmon) then
          nday=1
          if     (nmonth.eq.11) then
            nmonth=1
            nyear=nyear+1
          elseif (nmonth.eq.12) then
            nmonth=2
            nyear=nyear+1
          else
            nmonth=nmonth+2
          endif
        else
        endif
        call ncfnam(rngotr,'3d',time,otr3dfq/nstep_in_day,otr3dmon,
     .              fname)
        nyear=ny
        nmonth=nm
        nday=nd
        write (lp,'(2a)')
     .  'Reading 3D boundary conditions from: ',pthotr(1:lptotr)//fname
c
c --- - open netcdf file
        status=nf_open(pthotr(1:lptotr)//fname,nf_nowrite,ncid)
        if (status.ne.nf_noerr) then
          write(*,'(4a)') 'nf_open: ',pthotr(1:lptotr)//fname,',',
     .                    nf_strerror(status)
          stop
        endif
c
        call rdotr3(ncid,'temp',kk,temp_ow2,temp_oe2,temp_os2,temp_on2)
        call rdotr3(ncid,'saln',kk,saln_ow2,saln_oe2,saln_os2,saln_on2)
        call rdotr3(ncid,'utot',kk,utot_ow2,utot_oe2,utot_os2,utot_on2)
        call rdotr3(ncid,'vtot',kk,vtot_ow2,vtot_oe2,vtot_os2,vtot_on2)
        call rdotrz(ncid,depths_o,depths,z_ow,z_oe,z_os,z_on)
c
c --- - close netcdf file
        status=nf_close(ncid)
        if (status.ne.nf_noerr) then
          write(*,'(2a)') 'nf_close: ',nf_strerror(status)
          stop
        endif
c
c --- - get model time step for next 3D boundary update
        nst31=nst32
        if (otr3dmon) then
          nm=mod(nmonth,12)+1
          nst32=nint((time+.5*(nd_in_m(nmonth)+nd_in_m(nm)))
     .               *nstep_in_day)
        else
        endif
        write (*,*) 'nst31,nst32',nst31,nst32
c
c       open (10,file='rdotr.uf',form='unformatted')
c       write (10) idm,jdm,kdm,nbz
c       write (10) temp_ow1,temp_oe1,temp_os1,temp_on1,
c    .             temp_ow2,temp_oe2,temp_os2,temp_on2
c       write (10) saln_ow1,saln_oe1,saln_os1,saln_on1,
c    .             saln_ow2,saln_oe2,saln_os2,saln_on2
c       write (10) utot_ow1,utot_oe1,utot_os1,utot_on1,
c    .             utot_ow2,utot_oe2,utot_os2,utot_on2
c       write (10) vtot_ow1,vtot_oe1,vtot_os1,vtot_on1,
c    .             vtot_ow2,vtot_oe2,vtot_os2,vtot_on2
c
      endif
c
c --- read 2D boundary data
      if (otr2dmon.and.
     .    abs(nday+mod(time,1.)-.5*nd_in_m(nmonth)).lt.epsil) then
c
c --- - copy boundary data to old time level
c$OMP PARALLEL DO
        do k=1,kk
          do j=1,jj
            do ibz=1,nbz
              dp_ow1(ibz,j,k)=dp_ow2(ibz,j,k)
              dp_oe1(ibz,j,k)=dp_oe2(ibz,j,k)
            enddo
          enddo
          do ibz=1,nbz
            do i=1,ii
              dp_os1(i,ibz,k)=dp_os2(i,ibz,k)
              dp_on1(i,ibz,k)=dp_on2(i,ibz,k)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
        do j=1,jj
          ub_ow1(j)=ub_ow2(j)
          ub_oe1(j)=ub_oe2(j)
          vb_ow1(j)=vb_ow2(j)
          vb_oe1(j)=vb_oe2(j)
          pb_ow1(1,j)=pb_ow2(1,j)
          pb_oe1(1,j)=pb_oe2(1,j)
          pb_ow1(2,j)=pb_ow2(2,j)
          pb_oe1(2,j)=pb_oe2(2,j)
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
        do i=1,ii
          ub_os1(i)=ub_os2(i)
          ub_on1(i)=ub_on2(i)
          vb_os1(i)=vb_os2(i)
          vb_on1(i)=vb_on2(i)
          pb_os1(i,1)=pb_os2(i,1)
          pb_on1(i,1)=pb_on2(i,1)
          pb_os1(i,2)=pb_os2(i,2)
          pb_on1(i,2)=pb_on2(i,2)
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
        do j=1,jj
          do ibz=1,nbz
            ficem_ow1(ibz,j)=ficem_ow2(ibz,j)
            ficem_oe1(ibz,j)=ficem_oe2(ibz,j)
            hicem_ow1(ibz,j)=hicem_ow2(ibz,j)
            hicem_oe1(ibz,j)=hicem_oe2(ibz,j)
            hsnwm_ow1(ibz,j)=hsnwm_ow2(ibz,j)
            hsnwm_oe1(ibz,j)=hsnwm_oe2(ibz,j)
            iagem_ow1(ibz,j)=iagem_ow2(ibz,j)
            iagem_oe1(ibz,j)=iagem_oe2(ibz,j)
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
        do ibz=1,nbz
          do i=1,ii
            ficem_os1(i,ibz)=ficem_os2(i,ibz)
            ficem_on1(i,ibz)=ficem_on2(i,ibz)
            hicem_os1(i,ibz)=hicem_os2(i,ibz)
            hicem_on1(i,ibz)=hicem_on2(i,ibz)
            hsnwm_os1(i,ibz)=hsnwm_os2(i,ibz)
            hsnwm_on1(i,ibz)=hsnwm_on2(i,ibz)
            iagem_os1(i,ibz)=iagem_os2(i,ibz)
            iagem_on1(i,ibz)=iagem_on2(i,ibz)
          enddo
        enddo
c$OMP END PARALLEL DO
c
c --- - get file name for 2D netCDF file
        ny=nyear
        nm=nmonth
        nd=nday
        if (otr2dmon) then
          nday=1
          if     (nmonth.eq.11) then
            nmonth=1
            nyear=nyear+1
          elseif (nmonth.eq.12) then
            nmonth=2
            nyear=nyear+1
          else
            nmonth=nmonth+2
          endif
        else
        endif
        call ncfnam(rngotr,'2d',time,otr2dfq/nstep_in_day,otr2dmon,
     .              fname)
        nyear=ny
        nmonth=nm
        nday=nd
        write (lp,'(2a)')
     .  'Reading 2D boundary conditions from: ',pthotr(1:lptotr)//fname
c
c --- - open netcdf file
        status=nf_open(pthotr(1:lptotr)//fname,nf_nowrite,ncid)
        if (status.ne.nf_noerr) then
          write(*,'(4a)') 'nf_open: ',pthotr(1:lptotr)//fname,',',
     .                    nf_strerror(status)
          stop
        endif
c
        call rdotr3(ncid,'sealv',1,sealv_ow,sealv_oe,sealv_os,sealv_on)
        call rdotru(ncid,ub_ow2,ub_oe2,ub_os2,ub_on2,
     .                   vb_ow2,vb_oe2,vb_os2,vb_on2)
        call rdotr3(ncid,'fice',1,ficem_ow2,ficem_oe2,
     .                            ficem_os2,ficem_on2)
        call rdotr3(ncid,'hice',1,hicem_ow2,hicem_oe2,
     .                            hicem_os2,hicem_on2)
        call rdotr3(ncid,'hsnw',1,hsnwm_ow2,hsnwm_oe2,
     .                            hsnwm_os2,hsnwm_on2)
        call rdotr3(ncid,'iage',1,iagem_ow2,iagem_oe2,
     .                            iagem_os2,iagem_on2)
c
c --- - close netcdf file
        status=nf_close(ncid)
        if (status.ne.nf_noerr) then
          write(*,'(2a)') 'nf_close: ',nf_strerror(status)
          stop
        endif
c
c --- - reconstruct pressure interfaces and bottom pressure
c
c --- - western boundary
c$OMP PARALLEL DO PRIVATE(i,q,phiup,pup,philo,plo)
        do j=1,jj
          do ibz=1,nbz
            i=ibz
            if (ip(i,j).eq.1) then
              q=(sealv_ow(ibz,j)+depths(i,j)*1.e2)
     .          /(z_ow(ibz,j,1)+depths(i,j)*1.e2)
              phiup=sealv_ow(ibz,j)*g
              pup=0.
              do k=1,kk
                philo=(sealv_ow(ibz,j)
     .               +q*(z_ow(ibz,j,k+1)-z_ow(ibz,j,1)))*g
                plo=getlpi(temp_ow2(ibz,j,k),saln_ow2(ibz,j,k),
     .                     phiup,philo,pup)
                dp_ow2(ibz,j,k)=plo-pup
                phiup=philo
                pup=plo
              enddo
              if (ibz.le.2) then
                pb_ow2(ibz,j)=plo
              endif
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c
c --- - eastern boundary
c$OMP PARALLEL DO PRIVATE(i,q,phiup,pup,philo,plo)
        do j=1,jj
          do ibz=1,nbz
            i=ii-ibz
            if (ip(i,j).eq.1) then
              q=(sealv_oe(ibz,j)+depths(i,j)*1.e2)
     .          /(z_oe(ibz,j,1)+depths(i,j)*1.e2)
              phiup=sealv_oe(ibz,j)*g
              pup=0.
              do k=1,kk
                philo=(sealv_oe(ibz,j)
     .               +q*(z_oe(ibz,j,k+1)-z_oe(ibz,j,1)))*g
                plo=getlpi(temp_oe2(ibz,j,k),saln_oe2(ibz,j,k),
     .                     phiup,philo,pup)
                dp_oe2(ibz,j,k)=plo-pup
                phiup=philo
                pup=plo
              enddo
              if (ibz.le.2) then
                pb_oe2(ibz,j)=plo
              endif
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c
c --- - southern boundary
c$OMP PARALLEL DO PRIVATE(j,q,phiup,pup,philo,plo)
        do ibz=1,nbz
          j=ibz
          do i=1,ii
            if (ip(i,j).eq.1) then
              q=(sealv_os(i,ibz)+depths(i,j)*1.e2)
     .          /(z_os(i,ibz,1)+depths(i,j)*1.e2)
              phiup=sealv_os(i,ibz)*g
              pup=0.
              do k=1,kk
                philo=(sealv_os(i,ibz)
     .               +q*(z_os(i,ibz,k+1)-z_os(i,ibz,1)))*g
                plo=getlpi(temp_os2(i,ibz,k),saln_os2(i,ibz,k),
     .                     phiup,philo,pup)
                dp_os2(i,ibz,k)=plo-pup
                phiup=philo
                pup=plo
              enddo
              if (ibz.le.2) then
                pb_os2(i,ibz)=plo
              endif
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c
c --- - northern boundary
c$OMP PARALLEL DO PRIVATE(j,q,phiup,pup,philo,plo)
        do ibz=1,nbz
          j=jj-ibz
          do i=1,ii
            if (ip(i,j).eq.1) then
              q=(sealv_on(i,ibz)+depths(i,j)*1.e2)
     .          /(z_on(i,ibz,1)+depths(i,j)*1.e2)
              phiup=sealv_on(i,ibz)*g
              pup=0.
              do k=1,kk
                philo=(sealv_on(i,ibz)
     .               +q*(z_on(i,ibz,k+1)-z_on(i,ibz,1)))*g
                plo=getlpi(temp_on2(i,ibz,k),saln_on2(i,ibz,k),
     .                     phiup,philo,pup)
                dp_on2(i,ibz,k)=plo-pup
                phiup=philo
                pup=plo
              enddo
              if (ibz.le.2) then
                pb_on2(i,ibz)=plo
              endif
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c
c --- - get model time step for next 2D boundary update
        nst21=nst22
        if (otr2dmon) then
          nm=mod(nmonth,12)+1
          nst22=nint((time+.5*(nd_in_m(nmonth)+nd_in_m(nm)))
     .               *nstep_in_day)
        else
        endif
        write (*,*) 'nst21,nst22',nst21,nst22
c
c       write (10) dp_ow1,dp_oe1,dp_os1,dp_on1,
c    .             dp_ow2,dp_oe2,dp_os2,dp_on2
c       write (10) pb_ow1,pb_oe1,pb_os1,pb_on1,
c    .             pb_ow2,pb_oe2,pb_os2,pb_on2
c       write (10) ub_ow1,ub_oe1,ub_os1,ub_on1,
c    .             ub_ow2,ub_oe2,ub_os2,ub_on2
c       write (10) vb_ow1,vb_oe1,vb_os1,vb_on1,
c    .             vb_ow2,vb_oe2,vb_os2,vb_on2
c       write (10) ficem_ow1,ficem_oe1,ficem_os1,ficem_on1,
c    .             ficem_ow2,ficem_oe2,ficem_os2,ficem_on2
c       write (10) hicem_ow1,hicem_oe1,hicem_os1,hicem_on1,
c    .             hicem_ow2,hicem_oe2,hicem_os2,hicem_on2
c       write (10) hsnwm_ow1,hsnwm_oe1,hsnwm_os1,hsnwm_on1,
c    .             hsnwm_ow2,hsnwm_oe2,hsnwm_os2,hsnwm_on2
c       write (10) iagem_ow1,iagem_oe1,iagem_os1,iagem_on1,
c    .             iagem_ow2,iagem_oe2,iagem_os2,iagem_on2
c       close (10)
c
      endif
c
      return
      end
