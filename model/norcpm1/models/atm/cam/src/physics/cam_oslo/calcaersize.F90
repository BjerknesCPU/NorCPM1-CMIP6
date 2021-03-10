!     Ø Seland Calculates mean volume size and hygroscopic growth for use in 
!     dry deposition
      subroutine calcaersize( ncol, &
                 t, h2ommr, pmid, pdel,wetrad,wetrho,relhum)


      use shr_kind_mod,only: r8 => shr_kind_r8
      use ppgrid
      use wv_saturation, only: qsat_water
      use aerosoldef
      use physconst,     only: rhoh2o
      implicit none

      integer,  intent(in) :: ncol               ! number of columns


      real(r8), intent(in) :: t(pcols,pver)      ! layer temperatures (K)
      real(r8), intent(in) :: h2ommr(pcols,pver) ! layer specific humidity
      real(r8), intent(in) :: pmid(pcols,pver)   ! layer pressure (Pa)
      real(r8), intent(in) :: pdel(pcols,pver)  ! layer pressure thickness (Pa)

      real(r8), intent(out):: wetrad(pcols,pver,ncaer)  
                              ! wet radius of aerosol (m)
      real(r8), intent(out):: wetrho(pcols,pver,ncaer) ! wet aerosol density  
      real(r8), intent(out):: relhum(pcols,pver) ! Relative humidity  
!     local variables
      integer  :: i,k,m,irelh,mm
      real(r8) :: e
      integer  ::l ! species index
      real(r8) :: xrh(pcols,pver)
      real(r8) :: qs(pcols,pver)        ! saturation specific humidity
      real(r8) :: rmeanvol(pcols,pver) ! Mean radius with respect to volume 
      integer  :: irh1(pcols,pver),irh2(pcols,pver)
      integer  :: t_irh1,t_irh2
      real(r8) :: t_rh1,t_rh2,t_xrh,rr1,rr2
      real(r8) :: volaer(pcols,pver)

      e = 2.718281828_r8
       
      wetrad(:,:,:) =0._r8
      rmeanvol(:,:)=0._r8

      do k=1,pver
        do i=1,ncol
          qs(i,k)=qsat_water(t(i,k),pmid(i,k))
          xrh(i,k) = h2ommr(i,k)/qs(i,k)
          xrh(i,k) = max(xrh(i,k),0.0_r8)
          xrh(i,k) = min(xrh(i,k),1.0_r8)
          relhum(i,k)=xrh(i,k)
          xrh(i,k)=min(xrh(i,k),rhtab(10))                	
	end do
     end do
 

     do irelh=1,9
       do k=1,pver
	  do i=1,ncol
	    if(xrh(i,k).ge.rhtab(irelh).and. &
              xrh(i,k).le.rhtab(irelh+1)) then
              irh1(i,k)=irelh
              irh2(i,k)=irelh+1
            end if
	  end do
	end do
      end do
      do k=1,pver
	do i=1,ncol
	  t_irh1 = irh1(i,k)
          t_irh2 = irh2(i,k)
          t_rh1  = rhtab(t_irh1)
          t_rh2  = rhtab(t_irh2)
          t_xrh  = xrh(i,k)


	  do m=1,ncaer
            mm=m+ixae-1
!     if(species_class(m).eq.spec_class_aerosol) then
            rmeanvol(i,k)=effsize(mm)*(e**(1.5_r8*(log(sgpart(mm)))**2))
	    rr1=rdivr0(t_irh1,mm)
            rr2=rdivr0(t_irh2,mm)
	    wetrad(i,k,m)= (((t_rh2-t_xrh)*rr1+(t_xrh-t_rh1)*rr2)/ &
               (t_rh2-t_rh1))*rmeanvol(i,k)
            volaer(i,k)=(rmeanvol(i,k)/wetrad(i,k,m))**3._r8
            wetrho(i,k,m)=rhopart(mm)*volaer(i,k)+(1._r8-volaer(i,k))*rhoh2o 
!    endif


          end do
        end do	
      end do

      return
      end subroutine calcaersize


