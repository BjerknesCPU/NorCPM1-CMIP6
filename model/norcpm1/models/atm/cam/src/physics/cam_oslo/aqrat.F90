

module aqrat


use shr_kind_mod, only: r8 => shr_kind_r8

!
!
!  Purpose Module to calculate aqueous-phase reaction rates, and to tabulate 
! the effect of mixing of gases from interstitial air to cloud droplets.
! 

implicit none
private
!------------------------------Parameters-------------------------------

!     Physical parameters needed for wetphase reactions

  integer it0,im25,itm,ip40,i00,im10
  parameter(it0=273-100,im25=273-25,i00=273,itm=it0+150,  &
    im10=273-10)
!
  parameter(ip40=273+40)


  real (r8) :: xkch2(it0:itm),xkco3(it0:itm),fact(it0:itm), wt(it0:itm)

  real (r8) :: tabhet(30,30),tabeff(30,30)

public tabchem
public calchet

contains 

      subroutine tabchem
!-----------------------------------------------------------------------
!
! Establish a look-up table for effective aqueous phase reaction rate 
! for a interval of so2 and h2o2 concentrations.
!
!---------------------------Code history--------------------------------
!
! Original version:  O. Seland
!
!
!-----------------------------------------------------------------------

! Local variables

      real(r8) :: xkmin,ch2o2(30),co3,cso2(30)
      real(r8) :: t,f,td,tw,h,ppb,patm
      integer :: i,n,nn,nt,ns,nh
      real(r8) :: fm,kdif,kdifa,deltat,b,dt,k0
      real(r8) :: ko3,kh2	
      real(r8) :: so20,o30,h2o20
      real(r8) :: o3a,h2o2a,so2a
      real(r8) :: o3,h2o2,so2
      real(r8) :: o3ax,h2o2ax,so2ax
      real(r8) :: o3x,h2o2x,so2x
      real(r8) :: do3,dh2o2,dso2
      real(r8) :: xkap
      real(r8) :: aa	
      real(r8) :: ph,thion,r0,xkp
!    real(r8) :: sof,oxf	
!    real(r8) :: fact(it0:itm),
!    real(r8) :: w(it0:itm),
      real(r8) :: hs(it0:itm),                          &
         xs2h2o(it0:itm),xhso3(it0:itm),xso3(it0:itm), &
         ho3(it0:itm),hs1(it0:itm),hs2(it0:itm),       &
         hso2(it0:itm),hh2o2(it0:itm)                 
      real(r8) :: k1o3,k2o3,k3o3,kh2o2,roww,x,xtest

      real(r8) :: e298, dhr, tt, ek
      ek(e298, dhr, tt) = e298*exp(dhr*(1._r8/tt - 1._r8/298._r8))      


!    h=3600._r8
      r0 = 0.082_r8

      xkp = 13._r8
      roww = 1000._r8
!    w=0.5e-3_r8
      f = 0.15_r8
      ph = 4._r8
      thion = 10._r8**(-ph)

      ppb = 1.e-9_r8

!    patm=0.9_r8
!    oxf=ppb*patm
!    sof=ppb*patm/1.42_r8


      dt=10._r8
      deltat=1800._r8
      b=(1._r8/3600._r8)*(sqrt(0.1_r8)/(1._r8-sqrt(0.1_r8)))

      kdif=b*(1._r8-sqrt(f))/sqrt(f)
      fm=1._r8-f
!    write(6,*) ' 1-f=',fm
      kdifa=b*(1._r8-sqrt(fm))/sqrt(fm)


      nn=180	
!    cso2(1)=0.05_r8*sof
!    cso2(2)=0.2_r8*sof
!    cso2(3)=0.6_r8*sof
!    cso2(4)=1._r8*sof
!    cso2(5)=1.5_r8*sof
!    cso2(6)=3._r8*sof
!    cso2(7)=5._r8*sof
!    cso2(8)=7._r8*sof
      do i=1,30
	cso2(i)=((10._r8**(real(i,r8)/10._r8))/100._r8)*ppb
	ch2o2(i)=((10._r8**(real(i,r8)/10._r8))/100._r8)*ppb
	xtest=10._r8*log10((cso2(i)*100._r8)/ppb)
!	write(6,*) i,cso2(i),cso2(i)
      end do	
      co3=35._r8*ppb
      do i=it0,im25
         fact(i)=0.001_r8
      enddo
      do i=im25+1,i00-1
         x=real(i-im25,r8)/real(i00-im25,r8)
         fact(i) = x
      enddo
      do i=i00,itm
         fact(i)=1._r8
      enddo

!      
!..w

      do i=it0,im10
         wt(i)=0.2e-3_r8
      enddo
      do i=im10,itm
         x=real(i-im10,r8)/real(ip40-im10,r8)
         wt(i) = 0.2e-3_r8 + 0.6e-3_r8*x*x
      end do   
!    do i=1,8
!       tabs2(i)=cso2(i)
!	write(6,*) i,tabs2(i)
!    end do
      xkmin=1._r8/(3600._r8*50._r8)
      k0=xkmin
!   
!..xkcmax:
!  xkcmax=1._r8/3600._r8
!   
!..o3 
!    k1o3=2.4e4_r8
!    k2o3=3.7e5_r8
!    k3o3=1.5e9_r8

!   

!   

      do i=it0,itm
        t=real(i,r8)
        nt=nint(t)

        hh2o2(nt) = ek(7.4e4_r8, 7302._r8, t)
        kh2o2 = ek(7.5e7_r8,-4430._r8,t)
        hs1(nt) = ek(1.3e-2_r8,1960._r8,t)
        hs2(nt) = ek(6.6e-8_r8,1500._r8,t)
        ho3(nt) = ek(1.13e-2_r8,2538._r8,t)
        k1o3 = 2.4e4_r8
        k2o3 = ek(3.7e5_r8,-5530._r8,t)
        k3o3 = ek(1.5e9_r8,-5280._r8,t)
        hso2(nt) = ek(1.23_r8,3147._r8, t) 
!c
         hs(nt)=hso2(nt)*r0*t*(1._r8 + hs1(nt)/thion &
          + hs1(nt)*hs2(nt)/(thion*thion))
         xs2h2o(nt)=1._r8/(1._r8 + hs1(nt)/thion        &
          + hs1(nt)*hs2(nt)/(thion*thion))
         xhso3(nt)=1._r8/(1._r8 + thion/hs1(nt) + hs2(nt)/thion)
         xso3(nt)=1._r8/(1._r8 + thion/hs2(nt)          &
         + thion*thion/(hs1(nt)*hs2(nt)))



      xkco3(nt)=fact(nt)*(wt(nt)/roww)*hs(nt)*ho3(nt)*    &
          (k1o3*xs2h2o(nt) + k2o3*xhso3(nt) + k3o3*xso3(nt))
!   
!   
      xkch2(nt)=fact(nt)*(wt(nt)/roww)*hs(nt)                    &
          *thion*hh2o2(nt)*kh2o2                                & 
          *xhso3(nt)/((1._r8+xkp*thion)*(1._r8+(wt(nt)/roww)*hh2o2(nt)))
!   
      end do     
!  Use t=283 K in the table

      t=283._r8
      nt=nint(t)

      ko3=xkco3(nt)
      kh2=xkch2(nt)
!   
      o30=co3

      do ns=1,30
        so20=cso2(ns)
        do nh=1,30




          h2o20=ch2o2(nh)
!   


!   
          o3a=o30
          h2o2a=h2o20
          so2a=so20
!   
          o3=o30
          h2o2=h2o20
          so2=so20
!   
!  write(6,*) ' k0=',k0,' ko3=',ko3,' kh2=',kh2
!  write(6,*) ' o3=',o3/oxf,' h2o2=',h2o2/oxf
!   
          do n=1,nn
!   
!..so2:
            xkap = k0 + ko3*o3 + kh2*h2o2 + kdif
            aa = kdif*so2a
            so2x = aa/xkap + (so2 - aa/xkap)*exp(-xkap*dt)
!   
!..so2a:
           xkap = kdifa
           aa = kdifa*so2
           so2ax = aa/xkap + (so2a - aa/xkap)*exp(-xkap*dt)
!   
!..o3:
           xkap = ko3*so2 + kdif
           aa = kdif*o3a
           o3x = aa/xkap + (o3 - aa/xkap)*exp(-xkap*dt)
!   
!..o3a:
           xkap =  kdifa
           aa = kdifa*o3
           o3ax = aa/xkap + (o3a - aa/xkap)*exp(-xkap*dt)
!   
!..h2o2:
           xkap = kh2*so2 + kdif
           aa = kdif*h2o2a
           h2o2x = aa/xkap + (h2o2 - aa/xkap)*exp(-xkap*dt)
!   
!..h2o2a:
           xkap = kdifa
           aa = kdifa*h2o2
           h2o2ax = aa/xkap + (h2o2a - aa/xkap)*exp(-xkap*dt)
!   
           so2=so2x
           so2a=so2ax
           o3=o3x
           o3a=o3ax
           h2o2=h2o2x
           h2o2a=h2o2ax
!   
         enddo
!   
!   
         dso2=so2-so20
         do3=o3-o30
         dh2o2=h2o2-h2o20
!   
!   
!   
!   
         tabeff(nh,ns) = -dso2/(deltat*so20)
!   

!   

!          k00(nh,ns)=k0+ko3*o30+kh2*h2o20
         tabhet(nh,ns)=k0+ko3*o30+kh2*h2o20
       enddo   
      enddo   
      return

      end subroutine tabchem


      subroutine calchet(hso2,hs1,hs2,ho3,hh2o2,    &
          xkh2o2,xk1o3,xk2o3,xk3o3,xkhet,xcldv,    &
          temp,twat,cso2,ch2o2,co3,ilw,frh2o2)	

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid

  implicit none


    real(r8), intent(in) :: hso2(pcols),hs1(pcols),hs2(pcols),ho3(pcols),hh2o2(pcols)
    real(r8), intent(in) :: xkh2o2(pcols),xk1o3(pcols),xk2o3(pcols),xk3o3(pcols)
    real(r8), intent(in) :: temp(pcols),twat(pcols)
    real(r8), intent(in) :: co3(pcols),cso2(pcols)
    real(r8), intent(in) :: xcldv(pcols)
    real(r8), intent(inout) :: xkhet(pcols)
    integer, intent(in) :: ilw
    real(r8), intent(out) :: frh2o2(pcols)
    real(r8), intent(inout) :: ch2o2(pcols)
! Local      
    integer :: i2,n,nn,nt(pcols)
    real(r8) :: r0,roww,xkp,ppb,xkmin,b,dt,deltat,k0
    real(r8) :: f(pcols),fm(pcols),kdif(pcols),kdifa(pcols)
    real(r8) :: hs(pcols),xs2h2o(pcols),xhso3(pcols),xso3(pcols)
    real(r8) :: kco3(pcols),kch2(pcols),thion(pcols)
    real(r8) :: ko3(pcols),kh2(pcols),o30(pcols),so20(pcols),h2o20(pcols)
    real(r8) :: xkap(pcols),aa(pcols),o3(pcols),h2o2(pcols),so2(pcols)
    real(r8) :: o3a(pcols),h2o2a(pcols),so2a(pcols)
    real(r8) :: so2x(pcols),o3x(pcols),h2o2x(pcols)
    real(r8) :: so2ax(pcols),o3ax(pcols),h2o2ax(pcols)
    real(r8) :: dso2(pcols),do3(pcols),dh2o2(pcols)
    real(r8) :: keff(pcols),k00(pcols),so2fac(pcols)
    real(r8) :: fac
    integer :: ncat(pcols),nco(pcols)
    real(r8) :: xup,xdown
    integer :: ierr(pcols)
    integer :: nstep
!      real(r8) :: frh2o2(pcols)  ! Fraction of h2o2 used up in the timestep
!      real(r8) :: redch2(pcols)  ! Reduction factor of ch2 due to no h2o2 left.
   r0=0.082_r8
      roww=1000._r8
     
!       fac=1._r8
!      f=0.15_r8
      xkp=13._r8
      ppb=1.e-9_r8
       do i2=1,ilw
       thion(i2)=10._r8**(-4._r8)
!      patm(i2)=tprs(i2)/101325
!      sof(i2)=ppb*patm(i2)/1.42_r8


!         hs(nt)=hso2(nt)*r0*t*(1._r8 + hs1(nt)/thion 
!     &     + hs1(nt)*hs2(nt)/(thion*thion))
!         xs2h2o(nt)=1._r8/(1._r8 + hs1(nt)/thion 
!     &     + hs1(nt)*hs2(nt)/(thion*thion))
!         xhso3(nt)=1._r8/(1._r8 + thion/hs1(nt) + hs2(nt)/thion)
!         xso3(nt)=1._r8/(1._r8 + thion/hs2(nt) 
!     &    + thion*thion/(hs1(nt)*hs2(nt)))
         ch2o2(i2)=max(ch2o2(i2),1.e-14_r8)
         hs(i2)=hso2(i2)*r0*temp(i2)*(1._r8 + hs1(i2)/thion(i2)  & 
          + hs1(i2)*hs2(i2)/(thion(i2)*thion(i2)))
         xs2h2o(i2)=1._r8/(1._r8 + hs1(i2)/thion(i2) &
          + hs1(i2)*hs2(i2)/(thion(i2)*thion(i2)))
         xhso3(i2)=1._r8/(1._r8 + thion(i2)/hs1(i2) + hs2(i2)/thion(i2)) 
         xso3(i2)=1._r8/(1._r8 + thion(i2)/hs2(i2)              &
         + thion(i2)*thion(i2)/(hs1(i2)*hs2(i2)))
	  nt(i2)=nint(temp(i2))
!	  kco3(i2)=xkco3(nt(i2))
!          kch2(i2)=xkch2(nt(i2))
!	write(6,*) 'kch1 ',i2,nt(i2),kch2(i2),kco3(i2)
!         twat(i2)=w(nt(i2))
         kco3(i2)=(twat(i2)/roww)*hs(i2)*ho3(i2)*  &
          (xk1o3(i2)*xs2h2o(i2) + xk2o3(i2)*xhso3(i2) & 
        + xk3o3(i2)*xso3(i2))

	   kch2(i2)=0.8_r8*(twat(i2)/roww)*hs(i2)*  &
          hh2o2(i2)*thion(i2)*xkh2o2(i2)*       &
          xhso3(i2)/(1._r8+xkp*thion(i2))   
           kch2(i2)=kch2(i2)/(1._r8+0.8_r8*(twat(i2)/roww)*hh2o2(i2))

!           redch2(i2)=ch2o2(i2)/(mdt*min(kch2(i2),2.8e-4_r8)*max(cso2(i2),1.e-18_r8))
!           frh2o2(i2)=1._r8/redch2(i2)
!           frh2o2(i2)=min(frh2o2(i2),1._r8)
!           redch2(i2)=min(redch2(i2),1._r8)
!
!           kch2(i2)=redch2(i2)*kch2(i2)

!feil         kch2(i2)=fact(nt(i2))*twat(i2)*hs(i2)
!     -     *hh2o2(i2)*thion(i2)*xkh2o2(i2)
!     -     *xhso3(i2)/((1._r8+xkp*thion(i2))
!feil     -     *(1._r8+twat(i2)*hh2o2(i2)))
      end do
!            write(6,*) 'etter kch2 ',kch2(i2)
!	if (nstep.gt.454) then
!	write(6,*) nstep,ilw,' etter hs '
!        end if

     
!c..xkmin (supposed to take into account always
!c..   existing Mn and Fe in polluted air):
     

      xkmin=1._r8/(3600._r8*50._r8)

      b=(1._r8/3600._r8)*(sqrt(0.1_r8)/(1._r8-sqrt(0.1_r8)))

      dt=10._r8
      deltat=1800._r8
      nn=180

      k0=xkmin     

!      write(6,*) ' f=',f




      do i2=1,ilw

!      f(i2)=xcldv(i2)
!      kdif(i2)=b*(1._r8-sqrt(f(i2)))/sqrt(f(i2))
!      fm(i2)=1._r8-f(i2)
! !     fm(i2)=max(fm(i2),0.02_r8)
!      write(6,*) ' 1-f=',fm
!      kdifa(i2)=b*(1._r8-sqrt(fm(i2)))/sqrt(fm(i2))

!	if (nstep.gt.454) then
!	write(6,*) nstep,ilw,' etter kdifa '
!        end if

        ko3(i2)=kco3(i2)
        kh2(i2)=kch2(i2)
!     
!        o30(i2)=co3(i2)

!        so20(i2)=cso2(i2)
!        so20(i2)=max(so20(i2),1.e-20_r8)
!        h2o20(i2)=ch2o2(i2)	



!        o3a(i2)=o30(i2)
!        h2o2a(i2)=h2o20(i2)
!        so2a(i2)=so20(i2)
!     
!        o3(i2)=o30(i2)
!        h2o2(i2)=h2o20(i2)
!        so2(i2)=so20(i2)

!        cso2(i2)=max(cso2(i2),1.e-8_r8)
!	cso2(i2)=max(cso2(i2),1.e-20_r8)
!	xtest=10._r8*log10((cso2(i)*100._r8)/sof)
!	write(6,*) i2,sof(i2),cso2(i2),patm(i2),ppb
	ncat(i2)=nint(10._r8*log10((max(cso2(i2),1.e-20_r8)*100._r8)/ppb)) 
	ncat(i2)=max(ncat(i2),1)     	
        ncat(i2)=min(ncat(i2),30)
!        write(6,*) 'foer k00'
        k00(i2)=co3(i2)*kco3(i2)+ch2o2(i2)*kch2(i2)+ & 
              fact(nt(i2))*xkmin
!        write(6,*) 'etter k00 ',k00(i2)

        ierr(i2)=0 
        nco(i2)=0
!	write(6,*) 'foer while ',i2,nco(i2)
        do while(ierr(i2).eq.0)

        nco(i2)=nco(i2)+1
        if (k00(i2).le.tabhet(nco(i2),ncat(i2))) ierr(i2)=1
        if (nco(i2).ge.30) ierr(i2)=1
        end do
        if (nco(i2).gt.1) then         
        xup=tabhet(nco(i2),ncat(i2))-k00(i2)
        xdown=k00(i2)-tabhet(nco(i2)-1,ncat(i2))
!
          if ((xup-xdown).gt.0) then
            xkhet(i2)=tabeff(nco(i2)-1,ncat(i2))
          else   
            xkhet(i2)=tabeff(nco(i2),ncat(i2))
          end if
        else
         xkhet(i2)=tabeff(nco(i2),ncat(i2))
       end if 
!!ostest
!       xkhet(i2)=k00(i2)
!!ostest
       xkhet(i2)=min(xkhet(i2),2.8e-4_r8)
       xkhet(i2)=min(xkhet(i2),k00(i2))
!       redch2(i2)=(xkhet(i2)/min(k00(i2),2.8e-4_r8)
       frh2o2(i2)=((ch2o2(i2)*kch2(i2))/k00(i2))
!       if(frh2o2(i2).gt.1) &
!       write(6,*) 'calchet ',frh2o2(i2),redch2(i2)
!        write(6,*) i2,redch2(i2),frh2o2(i2),xkhet(i2)
!	if (ch2o2(i2).lt.1.e-10_r8.and.cso2(i2).gt.1.e-10_r8) then
!       write(6,*) 'calchet',i2,temp(i2),nt(i2)
!	write(6,*) 'calchet ',i2,nt(i2),cso2(i2),ch2o2(i2),xkhet(i2)
!        write(6,*)  k00(i2),ncat(i2),nco(i2)
!        write(6,*)  twat(i2),tabhet(nco(i2),ncat(i2)), &
!        tabeff(nco(i2),ncat(i2))
!	end if
       enddo






      return 
      end subroutine calchet


end module aqrat
