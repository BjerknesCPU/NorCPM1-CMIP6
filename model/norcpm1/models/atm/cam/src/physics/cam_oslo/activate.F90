      subroutine activate_ny(top,wbar,sigw,wdiab,wminf,wmaxf,tair,rhoair,  &
                          na,ntype,nmode,ma,lnsigman,         &
                          nu,phi,mwaer,rhodry,solubl,fn,fluxn,smax)

!      calculates number, surface, and mass fraction of aerosols activated as CCN
!      calculates flux of cloud droplets, surface area, and aerosol mass into cloud
!      assumes an internal mixture within each of up to pmode multiple aerosol modes
!      a gaussiam spectrum of updrafts can be treated.

!      cgs units

!      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
!      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

      use shr_kind_mod,	only: r8 => shr_kind_r8,r4 => shr_kind_r4
      use physconst, only: rair, epsilo, cpair, rh2o, latvap, gravit,   &
                                 rhoh2o, mwh2o, r_universal

      implicit none
	integer ptype, pmode
	parameter(pmode=12,ptype=5)
!      input

      logical, intent(in) :: top        ! true if cloud top, false if cloud base or new cloud
      real(r8), intent(in) :: wbar
      real(r8), intent(in) :: sigw
      real(r8), intent(in) :: wdiab         ! diabatic vertical velocity (0 if adiabatic)
      real(r8), intent(in) :: wminf         ! minimum updraft velocity for integration (cm/s)
      real(r8), intent(in) :: wmaxf         ! maximum updraft velocity for integration (cm/s)
      real(r8), intent(in) :: tair          ! air temperature (K)
      real(r8) tairc
      real(r8), intent(in) :: rhoair        ! air density (g/cm3)
      real(r8), intent(in) :: na(pmode)           ! aerosol number concentration (/cc)
      integer, intent(in) :: ntype(pmode)      ! number of aerosol types
      integer, intent(in) :: nmode      ! number of aerosol modes
      real(r8), intent(in) :: ma(ptype,pmode)     ! aerosol mass concentration (g/cc)
      real(r8), intent(in) :: lnsigman(pmode)  ! ln of geometric standard deviation of aerosol size distribution
      real(r8), intent(in) :: nu(ptype,pmode)     ! ions/molecule
      real(r8), intent(in) :: phi(ptype,pmode)    ! osmotic coefficient
      real(r8), intent(in) :: mwaer(ptype,pmode)  ! molecular weight of salt (grams/mole)
      real(r8), intent(in) :: solubl(ptype,pmode) ! mass fraction of soluble material in dry
!                                                 ! aerosol
      real(r8), intent(in) :: rhodry(ptype,pmode) ! density of dry aerosol (g/cm3)
!      output
      real(r8), intent(out) :: fn(pmode)      ! number fraction of aerosols activated
      real(r8), intent(out) :: fluxn(pmode)   ! flux of activated aerosol number fraction into cloud (cm/s)
      real(r8), intent(out) :: smax
      real(r8) fs(pmode)      ! surface fraction of aerosols activated
      real(r8) fm(pmode)      ! mass fraction of aerosols activated
      real(r8) fluxs(pmode)   ! flux of activated aerosol surface fraction (cm/s)
      real(r8) fluxm(pmode)   ! flux of activated aerosol mass fraction into cloud (cm/s)
!      local

      real(r8) derf,derfc

      integer nx
      parameter (nx=200)
      integer iquasisect_option, isectional
      real(r8) integ,integf
      real(r8) surften       ! surface tension of water w/respect to air (dyne/cm)
      data surften/76._r8/
      save surften
      real(r8) p0     ! reference pressure (dyne/cm2)
      real(r8) t0     ! reference temperature
      data p0/1013.25e3_r8/,t0/273.15_r8/
      save p0,t0
      real(r8) xmin(pmode),xmax(pmode) ! ln(r) at section interfaces
      real(r8) surfmin(pmode),surfmax(pmode) ! surface area at interfaces
      real(r8) volmin(pmode),volmax(pmode) ! volume at interfaces
      real(r8) surfc(pmode) ! surface concentration (cm2/cm3)
      real(r8) vol(pmode) ! total aerosol volume  (cm3)
      real(r8) volc(pmode) ! total aerosol volume  concentration (cm3/cm3)
      real(r8) solfrac ! mass fraction of soluble material in mixed aerosol
      real(r8) tmass ! total aerosol mass concentration (g/cm3)
      real(r8) sumphi ! mean osmotic coefficient of aerosol
      real(r8) sign(pmode)    ! geometric standard deviation of size distribution
      real(r8) alnsign(pmode) ! natl log of geometric standard dev of aerosol
      real(r8) am(pmode) ! number mode radius of dry aerosol (cm)
      real(r8) hygro(pmode) ! hygrscopicity
      real(r8) lnhygro(pmode) ! ln(b)
      real(r8) pres ! pressure (dynes/cm2)
      real(r8) diff0,conduct0
      real(r8) qs ! water vapor saturation mixing ratio
      real(r8) esatcgs
      real(r8) dqsdt ! change in qs with temperature
      real(r8) g ! thermodynamic function (cm2/s)
      real(r8) sm(pmode) ! critical supersaturation for number mode radius
      real(r8) zeta, eta(pmode)
      real(r8) lnsm(pmode) ! ln(sm)
      real(r8) lnsmax ! ln(smax)
      real(r8) alpha
      real(r8) gamma
      real(r8) beta
      real(r8) sqrtg
      real(r8) asub(pmode),bsub(pmode) ! coefficients of submode size distribution N=a+bx
      real(r8) totn ! total aerosol number concentration
      real(r8) aten ! surface tension parameter
      real(r8) gmrad ! geometric mean radius
      real(r8) gmradsq ! geometric mean of radius squared
      real(r8) gmlnsig ! geometric standard deviation
      real(r8) gmsm ! critical supersaturation at radius gmrad
      real(r8) sumflxn(pmode)
      real(r8) sumfn(pmode)
      real(r8) sumns
      real(r8) fnold(pmode)   ! number fraction activated
      real(r8) wold,gold
      real(r8) alog2,alog3,alogaten
      real(r8) alogam
      real(r8) ccc
      real(r8) wmin,wmax,w,dw,dwmax,dwmin,wnuc,dwnew,wb
      real(r8) dfmin,dfmax,fnew,fold,fnmin,fnbar
      real(r8) alw,sqrtalw
      real(r8) x,xmincoeff,xcut,volcut,surfcut,arg
      real(r8) z,z1,z2,wf1,wf2,zf1,zf2,gf1,gf2,gf
      integer m,n,k
!      numerical integration parameters
      real(r8) eps,fmax,sds
      data eps/0.3_r8/,fmax/0.99_r8/,sds/3._r8/
!      mathematical constants
      real(r8) third, twothird, sixth,zero,one,two,three
      data third/0.333333_r8/, twothird/0.66666667_r8/, sixth/0.166666667_r8/,zero/0._r8/,one/1._r8/,two/2._r8/,three/3._r8/
      save third, sixth,twothird,zero,one,two,three
      real(r8) sq2, sqpi, pi
      data sq2/1.4142136_r8/, sqpi/1.7724539_r8/,pi/3.14159_r8/
      save sq2,sqpi,pi

      real(r8) raircgs, cpaircgs, rh2ocgs, latvapcgs, gravitcgs, rhoh2ocgs,  &
             mwh2ocgs, r_universalcgs
      save raircgs, cpaircgs, rh2ocgs, latvapcgs, gravitcgs, rhoh2ocgs,  &
           mwh2ocgs, r_universalcgs

      real(r8) naa     ! aerosol number concentration (/cc)
      real(r8) namin   ! minimum aerosol number concentration (/cc)
      data namin/1._r8/

      integer ndist(nx)  ! accumulates frequency distribution of integration bins required
      data ndist/nx*0/
      save eps,fmax,sds,namin,ndist

      smax=0._r8
      wmin=0._r8
      wmax=0._r8
      tairc=0._r8
      integ=0._r8
      integf=0._r8
      solfrac=0._r8
      tmass=0._r8
      sumphi=0._r8 
      pres=0._r8   
      diff0=0._r8  
      conduct0=0._r8
      qs=0._r8      
      dqsdt=0._r8     
      esatcgs=0._r8   
      g=0._r8         
      zeta=0._r8      
      lnsmax=0._r8    
      alpha=0._r8     
      beta=0._r8      
      gamma=0._r8     
      sqrtg=0._r8     
      totn=0._r8      
      aten=0._r8      
      gmrad=0._r8     
      gmradsq=0._r8   
      gmlnsig=0._r8   
      gmsm=0._r8      
      sumns=0._r8     
      w=0._r8        
      wold=0._r8     
      gold=0._r8     
      alogam=0._r8   
      dw=0._r8       
      dwmin=0._r8    
      dwmax=0._r8    
      dwnew=0._r8    
      wnuc=0._r8     
      wb=0._r8       
      dfmin=0._r8    
      dfmax=0._r8    
      fnew=0._r8     
      fold=0._r8     
      fnmin=0._r8   
      fnbar=0._r8    
      alw=0._r8      
      sqrtalw=0._r8  
      x=0._r8        
      xmincoeff=0._r8
      xcut=0._r8     
      volcut=0._r8   
      surfcut=0._r8  
      arg=0._r8      
      z=0._r8        
      z1=0._r8       
      z2=0._r8       
      zf1=0._r8      
      zf2=0._r8      
      wf1=0._r8      
      wf2=0._r8      
      gf1=0._r8      
      gf2=0._r8      
      gf=0._r8       
      naa=0._r8      
!      do m=1,nmode
      do m=1,pmode    
         fn(m)=0._r8
         sumflxn(m)=0._r8
         sumfn(m)=0._r8
	 fnold(m)=0._r8
         fs(m)=0._r8
         fm(m)=0._r8
         fluxn(m)=0._r8
         fluxs(m)=0._r8
         fluxm(m)=0._r8
         sign(m)=0._r8   
         alnsign(m)=0._r8
         xmin(m)=0._r8   
         xmax(m)=0._r8   
         surfmin(m)=0._r8
         surfmax(m)=0._r8
         volmin(m)=0._r8 
         volmax(m)=0._r8 
         surfc(m)=0._r8  
         vol(m)=0._r8    
         volc(m)=0._r8   
         am(m)=0._r8     
         hygro(m)=0._r8  
         lnhygro(m)=0._r8
         sm(m)=0._r8     
         lnsm(m)=0._r8   
         eta(m)=0._r8    
         asub(m)=0._r8   
         bsub(m)=0._r8   
      end do 

!     for nmode>9, a sectional approach is used and isectional = iquasisect_option
!     activation fractions (fn,fs,fm) are computed as follows
!     iquasisect_option = 1,3 - each section treated as a narrow lognormal
!     iquasisect_option = 2,4 - within-section dn/dx = a + b*x,  x = ln(r)
!     smax is computed as follows (when explicit activation is OFF)
!     iquasisect_option = 1,2 - razzak-ghan modal parameterization with
!     single mode having same ntot, dgnum, sigmag as the combined sections
!     iquasisect_option = 3,4 - razzak-ghan sectional parameterization
!     for nmode=<9, a modal approach is used and isectional = 0
#undef ACTIVATE_QUASISECT_OPT
#ifdef ACTIVATE_QUASISECT_OPT
      iquasisect_option = ACTIVATE_QUASISECT_OPT
      if ((iquasisect_option .lt. 1) .or.(iquasisect_option .gt. 4)) then
         print *,'stopping in activate'
         print *,'ACTIVATE_QUASISECT_OPT must be 1, 2, 3, or 4'
      endif
#else
      iquasisect_option = 2
#endif
      isectional = 0
      if (nmode .ge. 13) isectional = iquasisect_option
      if (isectional .gt. 0) then
             print *,'stopping in activate'
             stop
      endif
      if(nmode.eq.1.and.na(1).lt.1.e-20_r8)then
         fn(1)=0._r8
         fs(1)=0._r8
         fm(1)=0._r8
         fluxn(1)=0._r8
         fluxs(1)=0._r8
         fluxm(1)=0._r8
         return
      endif
      r_universalcgs=r_universal*1.e4_r8
      raircgs=rair*1.e4_r8
      cpaircgs=cpair*1.e4_r8
      rh2ocgs=rh2o*1.e4_r8
      latvapcgs=latvap*1.e4_r8
      gravitcgs=gravit*100._r8
      rhoh2ocgs=rhoh2o*1.e-3_r8
      mwh2ocgs=mwh2o

      pres=raircgs*rhoair*tair
      tairc = tair - 273.0_r8
      diff0=0.211_r8*(p0/pres)*(tair/t0)**1.94_r8
      conduct0=(5.69_r8+0.017_r8*(tair-t0))*4.186e2_r8
      esatcgs = 611.2_r8*exp((17.67_r8*tairc)/(tairc+243.5_r8))
      qs=epsilo*esatcgs/pres
      dqsdt=latvapcgs/(rh2ocgs*tair*tair)*qs
      alpha=gravitcgs*(latvapcgs/(cpaircgs*rh2ocgs*tair*tair)-1._r8/(raircgs*tair))
      gamma=(1._r8+latvapcgs/cpaircgs*dqsdt)/(rhoair*qs)
      g=1._r8/(rhoh2ocgs/(diff0*rhoair*qs)                                    &
          +latvapcgs*rhoh2ocgs/(conduct0*tair)*(latvapcgs/(rh2ocgs*tair)-1._r8))
      sqrtg=sqrt(g)
      beta=4*pi*rhoh2ocgs*g*gamma
      totn=0
      gmrad=0._r8
      gmradsq=0._r8
      sumns=0._r8
      aten=2._r8*mwh2ocgs*surften/(r_universalcgs*tair*rhoh2ocgs)
      alogaten=log(aten)
      alog2=log(two)
      alog3=log(three)
      ccc=4._r8*pi*third

      do m=1,nmode
!cak_inputerror_4apr08         alnsign(m)=log(sigman(m))
         alnsign(m)=lnsigman(m)
!         internal mixture of aerosols
         volc(m)=0._r8
         solfrac=0._r8
         tmass=0._r8

         sumphi=0._r8
         do n=1,ntype(m)
            volc(m)=volc(m)+ma(n,m)/(rhodry(n,m))
            solfrac=solfrac+ma(n,m)*solubl(n,m)
            tmass=tmass+ma(n,m)
            sumphi=sumphi+nu(n,m)*ma(n,m)*phi(n,m)*solubl(n,m)/mwaer(n,m)
         enddo
!         convert concentrations to per particle
         vol(m)=volc(m)/na(m)
!    modified the following line by Yang Zhang, 1/26/99
!    because it gives sm(m)= 1.0, and causes overflow
!    in later calculation
         if(tmass.gt.1.e-30_r8.and.na(m).gt.1.e-30_r8)then
            solfrac=solfrac/tmass
            hygro(m)= mwh2ocgs*sumphi/(volc(m)*rhoh2ocgs)
            lnhygro(m)=log(hygro(m))
!            number mode radius (cm)
            am(m)=exp(-1.5_r8*alnsign(m)*alnsign(m))*              &
              (3._r8*volc(m)/(4._r8*pi*na(m)))**third
            if (isectional .gt. 0) then
!               sectional model.
!               need to use bulk properties because parameterization doesn't
!               work well for narrow bins.
!#ifdef CHEM
!  for explicit activation, am(m) is not used
!  for quasi-sect activation, when is am(m)=0.5_r8*dcen_sect(m) appropriate ???
!                am(m)=0.5_r8*dcen_sect(m)
!#endif
               totn=totn+na(m)
               alogam=log(am(m))
               gmrad=gmrad+na(m)*alogam
               gmradsq=gmradsq+na(m)*alogam*alogam
            endif

            if(hygro(m).gt.1.e-10_r8)then
               sm(m)=2._r8*aten/(3._r8*am(m))*sqrt(aten/(3._r8*hygro(m)*am(m)))
            else
               sm(m)=100._r8
            endif
         else
            sm(m)=1._r8
         endif
         lnsm(m)=log(sm(m))
         if ((isectional .eq. 3) .or. (isectional .eq. 4)) then
            sumns=sumns+na(m)/sm(m)**twothird
         endif
      enddo

!      sjg 7-16-98  upward
      naa=max(na(1),namin)
      if(sigw.gt.1.e-3_r8)then ! spectrum of updrafts
         if(top)then
           wmax=0._r8
           wmin=min(zero,-wdiab)
         else
           wmax=min(wmaxf,wbar+sds*sigw)
           wmin=max(wminf,-wdiab)
         endif
         wmin=max(wmin,wbar-sds*sigw)
         w=wmin
         dwmax=eps*sigw
         dw=dwmax
         dfmax=0.2_r8
         dfmin=0.1_r8
         if(wmax.le.w)then
            do m=1,nmode
               fluxn(m)=0._r8
               fn(m)=0._r8
               fluxs(m)=0._r8
               fs(m)=0._r8
               fluxm(m)=0._r8
               fm(m)=0._r8
            enddo
            return
         endif
         do m=1,nmode
            sumflxn(m)=0._r8
            sumfn(m)=0._r8
	    fnold(m)=0._r8
         enddo

!         z=(w-wbar)/(sigw*sq2)
!         sumg=sigw*0.5_r8*sq2*sqpi*(1+derf(z))
         fold=0._r8
         wold=0._r8
         gold=0._r8

         dwmin = min(dwmax,0.01_r8)
     do n=1,200
 100        wnuc=w+wdiab


            alw=alpha*wnuc
            sqrtalw=sqrt(alw)
            zeta=2._r8*sqrtalw*aten/(3._r8*sqrtg)
            if (isectional .gt. 0) then      ! "if" is not in use, only "else"
!              sectional model.
!              use bulk properties

              if(totn.gt.1.e-10_r8)then
                 eta(1)=2*alw*sqrtalw/(totn*beta*sqrtg)
              else
                 eta(1)=1.e10_r8
              endif
              call activem(zeta,eta,1,gmsm,gmlnsig,smax)
              lnsmax=log(smax)
              x=2*(log(gmsm)-lnsmax)/(3*sq2*gmlnsig)
              fnew=0.5_r8*(1._r8-derf(x))

            else
              do m=1,nmode
                 if(na(m).gt.1.e-10_r8)then
! Konkurranse-effekten
		     !if(na(m).gt.500._r8)then
                    !eta(m)=2*alw*sqrtalw/(500._r8*beta*sqrtg)
		!write(*,*)'ETA 1 & 2', eta(m) ,2*alw*sqrtalw/(na(m)*beta*sqrtg)
		!else
		!write(*,*)'NA ', na(m)
		   eta(m)=2*alw*sqrtalw/(na(m)*beta*sqrtg)
		!end if
                 else
                    eta(m)=1.e10_r8
                 endif
              enddo

              call activem(zeta,eta,nmode,sm,alnsign,smax)

              lnsmax=log(smax)

              x=2*(lnsm(nmode)-lnsmax)/(3*sq2*alnsign(nmode))
              fnew=0.5_r8*(1._r8-derf(x))


            endif


            dwnew = dw
            if(fnew-fold.gt.dfmax.and.n.gt.1)then
!              reduce updraft increment for greater accuracy in integration
	       if (dw .gt. 1.01_r8*dwmin) then
                  dw=0.7_r8*dw
                  dw=max(dw,dwmin)
                  w=wold+dw
                  go to 100
               else
                  dwnew = dwmin
               endif
            endif

            if(fnew-fold.lt.dfmin)then
!              increase updraft increment to accelerate integration
	       dwnew=min(1.5_r8*dw,dwmax)
            endif
            fold=fnew

            z=(w-wbar)/(sigw*sq2)
            g=exp(-z*z)
            fnmin=1._r8
            xmincoeff=alogaten-2._r8*third*(lnsmax-alog2)-alog3

            do m=1,nmode
               if ((isectional .eq. 2) .or. (isectional .eq. 4)) then
!                 sectional
!                  within-section dn/dx = a + b
                  xcut=xmincoeff-third*lnhygro(m)
!                   if(lnsmax.lt.lnsmn(m))then
                   if(xcut.gt.xmax(m))then
                     fn(m)=0._r8
!                   elseif(lnsmax.gt.lnsmx(m))then
                  elseif(xcut.lt.xmin(m))then
                     fn(m)=1._r8
                  else
                     volcut=exp(3._r8*xcut)
                     surfcut=exp(2._r8*xcut)
                     fn(m)=(asub(m)*(xmax(m)-xcut)                     &
                         +0.5_r8*bsub(m)*(xmax(m)*xmax(m)-xcut*xcut))/na(m)
                     if(fn(m).lt.0._r8.or.fn(m).gt.1._r8)then
                        print *,'fn(',m,')=',fn(m),' in activate'
                        fn(m)=max(fn(m),zero)
                        fn(m)=min(fn(m),one)
                     endif
                  endif
               else
!                 modal
                  x=2*(lnsm(m)-lnsmax)/(3*sq2*alnsign(m))
                  fn(m)=0.5_r8*(1._r8-derf(x))
               endif
               fnmin=min(fn(m),fnmin)
!               integration is second order accurate
!               assumes linear variation of f*g with w
               wb=(w+wold)
               fnbar=(fn(m)*g+fnold(m)*gold)

               if((top.and.w.lt.0._r8).or.(.not.top.and.w.gt.0._r8))then
                  sumflxn(m)=sumflxn(m)+sixth*(wb*fnbar           &
                      +(fn(m)*g*w+fnold(m)*gold*wold))*dw
               endif
               sumfn(m)=sumfn(m)+0.5_r8*fnbar*dw
               fnold(m)=fn(m)
            enddo
!            sumg=umg+0.5_r8*(g+gold)*dw
            gold=g
            wold=w
            dw=dwnew
            if(n .gt. 1 .and.(w.gt.wmax.or.fnmin.gt.fmax))go to 20
            w=w+dw
         enddo
         print *,'do loop is too short in activate'
         print *,'wmin=',wmin,' w=',w,' wmax=',wmax,' dw=',dw
         print *,'wbar=',wbar,' sigw=',sigw,' wdiab=',wdiab
         print *,'wnuc=',wnuc
         print *,'na=',(na(m),m=1,nmode)
!   dump all subr parameters to allow testing with standalone code
!   (build a driver that will read input and call activate)
         print *,'top,wbar,sigw,wdiab,tair,rhoair,nmode='
         print *, top,wbar,sigw,wdiab,tair,rhoair,nmode
         print *,'na='
         print *, na
         print *,'ntype='
         print *, ntype
         print *,'ma='
         print *, ma
         print *,'lnsigman='
         print *, lnsigman
         print *,'nu='
         print *, nu
         print *,'phi='
         print *, phi
         print *,'rhodry='
         print *, rhodry
         print *,'mwaer='
         print *, mwaer
         print *,'solubl='
         print *, solubl

         stop
   20    continue
         ndist(n)=ndist(n)+1
         if(.not.top.and.w.lt.wmaxf)then

!            contribution from all updrafts stronger than wmax
!            assuming constant f (close to fmax)
            wnuc=w+wdiab

            z1=(w-wbar)/(sigw*sq2)
            z2=(wmaxf-wbar)/(sigw*sq2)
            g=exp(-z1*z1)
            integ=sigw*0.5_r8*sq2*sqpi*(derfc(z1)-derfc(z2))
!            consider only upward flow into cloud base when estimating flux
            wf1=max(w,zero)
            zf1=(wf1-wbar)/(sigw*sq2)
            gf1=exp(-zf1*zf1)
            wf2=max(wmaxf,zero)
            zf2=(wf2-wbar)/(sigw*sq2)
            gf2=exp(-zf2*zf2)
            gf=(gf1-gf2)
            integf=wbar*sigw*0.5_r8*sq2*sqpi*(derfc(zf1)-derfc(zf2))+sigw*sigw*gf

            do m=1,nmode
               sumflxn(m)=sumflxn(m)+integf*fn(m)
               sumfn(m)=sumfn(m)+fn(m)*integ
            enddo
!            sumg=sumg+integ
         endif


         do m=1,nmode
            fn(m)=sumfn(m)/(sq2*sqpi*sigw)
            if(fn(m).gt.1.01_r8)then
               print *,'fn=',fn(m),' > 1 in activate'
            endif
            fluxn(m)=sumflxn(m)/(sq2*sqpi*sigw)
         enddo

      else

!         uniform updraft

         wnuc=wbar+wdiab
         if(wnuc.gt.0._r8)then

            w=wbar
            alw=alpha*wnuc
            sqrtalw=sqrt(alw)
            zeta=2._r8*sqrtalw*aten/(3._r8*sqrtg)

            if (isectional .gt. 0) then
!              sectional model.
!              use bulk properties

              if(totn.gt.1.e-10_r8)then
                 eta(1)=2*alw*sqrtalw/(totn*beta*sqrtg)
              else
                 eta(1)=1.e10_r8
              endif
              call activem(zeta,eta,1,gmsm,gmlnsig,smax)

            else

              do m=1,nmode

                 if(na(m).gt.1.e-10_r8)then
                    eta(m)=2*alw*sqrtalw/(na(m)*beta*sqrtg)
                 else
                    eta(m)=1.e10_r8
                 endif
              enddo

              call activem(zeta,eta,nmode,sm,alnsign,smax)

            endif
	   
            lnsmax=log(smax)

            xmincoeff=alogaten-2._r8*third*(lnsmax-alog2)-alog3


            do m=1,nmode
               if ((isectional .eq. 2) .or. (isectional .eq. 4)) then
!                 sectional
!                  within-section dn/dx = a + b*x
                  xcut=xmincoeff-third*lnhygro(m)
!                   if(lnsmax.lt.lnsmn(m))then
                   if(xcut.gt.xmax(m))then
                     fn(m)=0._r8

!                   elseif(lnsmax.gt.lnsmx(m))then
                  elseif(xcut.lt.xmin(m))then
                     fn(m)=1._r8
                  else
                     volcut=exp(3._r8*xcut)
                     surfcut=exp(2._r8*xcut)
                     fn(m)=(asub(m)*(xmax(m)-xcut)                      &
                         +0.5_r8*bsub(m)*(xmax(m)*xmax(m)-xcut*xcut))/na(m)
                     if(fn(m).lt.0._r8.or.fn(m).gt.1._r8)then
                        print *,'fn(',m,')=',fn(m),' in activate'
                        fn(m)=max(fn(m),zero)
                        fn(m)=min(fn(m),one)
                     endif
                     fs(m)=(0.5_r8*asub(m)*(surfmax(m)-surfcut)           &
                         +0.5_r8*bsub(m)*(xmax(m)*surfmax(m)-xcut*surfcut &
                              -0.5_r8*(surfmax(m)-surfcut) ) ) /surfc(m)
                        fs(m)=max(fs(m),zero)
                        fs(m)=min(fs(m),one)
                     fm(m)=ccc*(third*asub(m)*(volmax(m)-volcut)           &
                           +third*bsub(m)*(xmax(m)*volmax(m)-xcut*volcut   &
                               -third*(volmax(m)-volcut) ) ) /volc(m)
                     if(fm(m).lt.0._r8.or.fm(m).gt.1._r8)then
                        print *,'fm(',m,')=',fm(m),' in activate'
                        fm(m)=max(fm(m),zero)
                        fm(m)=min(fm(m),1._r8)
                     endif
                  endif
               else
!                 modal
		  x=2*(lnsm(m)-lnsmax)/(3*sq2*alnsign(m))
                  fn(m)=0.5_r8*(1._r8-derf(x))
		  arg=x-sq2*alnsign(m)
                  fs(m)=0.5_r8*(1._r8-derf(arg))
		  arg=x-1.5_r8*sq2*alnsign(m)
                  fm(m)=0.5_r8*(1._r8-derf(arg))
               endif
                if((top.and.wbar.lt.0._r8).or.(.not.top.and.wbar.gt.0._r8))then
                   fluxn(m)=fn(m)*w
                else
                   fluxn(m)=0._r8
                endif
            enddo
         else
            do m=1,nmode
               fn(m)=0._r8
               fluxn(m)=0._r8
            enddo
         endif
      endif
      return
      end


!cak----------------------------------------------------------------------------

      subroutine activem(zeta,eta,nmode,sm,alnsign,smax)

!      calculates maximum supersaturation for multiple
!      competing aerosol modes.

!      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
!      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.
      use shr_kind_mod,	only: r8 => shr_kind_r8,r4 => shr_kind_r4

      implicit none

      integer pmode
      parameter (pmode=12)
      integer, intent(in) :: nmode ! number of modes
      real(r8), intent(in) :: sm(pmode) ! critical supersaturation for number mode radius
      real(r8), intent(in) :: zeta
      real(r8), intent(in) :: eta(pmode)
      real(r8), intent(in) :: alnsign(pmode) ! ln(sigma)
      real(r8) f1(pmode)
      real(r8), intent(out) :: smax ! maximum supersaturation
      real(r8) g1, g2
      save f1

!t      logical first
!t      data first/.true./
!t      save first
      real(r8) twothird,sum
      data twothird/0.666666666_r8/
      save twothird
      integer m  ! mode index

!cak+t
         smax=0._r8
         do m=1,pmode
            f1(m)=0._r8
         enddo
!cak-t

!t      if(first)then
!         calculate and save f1(sigma). assumes sigma is invariant.
         if(nmode.gt.pmode)then
            print *,'nmode=',nmode,' pmode=',pmode,' in activem'
            call exit(1)
         endif
         do m=1,nmode
            f1(m)=0.5_r8*exp(2.5_r8*alnsign(m)*alnsign(m))
         enddo
!t         first=.false.
!t      endif

      do m=1,nmode
         if(zeta.gt.1.e5_r8*eta(m).or.sm(m)*sm(m).gt.1.e5_r8*eta(m))then
!            weak forcing. essentially none activated
            smax=1.e-20_r8
         else
!            significant activation of this mode. calc activation all modes.
            go to 1 
         endif
      enddo

      return

  1   continue

      sum=0._r8
      do m=1,nmode
         if(eta(m).gt.1.e-20_r8)then
            g1=sqrt(zeta/eta(m))
            g1=g1*g1*g1
            g2=sm(m)/sqrt(eta(m)+3*zeta)
            g2=sqrt(g2)
            g2=g2*g2*g2
            sum=sum+(f1(m)*g1+(1.+0.25_r8*alnsign(m))*g2)/(sm(m)*sm(m))
         else
            sum=1.e20_r8
         endif
      enddo

      smax=1._r8/sqrt(sum)
 
     end 


