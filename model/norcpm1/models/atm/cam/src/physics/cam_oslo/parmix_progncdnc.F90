
      subroutine parmix_progncdnc(lchnk,ncol,nlons,ind,lev,Nnatk,Cnso4,Cnbc,Cnoc,Cabce,Ca, &
                        f_c,f_bc,f_aq,fnbc,faitbc,relh,Ccn,cxstot, &
                        n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14, &
!TS++
		Caitso4, Caitbc, Caitocbc, rhoocbc, Cas75, C_dst2, C_dst3, &
		C_ss1, C_ss2, C_ss3, cam, facm, fbcm, faqm, lnsig, Massratio, &
                logsig, version2, lthick, fso4coat)
!TS-- 
!     This subroutine finds CCN concentrations for a composite aerosol,
!     where parameters for each external mode is found by linear interpolation
!     in the tables ccnk1.out-ccnk14.out. Optimized (CCM3) April 2002 by
!     Arild Burud/NoSerC, and updated with 2 new modes (OC) and rewritten for
!     CAM2 by Trude Storelvmo and Alf Kirkevag, November 2004. Updated and 
!     rewritten for CAM3 (CAM-Oslo, with prognostic sea-salt and mineral
!     aerosols) by Alf Kirkevag in August 2005. Modified for new aerosol schemes 
!     by Alf Kirkevaag in January 2006.  

        use shr_kind_mod, only: r8 => shr_kind_r8
	use ppgrid,	only: pcols, pver, pverp
        use opttab,   only: nbmodes, nmodes, nbmp1
        use aerosoldef
        use physconst,    only: pi
	use const

      implicit none

      logical, intent(in)  :: version2          ! version2= .true. --> version 2 of lognormal fit calculations used
      integer, intent(in)  :: lchnk             ! chunk identifier
      real(r8), intent(in) :: relh(pcols)       ! relative humidity (>1)
      real(r8), intent(in) :: Ca(pcols)	        ! total internally mixed conc. (ug/m3)	
      real(r8), intent(in) :: Cnso4(pcols)	! so4(n) (ug/m3)
!TS++   The parameters below are required for best lognormal fit calculations
      real(r8), intent(in) :: Caitso4(pcols)	! so4(ait) (ug/m3) !TS
      real(r8), intent(in) :: Caitbc(pcols)	! BC(ait) (ug/m3) !TS
      real(r8), intent(in) :: Caitocbc(pcols)	! OC+BC(ait) (ug/m3) !TS
      real(r8), intent(in) :: Cas75(pcols)	! Sulfate mode 0.75um (ug/m3) !TS
      real(r8), intent(in) :: C_dst2(pcols)	! Dust (ug/m3)	    
      real(r8), intent(in) :: C_dst3(pcols)       ! Dust (ug/m3)
      real(r8), intent(in) :: C_ss1(pcols)	! Sea salt 1 (ug/m3)	    
      real(r8), intent(in) :: C_ss2(pcols)       ! Sea salt 2 (ug/m3)
      real(r8), intent(in) :: C_ss3(pcols)	! Sea salt 3 (ug/m3)
      real(r8), intent(in) :: rhoocbc(pcols)	! OC+BC density
!TS--
      real(r8), intent(in) :: Cnbc(pcols)	! BC(n) (ug/m3)	
      real(r8), intent(in) :: Cabce(pcols)	! BC(ax) (ug/m3)	    
      real(r8), intent(in) :: Cnoc(pcols)       ! OC(n) (ug/m3)   
      real(r8), intent(in) :: f_c(pcols)        ! = (Cbc+Coc)/(Cbc+Coc+Cso4)
      real(r8), intent(in) :: f_bc(pcols)       ! = Cbc/(Cbc+Coc)
      real(r8), intent(in) :: f_aq(pcols)       ! = Cso4a2/(Cso4a1+Cso4a2)
      real(r8), intent(in) :: fnbc(pcols)       ! = Cbc/(Cbc+Coc) for BC&OC(n)
      real(r8), intent(in) :: faitbc(pcols)     ! = Cbc/(Cbc+Coc) for BC&OC(Ait)
      real(r8), intent(inout) :: Nnatk(pcols,0:nmodes) ! Modal number concentration
      real(r8), intent(inout) :: n1(pcols)
      real(r8), intent(inout) :: n2(pcols)
      real(r8), intent(inout) :: n3(pcols)
      real(r8), intent(inout) :: n4(pcols)
      real(r8), intent(inout) :: n5(pcols)
      real(r8), intent(inout) :: n6(pcols)
      real(r8), intent(inout) :: n7(pcols)
      real(r8), intent(inout) :: n8(pcols)
      real(r8), intent(inout) :: n9(pcols)
      real(r8), intent(inout) :: n10(pcols)
      real(r8), intent(inout) :: n11(pcols)
      real(r8), intent(inout) :: n12(pcols)
      real(r8), intent(inout) :: n13(pcols)
      real(r8), intent(inout) :: n14(pcols)
      real(r8), intent(inout) :: cxstot(pcols)
      real(r8), intent(inout) :: Ccn(pcols)
!cak+
      real(r8), intent(out) :: lthick(pcols,nbmodes)   ! coating leyer thickness
      real(r8), intent(out) :: fso4coat(pcols,nbmodes) ! SO4-bound mass fraction in the coating
!cak-
      integer, intent(in) :: ncol
      integer, intent(in) :: nlons
      integer, intent(in) :: ind(pcols)
      integer, intent(in) :: lev
      real(r8) eps, e , dndlrk(0:100) , dndlrkny(0:100) , dndlr(0:100), &
               dndlk(pcols,nbmodes,0:100) , nrd(pcols,nbmodes) , &
               nk12(pcols,nbmodes) , nk12oc(pcols,nbmodes) , nag(pcols,nbmodes) , &
               fcondk(nbmodes), fcoagk(nbmodes) , fcoagock(nbmodes) , faqk(nbmodes) , &
               nrdtot(pcols) , nk12tot(pcols) , nk12octot(pcols), nagtot(pcols)
      real(r8) cabck(nbmodes), caock(nbmodes) , caqsk(nbmodes) , ccondsk(nbmodes) , &
               cam(pcols,nbmodes) , facm(pcols,nbmodes) , fbcm(pcols,nbmodes) , &
               faqm(pcols,nbmodes) , fraccn(pcols,nmodes), rhonew(pcols,nmodes) 
      real(r8) cxs(pcols,nbmodes), cxsdummy(pcols,nbmodes), Newmass(pcols,nbmodes), &
               camnull(pcols), facmnull(pcols), Massratio(pcols,nbmodes), &
!TS++		
		lnsig(pcols,nmodes),	& !TS: ln of standard deviation for given mode
		logsig(pcols,nbmodes), 	& !TS: log of standard deviation for given mode
		rknew(pcols,nbmodes),   & !TS: New modal radius from look-up tables
!TS--
!cak+
                cams(pcols,nbmodes),    & !AK: cam*(process specific factor for sulfate)
                rhocams(pcols,nbmodes), & !AK: density of the (cams) sulfate mixture
                rhoc(pcols,nbmodes),    & !AK: density of the OC and BC mixture in modes 5-10
                athird,                 & !AK: 1/3
                rpart, rcore,           & !AK: dry radii of particle and nonsoluble core (um) 
                rhocore,                & !AK: density of the BC and DU mixture in modes 6 and 7 
                ctemp,                  & !AK: density (work variable)  
                constlthick               !AK: 10*(3*pi/4)**1/3
!cak-
      integer i , ic , kcomp , lon , long 

      parameter (e=2.718281828_r8, eps=1.0e-30_r8)
!cak28mai09      parameter (constlthick=13.3067004_r8, athird=0.333333333333333_r8)
      parameter (constlthick=6.203504909_r8, athird=0.333333333333333_r8)

 
!     Calculation of the apportionment of internally mixed SO4, BC and OC 
!     mass between the various mineral and sea-salt background modes.

!      call modalapp2d(lchnk,ncol,Nnatk,Ca,f_c,f_bc,f_aq,Cam,facm,fbcm,faqm)
      call modalapp2d(lchnk,ncol,nlons,ind,Nnatk,Ca,f_c,f_bc,f_aq,Cam,facm,fbcm,faqm)
 	
!     Initialize excess mass cxs and cxstot, wrt. maximum allowed internal mixing
      do lon=1,ncol
	do kcomp=1, 10
	 cxs(lon,kcomp) = 0.0_r8
	 logsig(lon,kcomp) = 0.0_r8
	 rknew(lon,kcomp) = 0.0_r8
	 Newmass(lon,kcomp) = 0.0_r8
	 Massratio(lon,kcomp) = 0.0_r8
	 cams(lon,kcomp) = 0.0_r8
	 rhocams(lon,kcomp) = 0.0_r8
 	 rhoc(lon,kcomp) = 0.0_r8
 	 lthick(lon,kcomp) = 0.0_r8
	end do
        cxstot(lon)   = 0.0_r8
        camnull(lon)  = 0.0_r8
        facmnull(lon) = 0.0_r8
        fraccn(lon,12)= 0.0_r8      ! BC(n) (ext.) is also assumed hydrophobic
        Ccn(lon)      = 0.0_r8      ! opt               
      enddo

!     Table look-up and interpolations of CCN for the background modes:

!     The BC(ax) mode (kcomp=0) is assumed hydrophobic, and is therefore skipped here

!     SO4(Ait), BC(Ait) and OC(Ait) modes: 
      do kcomp = 1,3  
!TS++
!cak30jun08	do lon=1,ncol
	do long=1,nlons
       	lon=ind(long)

        if(kcomp==1 .and. Nnatk(lon,kcomp) > 1.e-6_r8) then
!	--------------------------------------------------------------------------------------
!TS:    KCOMP=1:Calculating the standard deviation for aerosol mode SO4 aitken, simple method:
!	---------------------------------------------------------------------------------------
!TS:	Density of internal mixture (kcomp=1):
	rhonew(lon,kcomp)=rhopart(l_so4_na)                       ! = 1841. (sulfuric acid)
        cams(lon,kcomp)=cam(lon,kcomp)*(Msv/Mso4) ! H2SO4
	lnsig(lon,kcomp) = ((2._r8/9._r8)*log(6.e-9_r8*(Caitso4(lon)*Msv/Mso4+cams(lon,kcomp))/ &
	(pi*1.e6_r8*Nnatk(lon,kcomp)*rhopart(l_so4_na)*((2._r8*effsize(l_so4_n))**3))))**0.5_r8

	elseif(kcomp==2 .and. Nnatk(lon,kcomp) > 1.e-6_r8) then
!	--------------------------------------------------------------------------------------
!TS:    KCOMP=2:Calculating the standard deviation for aerosol mode BC aitken, simple method:
!	---------------------------------------------------------------------------------------
!TS:	Density for internal mixture (kcomp=2):
        cams(lon,kcomp)=cam(lon,kcomp)*(Msv/Mso4) ! H2SO4
	rhonew(lon,kcomp) = (Caitbc(lon)+cams(lon,kcomp))/ &
                (Caitbc(lon)/rhopart(l_bc_n)+cams(lon,kcomp)/rhopart(l_so4_na))
	lnsig(lon,kcomp) = ((2._r8/9._r8)*log(6.e-9_r8*(Caitbc(lon)+cams(lon,kcomp))/ &
	(pi*1.e6_r8*Nnatk(lon,kcomp)*rhonew(lon,kcomp)*((2._r8*effsize(l_bc_n))**3))))**0.5_r8
        rpart=constlthick*((Caitbc(lon)+cams(lon,kcomp))/(Nnatk(lon,kcomp)*rhonew(lon,kcomp)))**athird
        rcore=constlthick*(Caitbc(lon)/(Nnatk(lon,kcomp)*rhopart(l_bc_n)))**athird
        lthick(lon,kcomp)=rpart-rcore
!	--------------------------------------------------------------------------------------
!TS:    KCOMP=3: This is currently an empty mode, so no calculations necessary.
!	---------------------------------------------------------------------------------------
        end if
	end do

!---------------------------------------------------------------------------------------------
!TS:	intlogn1to3.F90 reads the modal radii and standard deviations for mode 1 to 3 from look-up 
!TS:	tables created by A. Kirkevag. This is a more sophisticated method of lognormal fiting.
!TS:    Using this method, mass is not necessarily conserved, so we have to calculate the "new" mass. 
!---------------------------------------------------------------------------------------------		
       if(version2) then
  	call intlog1to3(ncol,nlons,ind,lev,kcomp, &
          cam(1,kcomp),Nnatk(1,kcomp),cxs, logsig(1,kcomp), rknew(1,kcomp))	
	do long=1,nlons
       	lon=ind(long)
	if(Nnatk(lon,kcomp) > 1.e-6_r8) then
	  Newmass(lon,kcomp) = ((1._r8/6._r8)*pi*Nnatk(lon,kcomp)*rhonew(lon,kcomp)*(2*1.e-6_r8*rknew(lon,kcomp))**3)* &
		exp((9._r8/2._r8)*(log(10._r8**logsig(lon,kcomp)))**2)
	 if(kcomp==1 )then	
 	   Massratio(lon,kcomp) = 1.e15_r8*Newmass(lon,kcomp)/(Caitso4(lon)*Msv/Mso4+cams(lon,kcomp))
	 elseif(kcomp==2 )then
	   Massratio(lon,kcomp) = 1.e15_r8*Newmass(lon,kcomp)/(Caitbc(lon)+cams(lon,kcomp))
	 end if 
	end if
	end do
       endif

!TS: 	In this version of the model CCN is not obtained from look-up tables ----> 
!TS:	this call is not necessary
! 	call intccn1to3(ncol,nlons,ind,lev,kcomp,relh, &
!          cam(1,kcomp),Nnatk(1,kcomp),fraccn(1,kcomp),cxs)

      enddo  ! kcomp


!     BC&OC(Ait) mode:   ------ facm not valid here (=0). Use faitbc instead  
      do kcomp = 4,4
	do long=1,nlons
       	lon=ind(long)
	 if(Nnatk(lon,kcomp) > 1.e-6_r8) then
!TS	  Density for internal mixture (kcomp=4):
!         for a mixture of both H2SO4 and (NH4)2SO4:
          cams(lon,kcomp)=cam(lon,kcomp)*((1.0_r8-faqm(lon,kcomp))*Msv+faqm(lon,kcomp)*Ms)/Mso4
          rhocams(lon,kcomp)=((1.0_r8-faqm(lon,kcomp))*Msv+faqm(lon,kcomp)*Ms)/ &
            ((1.0_r8-faqm(lon,kcomp))*Msv/rhopart(l_so4_na)+faqm(lon,kcomp)*Ms/rhopart(l_so4_ac))
	  rhonew(lon,kcomp)= (Caitocbc(lon)+cams(lon,kcomp))/ &
            (Caitocbc(lon)/rhoocbc(lon)+cams(lon,kcomp)/rhocams(lon,kcomp)) 
	  lnsig(lon,kcomp) = ((2._r8/9._r8)*log(6.e-9_r8*(Caitocbc(lon)+cams(lon,kcomp))/ &
	  !(pi*1.e6*Nnatk(lon,kcomp)*rhonew(lon,kcomp)*((2._r8*effsize(l_om_n))**3))))**0.5_r8 !mode removed
	  (pi*1.e6_r8*Nnatk(lon,kcomp)*rhonew(lon,kcomp)*((2._r8*effsize(l_om_ni))**3))))**0.5_r8
          rpart=constlthick*((Caitocbc(lon)+cams(lon,kcomp))/(Nnatk(lon,kcomp)*rhonew(lon,kcomp)))**athird
          rcore=constlthick*(Caitocbc(lon)*faitbc(lon)/(Nnatk(lon,kcomp)*rhopart(l_bc_n)))**athird
          lthick(lon,kcomp)=rpart-rcore
          fso4coat(lon,kcomp)=cams(lon,kcomp)/((1.0_r8-faitbc(lon))*Caitocbc(lon)+cams(lon,kcomp)+eps)
          fso4coat(lon,kcomp)=max(0._r8,min(1._r8,fso4coat(lon,kcomp)))
         end if
	end do
!---------------------------------------------------------------------------------------------
!TS:	intlogn4.F90 reads the modal radii and standard deviations for mode 4 from look-up 
!TS:	tables created by A. Kirkevag. Improved, but not mass-conserving method. 
!---------------------------------------------------------------------------------------------	
        if(version2) then
 	 call intlog4(ncol,nlons,ind,lev,kcomp, cam(1,kcomp), &
	  Nnatk(1,kcomp),faitbc,faqm(1,kcomp),cxs, logsig(1,kcomp), rknew(1,kcomp))
	 do long=1,nlons
          lon=ind(long)
	  if(Nnatk(lon,kcomp) > 1.e-6_r8 ) then
	   Newmass(lon,kcomp) = ((1._r8/6._r8)*pi*Nnatk(lon,kcomp)*rhonew(lon,kcomp)*(2*1.e-6_r8*rknew(lon,kcomp))**3)* &
		exp((9._r8/2._r8)*(log(10._r8**logsig(lon,kcomp)))**2)
	   Massratio(lon,kcomp) = 1.e15_r8*Newmass(lon,kcomp)/(Caitocbc(lon)+cams(lon,kcomp))
	  end if 
	 end do
        endif

!	--------------------------------------------------------------------------------------
!TS:    KCOMP=4:Calculating the standard deviation for aitken mode OC/BC, simple method:
!	---------------------------------------------------------------------------------------
!TS: 	In this version of the model CCN is not obtained from look-up tables ----> 
!TS:	this call is not necessary
! 	call intccn4(ncol,nlons,ind,lev,kcomp,relh, &
!          cam(1,kcomp),Nnatk(1,kcomp),faitbc,faqm(1,kcomp),fraccn(1,kcomp),cxs)

      enddo  ! kcomp


!     SO4(Ait75) (5), mineral (6-7) and Sea-salt (8-10) modes:
      do kcomp = 5,10                                                    ! fjern kcomp-loekka og forenkle her?
	do long=1,nlons
       	lon=ind(long)

!         for a mixture of both H2SO4 and (NH4)2SO4 with OC and BC:
          cams(lon,kcomp)=cam(lon,kcomp)*(facm(lon,kcomp)+(1.0_r8-facm(lon,kcomp))* &
            ((1.0_r8-faqm(lon,kcomp))*Msv+faqm(lon,kcomp)*Ms)/Mso4)
          rhoc(lon,kcomp)=(fbcm(lon,kcomp)+1.0_r8)/(fbcm(lon,kcomp)/rhopart(l_bc_ac)+1.0_r8/rhopart(l_om_ac))
          rhocams(lon,kcomp)=(facm(lon,kcomp)*Mso4+(1.0_r8-faqm(lon,kcomp))*Msv+faqm(lon,kcomp)*Ms)/ &
            (facm(lon,kcomp)*Mso4/rhoc(lon,kcomp)+ &
            (1.0_r8-faqm(lon,kcomp))*Msv/rhopart(l_so4_na)+faqm(lon,kcomp)*Ms/rhopart(l_so4_ac))

	if(kcomp==5 .and. Nnatk(lon,kcomp) > 1.e-6_r8) then
        rhonew(lon,kcomp)=(Cas75(lon)*Ms/Mso4+cams(lon,kcomp))/ &
                (Cas75(lon)*(Ms/Mso4)/rhopart(l_so4_pr)+cams(lon,kcomp)/rhocams(lon,kcomp))
	lnsig(lon,kcomp) = ((2._r8/9._r8)*log(6.e-9_r8*(Cas75(lon)*Ms/Mso4+cams(lon,kcomp))/ &
	        (pi*1.e6_r8*Nnatk(lon,kcomp)*rhonew(lon,kcomp)*((2.*effsize(l_so4_pr))**3))))**0.5_r8
        rpart=constlthick*((Cas75(lon)*Ms/Mso4+cams(lon,kcomp))/(Nnatk(lon,kcomp)*rhonew(lon,kcomp)))**athird
        rcore=constlthick*(cam(lon,kcomp)*facm(lon,kcomp)*fbcm(lon,kcomp)/(Nnatk(lon,kcomp)*rhopart(l_bc_n)))**athird
        lthick(lon,kcomp)=rpart-rcore
        fso4coat(lon,kcomp)=(Cas75(lon)*Ms/Mso4+cams(lon,kcomp)-cam(lon,kcomp)*facm(lon,kcomp))&
                           /(Cas75(lon)*Ms/Mso4+cams(lon,kcomp)-cam(lon,kcomp)*facm(lon,kcomp)*fbcm(lon,kcomp)+eps)
        fso4coat(lon,kcomp)=max(0._r8,min(1._r8,fso4coat(lon,kcomp)))

	elseif(kcomp==6 .and. C_dst2(lon) > 1.e-50_r8 .and. Nnatk(lon,kcomp) > 1.e-50_r8) then
        rhonew(lon,kcomp)=(C_dst2(lon)+cams(lon,kcomp))/(C_dst2(lon)/rhopart(l_dst_a2)+cams(lon,kcomp)/rhocams(lon,kcomp))
	lnsig(lon,kcomp) = ((2._r8/9._r8)*log(6.e-9*(C_dst2(lon)+cams(lon,kcomp))/ &
	(pi*1.e6_r8*Nnatk(lon,kcomp)*rhonew(lon,kcomp)*((2._r8*effsize(l_dst_a2))**3))))**0.5_r8
         if(Nnatk(lon,kcomp) > 1.e-7_r8) then
        ctemp=C_dst2(lon)+cam(lon,kcomp)*facm(lon,kcomp)*fbcm(lon,kcomp)
        rhocore=ctemp/(C_dst2(lon)/rhopart(l_dst_a2)+cam(lon,kcomp)*facm(lon,kcomp)*fbcm(lon,kcomp)/rhopart(l_bc_n))
        rpart=constlthick*((C_dst2(lon)+cams(lon,kcomp))/(Nnatk(lon,kcomp)*rhonew(lon,kcomp)))**athird
        rcore=constlthick*((C_dst2(lon)+cam(lon,kcomp)*facm(lon,kcomp)*fbcm(lon,kcomp))/(Nnatk(lon,kcomp)*rhocore))**athird
        lthick(lon,kcomp)=rpart-rcore
         endif
        fso4coat(lon,kcomp)=(cams(lon,kcomp)-cam(lon,kcomp)*facm(lon,kcomp))&
                           /(cams(lon,kcomp)-cam(lon,kcomp)*facm(lon,kcomp)*fbcm(lon,kcomp)+eps)
        fso4coat(lon,kcomp)=max(0._r8,min(1._r8,fso4coat(lon,kcomp)))

	elseif(kcomp==7 .and. C_dst3(lon) > 1.e-50_r8 .and. Nnatk(lon,kcomp) > 1.e-50_r8) then
        rhonew(lon,kcomp)=(C_dst3(lon)+cams(lon,kcomp))/(C_dst3(lon)/rhopart(l_dst_a3)+cams(lon,kcomp)/rhocams(lon,kcomp))
	lnsig(lon,kcomp) = ((2./9._r8)*log(6.e-9_r8*(C_dst3(lon)+cams(lon,kcomp))/ &
	(pi*1.e6_r8*Nnatk(lon,kcomp)*rhonew(lon,kcomp)*((2._r8*effsize(l_dst_a3))**3))))**0.5_r8
         if(Nnatk(lon,kcomp) > 1.e-7_r8) then
        ctemp=C_dst3(lon)+cam(lon,kcomp)*facm(lon,kcomp)*fbcm(lon,kcomp)
        rhocore=ctemp/(C_dst3(lon)/rhopart(l_dst_a3)+cam(lon,kcomp)*facm(lon,kcomp)*fbcm(lon,kcomp)/rhopart(l_bc_n))
        rpart=constlthick*((C_dst3(lon)+cams(lon,kcomp))/(Nnatk(lon,kcomp)*rhonew(lon,kcomp)))**athird
        rcore=constlthick*((C_dst3(lon)+cam(lon,kcomp)*facm(lon,kcomp)*fbcm(lon,kcomp))/(Nnatk(lon,kcomp)*rhocore))**athird
        lthick(lon,kcomp)=rpart-rcore

        fso4coat(lon,kcomp)=(cams(lon,kcomp)-cam(lon,kcomp)*facm(lon,kcomp))&
                           /(cams(lon,kcomp)-cam(lon,kcomp)*facm(lon,kcomp)*fbcm(lon,kcomp)+eps)
        fso4coat(lon,kcomp)=max(0._r8,min(1._r8,fso4coat(lon,kcomp)))
         endif       
 
	elseif(kcomp==8 .and. C_ss1(lon) > 1.e-50_r8 .and. Nnatk(lon,kcomp) > 1.e-50_r8) then
        rhonew(lon,kcomp)=(C_ss1(lon)+cams(lon,kcomp))/(C_ss1(lon)/rhopart(l_ss_a1)+cams(lon,kcomp)/rhocams(lon,kcomp))
	lnsig(lon,kcomp) = ((2._r8/9._r8)*log(6.e-9_r8*(C_ss1(lon)+cams(lon,kcomp))/ &
	(pi*1.e6_r8*Nnatk(lon,kcomp)*rhonew(lon,kcomp)*((2._r8*effsize(l_ss_a1))**3))))**0.5_r8

	elseif(kcomp==9 .and. C_ss2(lon) > 1.e-50_r8 .and. Nnatk(lon,kcomp) > 1.e-50_r8) then
        rhonew(lon,kcomp)=(C_ss2(lon)+cams(lon,kcomp))/(C_ss2(lon)/rhopart(l_ss_a2)+cams(lon,kcomp)/rhocams(lon,kcomp))
	lnsig(lon,kcomp) = ((2._r8/9._r8)*log(6.e-9_r8*(C_ss2(lon)+cams(lon,kcomp))/ &
	(pi*1.e6_r8*Nnatk(lon,kcomp)*rhonew(lon,kcomp)*((2._r8*effsize(l_ss_a2))**3))))**0.5_r8

	elseif(kcomp==10 .and. C_ss3(lon) > 1.e-50_r8 .and. Nnatk(lon,kcomp) > 1.e-50_r8) then
        rhonew(lon,kcomp)=(C_ss3(lon)+cams(lon,kcomp))/(C_ss3(lon)/rhopart(l_ss_a3)+cams(lon,kcomp)/rhocams(lon,kcomp))
	lnsig(lon,kcomp) = ((2._r8/9._r8)*log(6.e-9_r8*(C_ss3(lon)+cams(lon,kcomp))/ &
	(pi*1.e6_r8*Nnatk(lon,kcomp)*rhonew(lon,kcomp)*((2._r8*effsize(l_ss_a3))**3))))**0.5_r8

	end if

       end do  ! long

       if(version2) then
	call intlog5to10(ncol,nlons,ind,lev,kcomp, &
          cam(1,kcomp),Nnatk(1,kcomp),facm(1,kcomp),fbcm(1,kcomp), &
          faqm(1,kcomp),cxs, logsig(1,kcomp), rknew(1,kcomp))
	do long=1,nlons
       	 lon=ind(long)
	 if(Nnatk(lon,kcomp) > 1.e-6_r8) then
	   Newmass(lon,kcomp) = ((1._r8/6._r8)*pi*Nnatk(lon,kcomp)*rhonew(lon,kcomp)*(2*1.e-6_r8*rknew(lon,kcomp))**3)* &
		exp((9._r8/2._r8)*(log(10._r8**logsig(lon,kcomp)))**2)
	   if(kcomp==5)then
	     Massratio(lon,kcomp) = 1.e15_r8*Newmass(lon,kcomp)/(Cas75(lon)*Ms/Mso4+cams(lon,kcomp))
	   elseif(kcomp==6) then
	     Massratio(lon,kcomp) = 1.e15_r8*Newmass(lon,kcomp)/(C_dst2(lon)+cams(lon,kcomp))
	   elseif(kcomp==7)then
	     Massratio(lon,kcomp) = 1.e15_r8*Newmass(lon,kcomp)/(C_dst3(lon)+cams(lon,kcomp))
	   elseif(kcomp==8)then
	     Massratio(lon,kcomp) = 1.e15_r8*Newmass(lon,kcomp)/(C_ss1(lon)+cams(lon,kcomp))
	   elseif(kcomp==9)then
	     Massratio(lon,kcomp) = 1.e15_r8*Newmass(lon,kcomp)/(C_ss2(lon)+cams(lon,kcomp))
	   elseif(kcomp==10)then
	     Massratio(lon,kcomp) = 1.e15_r8*Newmass(lon,kcomp)/(C_ss3(lon)+cams(lon,kcomp))
	   end if 
	 end if
	end do
       endif

!TS: 	In this version of the model CCN is not obtained from look-up tables ----> 
!TS:	these calls are not necessary
! 	call intccn5to10(ncol,nlons,ind,lev,kcomp,relh, &
!          cam(1,kcomp),Nnatk(1,kcomp),facm(1,kcomp),fbcm(1,kcomp), &
!          faqm(1,kcomp),fraccn(1,kcomp),cxs)

      enddo  ! kcomp

!     SO4(n) and OC(n) modes:
!      do kcomp = 11,13,2
! 	call intccn1to3(ncol,nlons,ind,lev,kcomp-10,relh, &
!          camnull,Nnatk(1,kcomp),fraccn(1,kcomp),cxsdummy)
!      enddo
!     BC&OC(n) mode:   ------ facm not valid here (=0). Use fnbc instead  
!      do kcomp = 14,14
! 	call intccn4(ncol,nlons,ind,lev,kcomp-10,relh, &
!          camnull,Nnatk(1,kcomp),fnbc,faqm(1,kcomp-10),fraccn(1,kcomp),cxsdummy)
!      enddo

      if(version2) then          ! boer det vaere slik for version2=.true.?
!     Calculate excess mass cxs summed over all modes, cxstot:
       do kcomp = 1 , nbmodes
	do long=1,nlons
       	lon=ind(long)
          cxstot(lon)=cxstot(lon)+cxs(lon,kcomp)
        enddo
       enddo 
!     Lumping excess internally mixed mass over to external aitken modes:
	do long=1,nlons
       	lon=ind(long)
        Nnatk(lon,11) = efact_so4n*(Cnso4(lon)+cxstot(lon)*(1.0_r8-f_c(lon)))
        Nnatk(lon,12) = efact_bcn *(Cnbc(lon) +cxstot(lon)*f_c(lon)*f_bc(lon))
        Nnatk(lon,13) = efact_omn *(Cnoc(lon) +cxstot(lon)*f_c(lon)*(1.0_r8-f_bc(lon)))
       enddo
      endif

      do long=1,nlons
      lon=ind(long)
	lnsig(lon,11) = log(sgpart(l_so4_n))
	lnsig(lon,12) = log(sgpart(l_bc_n))
 	!lnsig(lon,13) = log(sgpart(l_om_n)) !mode removed
 	lnsig(lon,13) = log(sgpart(l_om_ni))
 	lnsig(lon,14) = log(sgpart(l_om_ni))
        n1(lon) =Nnatk(lon,1)
        n2(lon) =Nnatk(lon,2)
        n3(lon) =Nnatk(lon,3)
        n4(lon) =Nnatk(lon,4)
        n5(lon) =Nnatk(lon,5)
        n6(lon) =Nnatk(lon,6)
        n7(lon) =Nnatk(lon,7)
        n8(lon) =Nnatk(lon,8)
        n9(lon) =Nnatk(lon,9)
        n10(lon)=Nnatk(lon,10)
        n11(lon)=Nnatk(lon,11)
        n12(lon)=Nnatk(lon,12)
        n13(lon)=Nnatk(lon,13)
        n14(lon)=Nnatk(lon,14)
      enddo

      return
      end 
