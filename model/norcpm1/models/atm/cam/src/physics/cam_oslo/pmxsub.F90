subroutine pmxsub(lchnk   ,ncol    ,pint    ,pmid    ,iant    , coszrs , &
                  t       ,qm1     ,rh_temp ,Nnatk   , &
!cakx             ssatot  ,asymtot ,betot   ,Ctot    ,deltah_km, cltot)
                  ssatot  ,asymtot ,betot   ,Ctot    ,deltah_km, & 
#ifdef CMIP6
                  volc_ext_sun, volc_omega_sun, volc_g_sun, &
                  volc_ext_earth, volc_omega_earth, & 
#endif
                  taue550, taua550 ,volcmmr)

! Optical parameters for a composite aerosol is calculated by interpolation  
! from the tables kcomp1.out-kcomp14.out.
! Optimized June 2002 byrild Burud/NoSerC 
! Optimized July 2002 by Egil Storen/NoSerC (ces)
! Revised for inclusion of OC and modified aerosol backgeound aerosol 
! by Alf Kirkevaag in 2003, and finally rewritten for CAM3 February 2005.
! Modified for new aerosol schemes by Alf Kirkevaag in January 2006.

   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_history,  only: outfld
   use constituents, only: pcnst
   use physconst,    only: rair,pi
   use radconstants, only: nlwbands
   use opttab
   use const
   use aerosoldef
   use optinterpol,       only: interpol0,interpol1to3,interpol4,interpol5to10
   implicit none

!
! Input arguments
!
   integer, intent(in) :: lchnk                   ! chunk identifier
   integer, intent(in) :: ncol                    ! number of atmospheric columns
   integer, intent(in) :: iant                    ! 0=natural, 1=total aerosol
   real(r8), intent(in) :: coszrs(pcols)          ! Cosine solar zenith angle
   real(r8), intent(in) :: pint(pcols,pverp)      ! Model interface pressures (10*Pa)
   real(r8), intent(in) :: pmid(pcols,pver)       ! Model level pressures (Pa)
   real(r8), intent(in) :: t(pcols,pver)          ! Model level temperatures (K)
   real(r8), intent(in) :: qm1(pcols,pver,pcnst) ! Specific humidity and tracers (kg/kg)
#ifdef CMIP6
   real(r8), intent(in) :: volc_ext_sun(pcols,1:pver,nbands) ! volcanic aerosol extinction for solar bands, CMIP6
   real(r8), intent(in) :: volc_omega_sun(pcols,1:pver,nbands) ! volcanic aerosol SSA for solar bands, CMIP6
   real(r8), intent(in) :: volc_g_sun(pcols,1:pver,nbands) ! volcanic aerosol g for solar bands, CMIP6
   real(r8), intent(in) :: volc_ext_earth(pcols,1:pver,nlwbands) ! volcanic aerosol extinction for terrestrial bands, CMIP6
   real(r8), intent(in) :: volc_omega_earth(pcols,1:pver,nlwbands) ! volcanic aerosol SSA for terrestrial bands, CMIP6
#endif
   real(r8), intent(in) :: rh_temp(pcols,pver)    ! level relative humidity (fraction)
!cakx   real(r8), intent(in) :: cltot(pcols)       ! Total random overlap cloud cover
!
! Input-output arguments
!
   real(r8), intent(inout) :: Nnatk(pcols,pver,0:nmodes)! aerosol mode number concentration  
!
! Output arguments
!
   real(r8), intent(out) :: ssatot(pcols,pver,nbands) ! spectral aerosol single scattering albedo
   real(r8), intent(out) :: asymtot(pcols,pver,nbands)! spectral aerosol asymmetry factor
   real(r8), intent(out) :: betot(pcols,pver,nbands)  ! spectral aerosol extinction coefficient
   real(r8), intent(out) :: Ctot(pcols,pver)          ! aerosol mass concentration
   real(r8), intent(out) :: deltah_km(pcols,pver)     ! Layer thickness, unit km
   real(r8), intent(out) :: taue550(pcols)            ! AOD vis
   real(r8), intent(out) :: taua550(pcols)            ! absorptive AOD vis
   real(r8), intent(in)  :: volcmmr(pcols,pver)
!
!---------------------------Local variables-----------------------------
!
   integer  i, k, ib, kcomp, icol, mplus10
   logical  daylight(pcols)        ! only daylight calculations in interpol* if .true.

   real(r8) rhum(pcols,pver)       ! (trimmed) relative humidity for the aerosol calculations
#ifdef CMIP6
   real(r8) aodvisvolc(pcols)         ! AOD vis for CMIP6 volcanic aerosol
   real(r8) absvisvolc(pcols)         ! AAOD vis for CMIP6 volcanic aerosol
#endif 
   real(r8) deltah, airmass(pcols,pver) 
   real(r8) Ca(pcols,pver), f_c(pcols,pver), f_bc(pcols,pver), f_aq(pcols,pver)
   real(r8) fnbc(pcols,pver), faitbc(pcols,pver), vnbc, vaitbc
   real(r8) Cnso4(pcols,pver), Caso4(pcols,pver), Cnbc(pcols,pver), Cabc(pcols,pver), & 
            Cabce(pcols,pver), Cnoc(pcols,pver), Caoc(pcols,pver), dCtot(pcols,pver)
   real(r8) Cas75(pcols,pver)
   real(r8) Caitso4(pcols,pver), Cnbcioc(pcols,pver), Caitbc(pcols,pver), &
            Caitbcioc(pcols,pver), Cnocibc(pcols,pver), Caitoc(pcols,pver), &
            Caitocibc(pcols,pver), rhobcocn(pcols,pver), rhobcocait(pcols,pver)
   real(r8) Cas1, Cas2, Casc, totant, efact_bcocait, efact_bcocn 
   real(r8) cabck(nbmodes), caock(nbmodes), caqsk(nbmodes), ccondsk(nbmodes)
   real(r8) Cam(pcols,pver,nbmodes), fbcm(pcols,pver,nbmodes), fcm(pcols,pver,nbmodes), &
            faqm(pcols,pver,nbmodes), camnull(pcols,pver,nbmodes), faqm4(pcols,pver) 
   real(r8) C_ss1(pcols,pver), C_ss2(pcols,pver), C_ss3(pcols,pver), &
            C_dst1(pcols,pver), C_dst2(pcols,pver), C_dst3(pcols,pver) 
   real(r8) akso4c(pcols), akbcc(pcols), akocc(pcols)
   real(r8) ssa(pcols,pver,0:nmodes,nbands), asym(pcols,pver,0:nmodes,nbands), & 
            be(pcols,pver,0:nmodes,nbands), ke(pcols,pver,0:nmodes,nbands),    &
            betot550(pcols,pver), batot550(pcols,pver)
!cakx            betot550(pcols,pver), batot550(pcols,pver), taue550(pcols), taua550(pcols)
   real(r8) bekc1(pcols,pver), bekc2(pcols,pver), bekc3(pcols,pver), bekc4(pcols,pver), & 
            bekc5(pcols,pver), bekc6(pcols,pver), bekc7(pcols,pver), bekc8(pcols,pver), &
            bekc9(pcols,pver), bekc10(pcols,pver), bekc11(pcols,pver), &
            bekc12(pcols,pver), bekc13(pcols,pver), bekc14(pcols,pver), bekc0(pcols,pver)   
   real(r8) taukc1(pcols), taukc2(pcols), taukc3(pcols), taukc4(pcols), taukc5(pcols),  & 
            taukc6(pcols), taukc7(pcols), taukc8(pcols), taukc9(pcols), taukc10(pcols), & 
            taukc11(pcols), taukc12(pcols), taukc13(pcols), taukc14(pcols), taukc0(pcols)   
   real(r8) rh0(pcols,pver), rhoda(pcols,pver), cxstot(pcols,pver), akcxs(pcols) 
   real(r8) wak(pcols,pver), gak(pcols,pver), bak(pcols,pver), dayfoc(pcols,pver)
   real(r8) n_aerorig(pcols,pver), n_aer(pcols,pver)
!2oct07
   real(r8) rhobcbyoc
!2oct07
!cak_sept09
!cakx   real(r8) clearodvis(pcols), clearabsvis(pcols), cloudfree(pcols)
!cak_sept09

#undef AEROCOM
#ifdef AEROCOM 
   real(r8) cam4cond(pcols,pver), cam4aq(pcols,pver), cam4oc(pcols,pver), cam4bc(pcols,pver), &
            cam5cond(pcols,pver), cam5aq(pcols,pver), cam5oc(pcols,pver), cam5bc(pcols,pver)
   real(r8) Ctotdry(pcols,pver), Cwater(pcols,pver), mmr_aerh2o(pcols,pver), &
            dod550dry(pcols), abs550dry(pcols)
   real(r8) daerh2o(pcols),  dload(pcols,0:nmodes), dload3d(pcols,pver,0:nmodes), &
            dload_mi(pcols), dload_ss(pcols)
   real(r8) cmin(pcols,pver), cseas(pcols,pver)
   real(r8) nnat_1(pcols,pver), nnat_2(pcols,pver), nnat_3(pcols,pver), &
            nnat_4(pcols,pver), nnat_5(pcols,pver), nnat_6(pcols,pver), &
            nnat_7(pcols,pver), nnat_8(pcols,pver), nnat_9(pcols,pver), &
            nnat_10(pcols,pver),nnat_11(pcols,pver), nnat_12(pcols,pver), &
            nnat_13(pcols,pver), nnat_14(pcols,pver), nnat_0(pcols,pver)
   real(r8) ck(pcols,pver,0:nmodes), cknorm(pcols,pver,0:nmodes), &
            cknlt05(pcols,pver,0:nmodes), ckngt125(pcols,pver,0:nmodes)
   real(r8) aaerosn(pcols,pver,nbmp1:nmodes), aaeroln(pcols,pver,nbmp1:nmodes), &
            vaerosn(pcols,pver,nbmp1:nmodes), vaeroln(pcols,pver,nbmp1:nmodes), &
            aaeros(pcols,pver,0:nbmodes), aaerol(pcols,pver,0:nbmodes), & 
            vaeros(pcols,pver,0:nbmodes), vaerol(pcols,pver,0:nbmodes) 
   real(r8) cintbg(pcols,pver,0:nbmodes), &
            cintbg05(pcols,pver,0:nbmodes), cintbg125(pcols,pver,0:nbmodes), &
            cintbc(pcols,pver,0:nbmodes), &
            cintbc05(pcols,pver,0:nbmodes), cintbc125(pcols,pver,0:nbmodes), &  
            cintoc(pcols,pver,0:nbmodes), &
            cintoc05(pcols,pver,0:nbmodes), cintoc125(pcols,pver,0:nbmodes), &
            cintsc(pcols,pver,0:nbmodes), &
            cintsc05(pcols,pver,0:nbmodes), cintsc125(pcols,pver,0:nbmodes), &        
            cintsa(pcols,pver,0:nbmodes), &
            cintsa05(pcols,pver,0:nbmodes), cintsa125(pcols,pver,0:nbmodes)
   real(r8) c_mi(pcols,pver), c_mi05(pcols,pver), c_mi125(pcols,pver), &
            c_ss(pcols,pver), c_ss05(pcols,pver), c_ss125(pcols,pver), &
            c_bc(pcols,pver), c_bc05(pcols,pver), c_bc125(pcols,pver), &
            c_oc(pcols,pver), c_oc05(pcols,pver), c_oc125(pcols,pver), &
            c_sa(pcols,pver), c_sa05(pcols,pver), c_sa125(pcols,pver), &
            c_sc(pcols,pver), c_sc05(pcols,pver), c_sc125(pcols,pver), &
            c_s4(pcols,pver), c_s405(pcols,pver), c_s4125(pcols,pver)
   real(r8) aaeros_tot(pcols,pver), aaerol_tot(pcols,pver), vaeros_tot(pcols,pver), &
            vaerol_tot(pcols,pver), aaercols(pcols), aaercoll(pcols), vaercols(pcols), & 
            vaercoll(pcols), derlt05(pcols), dergt05(pcols), der(pcols), &
            erlt053d(pcols,pver), ergt053d(pcols,pver), er3d(pcols,pver)
!old   real(r8) bext1(pcols,pver,0:nbmodes), babs1(pcols,pver,0:nbmodes), bext2(pcols,pver,0:nbmodes), & 
!old            babs2(pcols,pver,0:nbmodes), bebg1(pcols,pver,0:nbmodes), babg1(pcols,pver,0:nbmodes), & 
!old            bebg2(pcols,pver,0:nbmodes), babg2(pcols,pver,0:nbmodes), bebc1(pcols,pver,0:nbmodes), & 
!old            babc1(pcols,pver,0:nbmodes), bebc2(pcols,pver,0:nbmodes), babc2(pcols,pver,0:nbmodes), & 
!old            beoc1(pcols,pver,0:nbmodes), baoc1(pcols,pver,0:nbmodes), beoc2(pcols,pver,0:nbmodes), & 
!old            baoc2(pcols,pver,0:nbmodes), bes41(pcols,pver,0:nbmodes), bas41(pcols,pver,0:nbmodes), & 
!old            bes42(pcols,pver,0:nbmodes), bas42(pcols,pver,0:nbmodes)
   real(r8) bext440(pcols,pver,0:nbmodes), babs440(pcols,pver,0:nbmodes), &
            bext500(pcols,pver,0:nbmodes), babs500(pcols,pver,0:nbmodes), &
            bext550(pcols,pver,0:nbmodes), babs550(pcols,pver,0:nbmodes), &
            bext670(pcols,pver,0:nbmodes), babs670(pcols,pver,0:nbmodes), &
            bext870(pcols,pver,0:nbmodes), babs870(pcols,pver,0:nbmodes), &
            bebg440(pcols,pver,0:nbmodes), babg440(pcols,pver,0:nbmodes), &
            bebg500(pcols,pver,0:nbmodes), babg500(pcols,pver,0:nbmodes), &
            bebg550(pcols,pver,0:nbmodes), babg550(pcols,pver,0:nbmodes), &
            bebg670(pcols,pver,0:nbmodes), babg670(pcols,pver,0:nbmodes), &
            bebg870(pcols,pver,0:nbmodes), babg870(pcols,pver,0:nbmodes), &
            bebc440(pcols,pver,0:nbmodes), babc440(pcols,pver,0:nbmodes), &
            bebc500(pcols,pver,0:nbmodes), babc500(pcols,pver,0:nbmodes), &
            bebc550(pcols,pver,0:nbmodes), babc550(pcols,pver,0:nbmodes), &
            bebc670(pcols,pver,0:nbmodes), babc670(pcols,pver,0:nbmodes), &
            bebc870(pcols,pver,0:nbmodes), babc870(pcols,pver,0:nbmodes), &
            beoc440(pcols,pver,0:nbmodes), baoc440(pcols,pver,0:nbmodes), &
            beoc500(pcols,pver,0:nbmodes), baoc500(pcols,pver,0:nbmodes), &
            beoc550(pcols,pver,0:nbmodes), baoc550(pcols,pver,0:nbmodes), &
            beoc670(pcols,pver,0:nbmodes), baoc670(pcols,pver,0:nbmodes), &
            beoc870(pcols,pver,0:nbmodes), baoc870(pcols,pver,0:nbmodes), &
            besu440(pcols,pver,0:nbmodes), basu440(pcols,pver,0:nbmodes), &
            besu500(pcols,pver,0:nbmodes), basu500(pcols,pver,0:nbmodes), &
            besu550(pcols,pver,0:nbmodes), basu550(pcols,pver,0:nbmodes), &
            besu670(pcols,pver,0:nbmodes), basu670(pcols,pver,0:nbmodes), &
            besu870(pcols,pver,0:nbmodes), basu870(pcols,pver,0:nbmodes)
   real(r8) bebglt1(pcols,pver,0:nbmodes), bebggt1(pcols,pver,0:nbmodes), &
            bebclt1(pcols,pver,0:nbmodes), bebcgt1(pcols,pver,0:nbmodes), & 
            beoclt1(pcols,pver,0:nbmodes), beocgt1(pcols,pver,0:nbmodes), & 
            bes4lt1(pcols,pver,0:nbmodes), bes4gt1(pcols,pver,0:nbmodes), &
            backsc550(pcols,pver,0:nbmodes), backsc550x(pcols,pver,nbmp1:nmodes), &
            backsc550tot(pcols,pver), ec550_aer(pcols,pver), abs550_aer(pcols,pver), &
            bs550_aer(pcols,pver), ec550dry_aer(pcols,pver), abs550dry_aer(pcols,pver)  
!old   real(r8) bext1n(pcols,pver,0:nbmodes), babs1n(pcols,pver,0:nbmodes), bext2n(pcols,pver,0:nbmodes), & 
!old            babs2n(pcols,pver,0:nbmodes), bebg1n(pcols,pver,0:nbmodes), babg1n(pcols,pver,0:nbmodes), & 
!old            bebg2n(pcols,pver,0:nbmodes), babg2n(pcols,pver,0:nbmodes), bebc1n(pcols,pver,0:nbmodes), & 
!old            babc1n(pcols,pver,0:nbmodes), bebc2n(pcols,pver,0:nbmodes), babc2n(pcols,pver,0:nbmodes), & 
!old            beoc1n(pcols,pver,0:nbmodes), baoc1n(pcols,pver,0:nbmodes), beoc2n(pcols,pver,0:nbmodes), & 
!old            baoc2n(pcols,pver,0:nbmodes), bes41n(pcols,pver,0:nbmodes), bas41n(pcols,pver,0:nbmodes), & 
!old            bes42n(pcols,pver,0:nbmodes), bas42n(pcols,pver,0:nbmodes)
   real(r8) bext440n(pcols,pver,0:nbmodes), babs440n(pcols,pver,0:nbmodes), &
            bext500n(pcols,pver,0:nbmodes), babs500n(pcols,pver,0:nbmodes), &
            bext550n(pcols,pver,0:nbmodes), babs550n(pcols,pver,0:nbmodes), &
            bext670n(pcols,pver,0:nbmodes), babs670n(pcols,pver,0:nbmodes), &
            bext870n(pcols,pver,0:nbmodes), babs870n(pcols,pver,0:nbmodes), &
            bebg440n(pcols,pver,0:nbmodes), babg440n(pcols,pver,0:nbmodes), &
            bebg500n(pcols,pver,0:nbmodes), babg500n(pcols,pver,0:nbmodes), &
            bebg550n(pcols,pver,0:nbmodes), babg550n(pcols,pver,0:nbmodes), &
            bebg670n(pcols,pver,0:nbmodes), babg670n(pcols,pver,0:nbmodes), &
            bebg870n(pcols,pver,0:nbmodes), babg870n(pcols,pver,0:nbmodes), &
            bebc440n(pcols,pver,0:nbmodes), babc440n(pcols,pver,0:nbmodes), &
            bebc500n(pcols,pver,0:nbmodes), babc500n(pcols,pver,0:nbmodes), &
            bebc550n(pcols,pver,0:nbmodes), babc550n(pcols,pver,0:nbmodes), &
            bebc670n(pcols,pver,0:nbmodes), babc670n(pcols,pver,0:nbmodes), &
            bebc870n(pcols,pver,0:nbmodes), babc870n(pcols,pver,0:nbmodes), &
            beoc440n(pcols,pver,0:nbmodes), baoc440n(pcols,pver,0:nbmodes), &
            beoc500n(pcols,pver,0:nbmodes), baoc500n(pcols,pver,0:nbmodes), &
            beoc550n(pcols,pver,0:nbmodes), baoc550n(pcols,pver,0:nbmodes), &
            beoc670n(pcols,pver,0:nbmodes), baoc670n(pcols,pver,0:nbmodes), &
            beoc870n(pcols,pver,0:nbmodes), baoc870n(pcols,pver,0:nbmodes), &
            besu440n(pcols,pver,0:nbmodes), basu440n(pcols,pver,0:nbmodes), &
            besu500n(pcols,pver,0:nbmodes), basu500n(pcols,pver,0:nbmodes), &
            besu550n(pcols,pver,0:nbmodes), basu550n(pcols,pver,0:nbmodes), &
            besu670n(pcols,pver,0:nbmodes), basu670n(pcols,pver,0:nbmodes), &
            besu870n(pcols,pver,0:nbmodes), basu870n(pcols,pver,0:nbmodes)
   real(r8) bebglt1n(pcols,pver,0:nbmodes), bebggt1n(pcols,pver,0:nbmodes), &
            bebclt1n(pcols,pver,0:nbmodes), bebcgt1n(pcols,pver,0:nbmodes), & 
            beoclt1n(pcols,pver,0:nbmodes), beocgt1n(pcols,pver,0:nbmodes), & 
            bes4lt1n(pcols,pver,0:nbmodes), bes4gt1n(pcols,pver,0:nbmodes), &
            backsc550n(pcols,pver,0:nbmodes) 
!old
!   real(r8) bext1tot(pcols,pver), babs1tot(pcols,pver), bext2tot(pcols,pver), & 
!            babs2tot(pcols,pver), bebg1tot(pcols,pver), babg1tot(pcols,pver), &
!            bebg2tot(pcols,pver), babg2tot(pcols,pver), bebc1tot(pcols,pver), &
!            babc1tot(pcols,pver), bebc2tot(pcols,pver), babc2tot(pcols,pver), &
!            beoc1tot(pcols,pver), baoc1tot(pcols,pver), beoc2tot(pcols,pver), &
!            baoc2tot(pcols,pver), bes41tot(pcols,pver), bas41tot(pcols,pver), & 
!            bes42tot(pcols,pver), bas42tot(pcols,pver)
!old
   real(r8) bext440tot(pcols,pver), babs440tot(pcols,pver), &
            bext500tot(pcols,pver), babs500tot(pcols,pver), &
            bext550tot(pcols,pver), babs550tot(pcols,pver), &
            bext670tot(pcols,pver), babs670tot(pcols,pver), &
            bext870tot(pcols,pver), babs870tot(pcols,pver), &
            bebg440tot(pcols,pver), babg440tot(pcols,pver), &
            bebg500tot(pcols,pver), babg500tot(pcols,pver), &
            bebg550tot(pcols,pver), babg550tot(pcols,pver), &
            bebg670tot(pcols,pver), babg670tot(pcols,pver), &
            bebg870tot(pcols,pver), babg870tot(pcols,pver), &
            bebc440tot(pcols,pver), babc440tot(pcols,pver), &
            bebc500tot(pcols,pver), babc500tot(pcols,pver), &
            bebc550tot(pcols,pver), babc550tot(pcols,pver), &
            bebc670tot(pcols,pver), babc670tot(pcols,pver), &
            bebc870tot(pcols,pver), babc870tot(pcols,pver), &
            beoc440tot(pcols,pver), baoc440tot(pcols,pver), &
            beoc500tot(pcols,pver), baoc500tot(pcols,pver), &
            beoc550tot(pcols,pver), baoc550tot(pcols,pver), &
            beoc670tot(pcols,pver), baoc670tot(pcols,pver), &
            beoc870tot(pcols,pver), baoc870tot(pcols,pver), &
            besu440tot(pcols,pver), basu440tot(pcols,pver), &
            besu500tot(pcols,pver), basu500tot(pcols,pver), &
            besu550tot(pcols,pver), basu550tot(pcols,pver), &
            besu670tot(pcols,pver), basu670tot(pcols,pver), &
            besu870tot(pcols,pver), basu870tot(pcols,pver)
   real(r8) bebglt1t(pcols,pver), bebggt1t(pcols,pver), bebclt1t(pcols,pver), & 
            bebcgt1t(pcols,pver), beoclt1t(pcols,pver), beocgt1t(pcols,pver), &
            bes4lt1t(pcols,pver), bes4gt1t(pcols,pver)
!old   real(r8) be1x(pcols,pver,nbmp1:nmodes), ba1x(pcols,pver,nbmp1:nmodes), &
!old            be2x(pcols,pver,nbmp1:nmodes), ba2x(pcols,pver,nbmp1:nmodes),   &
   real(r8) be440x(pcols,pver,nbmp1:nmodes), ba440x(pcols,pver,nbmp1:nmodes), &
            be500x(pcols,pver,nbmp1:nmodes), ba500x(pcols,pver,nbmp1:nmodes), &
            be550x(pcols,pver,nbmp1:nmodes), ba550x(pcols,pver,nbmp1:nmodes), &
            be670x(pcols,pver,nbmp1:nmodes), ba670x(pcols,pver,nbmp1:nmodes), &
            be870x(pcols,pver,nbmp1:nmodes), ba870x(pcols,pver,nbmp1:nmodes), &
            belt1x(pcols,pver,nbmp1:nmodes), begt1x(pcols,pver,nbmp1:nmodes)
!old   real(r8) bes41xt(pcols,pver), bas41xt(pcols,pver), bes42xt(pcols,pver), &
!old            bas42xt(pcols,pver), bebc1xt(pcols,pver), babc1xt(pcols,pver), &
!old            bebc2xt(pcols,pver), babc2xt(pcols,pver), beoc1xt(pcols,pver), &
!old            baoc1xt(pcols,pver), beoc2xt(pcols,pver), baoc2xt(pcols,pver) 
   real(r8) besu440xt(pcols,pver),basu440xt(pcols,pver), &
            besu500xt(pcols,pver),basu500xt(pcols,pver), &
            besu550xt(pcols,pver),basu550xt(pcols,pver), &
            besu670xt(pcols,pver),basu670xt(pcols,pver), &
            besu870xt(pcols,pver),basu870xt(pcols,pver), &
            bebc440xt(pcols,pver),babc440xt(pcols,pver), &
            bebc500xt(pcols,pver),babc500xt(pcols,pver), &
            bebc550xt(pcols,pver),babc550xt(pcols,pver), &
            bebc670xt(pcols,pver),babc670xt(pcols,pver), &
            bebc870xt(pcols,pver),babc870xt(pcols,pver), &
            beoc440xt(pcols,pver),baoc440xt(pcols,pver), &
            beoc500xt(pcols,pver),baoc500xt(pcols,pver), &
            beoc550xt(pcols,pver),baoc550xt(pcols,pver), &
            beoc670xt(pcols,pver),baoc670xt(pcols,pver), &
            beoc870xt(pcols,pver),baoc870xt(pcols,pver) 
   real(r8) bs4lt1xt(pcols,pver), bs4gt1xt(pcols,pver), bbclt1xt(pcols,pver), &
            bbcgt1xt(pcols,pver), boclt1xt(pcols,pver), bocgt1xt(pcols,pver)
!old   real(r8) bint1du(pcols,pver), bint2du(pcols,pver), bint1mi(pcols,pver), &
!old            bint2mi(pcols,pver), bint1ss(pcols,pver), bint2ss(pcols,pver), &
!old            baintdu(pcols,pver), baintmi(pcols,pver), baintss(pcols,pver)
   real(r8) bint440du(pcols,pver), bint500du(pcols,pver), bint550du(pcols,pver), &
            bint670du(pcols,pver), bint870du(pcols,pver), &
            bint440ss(pcols,pver), bint500ss(pcols,pver), bint550ss(pcols,pver), &
            bint670ss(pcols,pver), bint870ss(pcols,pver), &
            baint550du(pcols,pver), baint550ss(pcols,pver)
   real(r8) bedustlt1(pcols,pver), bedustgt1(pcols,pver), &
            besslt1(pcols,pver), bessgt1(pcols,pver)
   real(r8) dod4403d(pcols,pver), abs4403d(pcols,pver), &
            dod4403d_ss(pcols,pver),   & ! abs4403d_ss(pcols,pver), & 
            dod4403d_dust(pcols,pver), & ! abs4403d_dust(pcols,pver), &
            dod4403d_so4(pcols,pver),  & ! abs4403d_so4(pcols,pver), &
            dod4403d_bc(pcols,pver),   & ! abs4403d_bc(pcols,pver), &
            dod4403d_pom(pcols,pver),  & ! abs4403d_pom(pcols,pver), &
            dod5003d(pcols,pver), abs5003d(pcols,pver), &
            dod5003d_ss(pcols,pver),   & ! abs5003d_ss(pcols,pver), & 
            dod5003d_dust(pcols,pver), & ! abs5003d_dust(pcols,pver), &
            dod5003d_so4(pcols,pver),  & ! abs5003d_so4(pcols,pver), &
            dod5003d_bc(pcols,pver),   & ! abs5003d_bc(pcols,pver), &
            dod5003d_pom(pcols,pver),  & ! abs5003d_pom(pcols,pver), &
            dod5503d(pcols,pver), abs5503d(pcols,pver), abs5503dalt(pcols,pver), &
            dod5503d_ss(pcols,pver), abs5503d_ss(pcols,pver), & 
            dod5503d_dust(pcols,pver), abs5503d_dust(pcols,pver), &
            dod5503d_so4(pcols,pver), abs5503d_so4(pcols,pver), &
            dod5503d_bc(pcols,pver), abs5503d_bc(pcols,pver), &
            dod5503d_pom(pcols,pver), abs5503d_pom(pcols,pver), &
            dod6703d(pcols,pver), abs6703d(pcols,pver), &
            dod6703d_ss(pcols,pver),   & ! abs6703d_ss(pcols,pver), & 
            dod6703d_dust(pcols,pver), & ! abs6703d_dust(pcols,pver), &
            dod6703d_so4(pcols,pver),  & ! abs6703d_so4(pcols,pver), &
            dod6703d_bc(pcols,pver),   & ! abs6703d_bc(pcols,pver), &
            dod6703d_pom(pcols,pver),  & ! abs6703d_pom(pcols,pver), &
            dod8703d(pcols,pver), abs8703d(pcols,pver), &
            dod8703d_ss(pcols,pver),   & ! abs8703d_ss(pcols,pver), & 
            dod8703d_dust(pcols,pver), & ! abs8703d_dust(pcols,pver), &
            dod8703d_so4(pcols,pver),  & ! abs8703d_so4(pcols,pver), &
            dod8703d_bc(pcols,pver),   & ! abs8703d_bc(pcols,pver), &
            dod8703d_pom(pcols,pver) ! abs8703d_pom(pcols,pver)
!old            dod8653d_ss(pcols,pver), dod8653d_dust(pcols,pver), &
!old            dod8653d_so4(pcols,pver), &
!old            dod8653d_bc(pcols,pver), dod8653d_pom(pcols,pver)
   real(r8) dod5503dlt1_ss(pcols,pver), dod5503dgt1_ss(pcols,pver), &
            dod5503dlt1_dust(pcols,pver), dod5503dgt1_dust(pcols,pver), &
            dod5503dlt1_so4(pcols,pver), dod5503dgt1_so4(pcols,pver), &
            dod5503dlt1_bc(pcols,pver), dod5503dgt1_bc(pcols,pver), &
            dod5503dlt1_pom(pcols,pver), dod5503dgt1_pom(pcols,pver)
!old   real(r8) dod550(pcols), abs550(pcols), dod865_ss(pcols), &
!old            dod865_dust(pcols), &
!old            dod865_so4(pcols), dod865_bc(pcols), dod865_pom(pcols), &
   real(r8) dod440(pcols), abs440(pcols), dod500(pcols), abs500(pcols),  &
            dod550(pcols), abs550(pcols), abs550alt(pcols), dod670(pcols),& 
            abs670(pcols), dod870(pcols), abs870(pcols),                 &
            dod440_ss(pcols), dod440_dust(pcols), dod440_so4(pcols),     &
            dod440_bc(pcols), dod440_pom(pcols),                         & 
            dod500_ss(pcols), dod500_dust(pcols), dod500_so4(pcols),     &
            dod500_bc(pcols), dod500_pom(pcols),                         &        
            dod550_ss(pcols), dod550_dust(pcols), dod550_so4(pcols),     &
            dod550_bc(pcols), dod550_pom(pcols),                         &        
            dod670_ss(pcols), dod670_dust(pcols), dod670_so4(pcols),     &
            dod670_bc(pcols), dod670_pom(pcols),                         &        
            dod870_ss(pcols), dod870_dust(pcols), dod870_so4(pcols),     &
            dod870_bc(pcols), dod870_pom(pcols),                         &        


            dod550lt1_ss(pcols), dod550gt1_ss(pcols), dod550lt1_dust(pcols), & 
            dod550gt1_dust(pcols), dod550lt1_so4(pcols), &
            dod550gt1_so4(pcols), dod550lt1_bc(pcols), dod550gt1_bc(pcols), &
            dod550lt1_pom(pcols), dod550gt1_pom(pcols)
!old   real(r8) dod550_ss(pcols), dod550_dust(pcols), &
!old            dod550_so4(pcols), dod550_bc(pcols), dod550_pom(pcols)
! Added 021007 OS
   real(r8) abs550_ss(pcols), abs550_dust(pcols), &
            abs550_so4(pcols), abs550_bc(pcols), abs550_pom(pcols)
!old   real(r8) clearod550(pcols),clearod865(pcols),clearabs550(pcols),cloudfree(pcols)
   real(r8) clearod550(pcols),clearod870(pcols),clearabs550(pcols),clearabs550alt(pcols)
#endif  ! AEROCOM


!
!-------------------------------------------------------------------------
!
      do k=1,pver
         do icol=1,ncol
!         Set upper and lower relative humidity for the aerosol calculations
!old          rhum(icol,k) = min(0.98_r8, max(rh_temp(icol,k), 0.01_r8))
          rhum(icol,k) = min(0.995_r8, max(rh_temp(icol,k), 0.01_r8))
!test          rhum(icol,k) = 0.01_r8
          rhoda(icol,k) = pmid(icol,k)/(rair*t(icol,k))      ! unit kg/m^3
         end do
      end do

!     interpol-calculations only when daylight or not:
#ifdef AEROCOM
      do icol=1,ncol
        daylight(icol) = .true.
      end do
#else
      do icol=1,ncol
        if (coszrs(icol) > 0.0_r8) then
          daylight(icol) = .true.
        else
          daylight(icol) = .false.
        endif
      end do
#endif  ! AEROCOM

!     Set SO4, BC and OC concentrations:

!     initialize concentration fields
      do i=0,nmodes
        do k=1,pver
          do icol=1,ncol
            Nnatk(icol,k,i)  = 0.0_r8
          end do
        end do 
      end do 
      do k=1,pver
        do icol=1,ncol
          Cnso4(icol,k)     = 0.0_r8
          Cas75(icol,k)     = 0.0_r8
          Caso4(icol,k)     = 0.0_r8
   	  Caitso4(icol,k)   = 0.0_r8
          Cnbc(icol,k)      = 0.0_r8
          Caitbc(icol,k)    = 0.0_r8
  	  Cnbcioc(icol,k)   = 0.0_r8 
     	  Caitbcioc(icol,k) = 0.0_r8 
          Cabc(icol,k)      = 0.0_r8
          Cabce(icol,k)     = 0.0_r8
          Cnoc(icol,k)      = 0.0_r8
          Caoc(icol,k)      = 0.0_r8
	  Cnocibc(icol,k)   = 0.0_r8 
     	  Caitoc(icol,k)    = 0.0_r8
     	  Caitocibc(icol,k) = 0.0_r8
          f_c(icol,k)       = 0.0_r8
          f_bc(icol,k)      = 0.0_r8
          f_aq(icol,k)      = 0.0_r8
          Ca(icol,k)        = 0.0_r8
          n_aerorig(icol,k) = 0.0_r8
          n_aer(icol,k)     = 0.0_r8
        end do
      end do


!      and then calculate the concentration fields,
!      first sea-salt and mineral background:
       do k=1,pver
         do icol=1,ncol

          C_ss1(icol,k)  = 1.e9_r8*qm1(icol,k,l_ss_a1)*rhoda(icol,k)
          C_ss2(icol,k)  = 1.e9_r8*qm1(icol,k,l_ss_a2)*rhoda(icol,k)
          C_ss3(icol,k)  = 1.e9_r8*qm1(icol,k,l_ss_a3)*rhoda(icol,k)
          C_dst1(icol,k) = 0.0_r8  ! non-existent in the aerocomB emissions
          C_dst2(icol,k) = 1.e9_r8*qm1(icol,k,l_dst_a2)*rhoda(icol,k)
          C_dst3(icol,k) = 1.e9_r8*qm1(icol,k,l_dst_a3)*rhoda(icol,k)
         end do
       end do
!      and then added SO4, BC and OC (from the life cycle module): 
       if(iant==1) then
        call outfld('N_AERORG',n_aerorig,pcols,lchnk)
       do k=1,pver
         do icol=1,ncol
!         Note: concentrations of SO4 = 3 * S, here with unit ug/m^3:
!SO4
          Cnso4(icol,k)     = 3.e9_r8*qm1(icol,k,l_so4_n)*rhoda(icol,k)   ! SO4(n) mass
 	  Cas75(icol,k)     = 3.e9_r8*qm1(icol,k,l_so4_pr)*rhoda(icol,k)  + &
          1.e9_r8*volcmmr(icol,k)*rhoda(icol,k)! SO4(Ait75)
   	  Caitso4(icol,k)   = 3.e9_r8*qm1(icol,k,l_so4_na)*rhoda(icol,k)  ! SO4(Ait) mass (without cond. SO4) 
          Cas1              = 3.e9_r8*qm1(icol,k,l_so4_a1)*rhoda(icol,k)  ! SO4 condensation mass
          Cas2              = 3.e9_r8*qm1(icol,k,l_so4_a2)*rhoda(icol,k)  ! SO4 from cloud processing
          Casc              = 3.e9_r8*qm1(icol,k,l_so4_ac)*rhoda(icol,k)  ! SO4 coagulation mass
          Caso4(icol,k)     = Cas1+Cas2+Casc                           ! internally mixed SO4 mass
!BC
          Cnbc(icol,k)      = 1.e9_r8*qm1(icol,k,l_bc_n)*rhoda(icol,k)    ! BC(n) mass
          Cabce(icol,k)     = 1.e9_r8*qm1(icol,k,l_bc_ax)*rhoda(icol,k)   ! BC(ax) mass
  	  Cnbcioc(icol,k)   = 1.e9_r8*qm1(icol,k,l_bc_ni)*rhoda(icol,k)   ! BC mass for BC&OC(n) mode
          Caitbc(icol,k)    = 1.e9_r8*qm1(icol,k,l_bc_a)*rhoda(icol,k)    ! BC(Ait) mass
     	  Caitbcioc(icol,k) = 1.e9_r8*qm1(icol,k,l_bc_ai)*rhoda(icol,k)   ! BC mass for BC&OC(Ait) mode
          Cabc(icol,k)      = 1.e9_r8*qm1(icol,k,l_bc_ac)*rhoda(icol,k)   ! BC coagulation mass
!OC
!          Cnoc(icol,k)      = 1.e9_r8*qm1(icol,k,l_om_n)*rhoda(icol,k)    ! OC(n) mass
          Cnoc(icol,k)= eps
	  Cnocibc(icol,k)   = 1.e9_r8*qm1(icol,k,l_om_ni)*rhoda(icol,k)   ! OC mass for BC&OC(n) mode 
     	  Caitoc(icol,k)    = 0._r8    ! Aitken mode ff om included together with biomass
!1.e9_r8*qm1(icol,k,l_om_a)*rhoda(icol,k)    ! OC(Ait) mass 
     	  Caitocibc(icol,k) = 1.e9_r8*qm1(icol,k,l_om_ai)*rhoda(icol,k)   ! OC mass for BC&OC(Ait) mode
          Caoc(icol,k)      = 1.e9_r8*qm1(icol,k,l_om_ac)*rhoda(icol,k)   ! OC coagulation mass
!BC&OC
          rhobcocn(icol,k)  = (Cnbcioc(icol,k)+Cnocibc(icol,k)) &      ! mass density for the BC&OC(n) mode 
                              /(Cnbcioc(icol,k)/rhopart(l_bc_ni)+Cnocibc(icol,k)/rhopart(l_om_ni)+eps)
          rhobcocait(icol,k)= (Caitbcioc(icol,k)+Caitocibc(icol,k)) &  ! mass density for the BC&OC(Ait) mode 
                              /(Caitbcioc(icol,k)/rhopart(l_bc_ni)+Caitocibc(icol,k)/rhopart(l_om_ni)+eps)
!ALL 
          Ca(icol,k)        = Caso4(icol,k)+Cabc(icol,k)+Caoc(icol,k)
          f_c(icol,k)       = min((Cabc(icol,k)+Caoc(icol,k))/(Ca(icol,k)+eps),0.99999_r8)
	  f_bc(icol,k)      = min(Cabc(icol,k)/(Cabc(icol,k)+Caoc(icol,k)+eps),0.99999_r8)
          f_aq(icol,k)      = min(Cas2/(Caso4(icol,k)+eps),0.99999_r8)
          fnbc(icol,k)      = min(Cnbcioc(icol,k)/(Cnbcioc(icol,k)+Cnocibc(icol,k)+eps),0.99999_r8)         
          faitbc(icol,k)    = min(Caitbcioc(icol,k)/(Caitbcioc(icol,k)+Caitocibc(icol,k)+eps),0.99999_r8)
         end do
       end do
      else  ! iant==0
       do k=1,pver
          do icol=1,ncol
           Cnso4(icol,k)     = 0.0_r8
  	   Cas75(icol,k)     = 0.0_r8
           Caso4(icol,k)     = 0.0_r8
    	   Caitso4(icol,k)   = 0.0_r8
           Cnbc(icol,k)      = 0.0_r8
           Cabc(icol,k)      = 0.0_r8
   	   Cnbcioc(icol,k)   = 0.0_r8
           Caitbc(icol,k)    = 0.0_r8
     	   Caitbcioc(icol,k) = 0.0_r8
           Cabce(icol,k)     = 0.0_r8
           Cnoc(icol,k)      = 0.0_r8
           Caoc(icol,k)      = 0.0_r8
  	   Cnocibc(icol,k)   = 0.0_r8
     	   Caitoc(icol,k)    = 0.0_r8
     	   Caitocibc(icol,k) = 0.0_r8
           rhobcocn(icol,k)  = 1500._r8
           rhobcocait(icol,k)= 1500._r8
           f_c(icol,k)       = 0.0_r8
           f_bc(icol,k)      = 0.0_r8
           f_aq(icol,k)      = 0.75_r8
           Ca(icol,k)        = 0.0_r8
           fnbc(icol,k)      = 0.0_r8
           faitbc(icol,k)    = 0.0_r8
         end do
       end do
      endif

!     Avoid very small numbers
      do k=1,pver
         do icol=1,ncol
          Ca(icol,k)        = max(eps,Ca(icol,k))
          f_c(icol,k)       = max(eps,f_c(icol,k))
	  f_bc(icol,k)      = max(eps,f_bc(icol,k)) 
          f_aq(icol,k)      = max(eps,f_aq(icol,k)) 
          fnbc(icol,k)      = max(eps,fnbc(icol,k))
          faitbc(icol,k)    = max(eps,faitbc(icol,k)) 
         end do
       end do

!     Determine number concentration of the BC(ax) mode (0), Aitken modes SO4, BC and OC (1-4)
!     and mineral (6-7) and sea-salt modes (8-10):
!      write(6,*) 'writebcax',efact_bcax
      do k=1,pver
         do icol=1,ncol
           Nnatk(icol,k,0)  = Cabce(icol,k)*efact_bcax
           Nnatk(icol,k,1)  = Caitso4(icol,k)*efact_so4na
           Nnatk(icol,k,2)  = Caitbc(icol,k)*efact_bca
!          Nnatk(icol,k,3)  = Caitoc(icol,k)*efact_oma
           Nnatk(icol,k,3)  = 0.0_r8
           efact_bcocait    = 1e-15/(e**(4.5_r8*(log(sgpart(l_om_ai)))**2) &
             *(4.0_r8/3.0_r8)*pi*(effsize(l_om_ai))**3*rhobcocait(icol,k)+eps)
           Nnatk(icol,k,4)  = (Caitbcioc(icol,k)+Caitocibc(icol,k))*efact_bcocait
           Nnatk(icol,k,5)  = Cas75(icol,k)*efact_so4pr
           Nnatk(icol,k,6)  = C_dst2(icol,k)*efact_dst2
           Nnatk(icol,k,7)  = C_dst3(icol,k)*efact_dst3
           Nnatk(icol,k,8)  = C_ss1(icol,k)*efact_ss1
           Nnatk(icol,k,9)  = C_ss2(icol,k)*efact_ss2
           Nnatk(icol,k,10) = C_ss3(icol,k)*efact_ss3
         end do
      end do
!     and for the externally mixed SO4 and BC modes:
      do k=1,pver
         do icol=1,ncol
           Nnatk(icol,k,11)=Cnso4(icol,k)*efact_so4n
           Nnatk(icol,k,12)=Cnbc(icol,k)*efact_bcn
           Nnatk(icol,k,13)=Cnoc(icol,k)*efact_omn
           efact_bcocn     = 1e-15/(e**(4.5_r8*(log(sgpart(l_om_ni)))**2) &
              *(4.0_r8/3.0_r8)*pi*(effsize(l_om_ni))**3*rhobcocn(icol,k)+eps) 
           Nnatk(icol,k,14) =(Cnbcioc(icol,k)+Cnocibc(icol,k))*efact_bcocn
         end do
      end do
      do i=0,nmodes    ! mode 0 to 14 
        do k=1,pver
          do icol=1,ncol
!later (after modalapp-redist.): n_aer(icol,k)=n_aer(icol,k)+Nnatk(icol,k,i)
            n_aerorig(icol,k)=n_aerorig(icol,k)+Nnatk(icol,k,i)
          end do
        enddo  
      enddo

!     Calculation of the apportionment of internally mixed SO4, BC and OC 
!     mass between the various mineral and sea-salt background modes.

      call modalapp (lchnk,ncol,Nnatk,Ca,f_c,f_bc,f_aq,Cam,fcm,fbcm,faqm)
 
!#ifdef AEROCOM      
!      do k=1,pver
!        do icol=1,ncol
!          cam4cond(icol,k) = (1.0_r8-fcm(icol,k,4))*(1.0_r8-faqm(icol,k,4))*Cam(icol,k,4)
!          cam4aq(icol,k)   = (1.0_r8-fcm(icol,k,4))*faqm(icol,k,4)*Cam(icol,k,4)
!          cam4oc(icol,k)   = fcm(icol,k,4)*(1.0_r8-fbcm(icol,k,4))*Cam(icol,k,4)
!          cam4bc(icol,k)   = fcm(icol,k,4)*fbcm(icol,k,4)*Cam(icol,k,4)
!          cam5cond(icol,k) = (1.0_r8-fcm(icol,k,5))*(1.0_r8-faqm(icol,k,5))*Cam(icol,k,5)
!          cam5aq(icol,k)   = (1.0_r8-fcm(icol,k,5))*faqm(icol,k,5)*Cam(icol,k,5)
!          cam5oc(icol,k)   = fcm(icol,k,5)*(1.0_r8-fbcm(icol,k,5))*Cam(icol,k,5)
!          cam5bc(icol,k)   = fcm(icol,k,5)*fbcm(icol,k,5)*Cam(icol,k,5)
!        end do
!      enddo  
!#endif

      do k=1,pver
        do icol=1,ncol
          faqm4(icol,k) = faqm(icol,k,4)
        end do
      enddo  
      do i=1,nbmodes
        do k=1,pver
          do icol=1,ncol
            camnull(icol,k,i) = 0.0_r8
          enddo
        enddo
      enddo

!     AeroCom: Find dry aerosol mass for subsequent calculation of
!     condensed water mass below...


#ifdef AEROCOM      

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      if(iant==1) then 

      do k=1,pver
        do icol=1,ncol
          Ctotdry(icol,k)=0.0_r8
          rh0(icol,k)=0.0_r8
          cxstot(icol,k) = 0.0_r8
        end do
      enddo  

!     BC(ax) mode (dry only): 
      call interpol0 (lchnk, ncol, daylight, Nnatk, ssa, asym, be, ke)

!     SO4(Ait), BC(Ait) and OC(Ait) modes:
      mplus10=0
      call interpol1to3 (lchnk, ncol, daylight, rh0, mplus10, Nnatk, Cam, &
                         ssa, asym, be, ke, cxstot)
!     BC&OC(Ait) mode:   ------ fcm not valid here (=0). Use faitbc instead
      mplus10=0
      call interpol4 (lchnk, ncol, daylight, rh0, mplus10, Nnatk, Cam, faitbc, faqm4, &
                      ssa, asym, be, ke, cxstot)

!     SO4(Ait75) (5), Mineral (6-7) and Sea-salt (8-10) modes:
      call interpol5to10 (lchnk, ncol, daylight, rh0, Nnatk, Cam, fcm, fbcm, faqm, &
                      ssa, asym, be, ke, cxstot)

!cak+   lumping excess internally mixed mass over to external aitken modes:
       do k=1,pver
!          write(*,*) 'dir: k, cxstot =', k, cxstot(1,k)
          do icol=1,ncol
           Nnatk(icol,k,11)=efact_so4n* &
              (Cnso4(icol,k)+cxstot(icol,k)*(1.0_r8-f_c(icol,k)))
           Nnatk(icol,k,12)=efact_bcn * &
              (Cnbc(icol,k) +cxstot(icol,k)*f_c(icol,k)*f_bc(icol,k))
           Nnatk(icol,k,13)=efact_omn * &
              (Cnoc(icol,k) +cxstot(icol,k)*f_c(icol,k)*(1.0_r8-f_bc(icol,k)))
          end do
       end do
!cak-

!     SO4(n), BC(n) and OC(n) modes:
      mplus10=1
      call interpol1to3 (lchnk, ncol, daylight, rh0, mplus10, Nnatk, camnull, &
                         ssa, asym, be, ke, cxstot)
!     BC&OC(n) mode:   ------ fcm not valid here (=0). Use fnbc instead
      mplus10=1
      call interpol4 (lchnk, ncol, daylight, rh0, mplus10, Nnatk, camnull, fnbc, faqm4, &
                      ssa, asym, be, ke, cxstot)

       do i=0,nmodes    ! mode 0 to 14 
         do k=1,pver
            do icol=1,ncol
!g             if(Nnatk(icol,k,i)>0.0_r8) then
               dCtot(icol,k)=1.e3_r8*be(icol,k,i,8)/(ke(icol,k,i,8)+eps)
               Ctotdry(icol,k)=Ctotdry(icol,k)+dCtot(icol,k)*Nnatk(icol,k,i) 
!g             endif
           end do
         enddo  
       enddo

      endif ! iant==1
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
#endif  ! AEROCOM

      do k=1,pver
        do icol=1,ncol
          cxstot(icol,k) = 0.0_r8
        end do
      enddo  
!     (Wet) Optical properties for each of the aerosol modes:

!     BC(ax) mode (dry only):
      call interpol0 (lchnk, ncol, daylight, Nnatk, ssa, asym, be, ke)

!     SO4(Ait), BC(Ait) and OC(Ait) modes:
      mplus10=0
      call interpol1to3 (lchnk, ncol, daylight, rhum, mplus10, Nnatk, Cam, &
                         ssa, asym, be, ke, cxstot)
!     BC&OC(Ait) mode:   ------ fcm invalid here (=0). Using faitbc instead

      mplus10=0
      call interpol4 (lchnk, ncol, daylight, rhum, mplus10, Nnatk, Cam, faitbc, faqm4,&
                      ssa, asym, be, ke, cxstot)


!     SO4(Ait75) (5), Mineral (6-7) and Sea-salt (8-10) modes:
      call interpol5to10 (lchnk, ncol, daylight, rhum, Nnatk, Cam, fcm, fbcm, faqm, &
                      ssa, asym, be, ke, cxstot)


#ifndef AEROCOM      
!cak+   lumping excess internally mixed mass over to external aitken modes:
       do k=1,pver
!          write(*,*) 'dir: k, cxstot =', k, cxstot(1,k)
          do icol=1,ncol
           Nnatk(icol,k,11)=efact_so4n* &
              (Cnso4(icol,k)+cxstot(icol,k)*(1.0_r8-f_c(icol,k)))
           Nnatk(icol,k,12)=efact_bcn * &
              (Cnbc(icol,k) +cxstot(icol,k)*f_c(icol,k)*f_bc(icol,k))
           Nnatk(icol,k,13)=efact_omn * &
              (Cnoc(icol,k) +cxstot(icol,k)*f_c(icol,k)*(1.0_r8-f_bc(icol,k)))
          end do
       end do
!cak-
#endif

!     SO4(n), BC(n) and OC(n) modes:
      mplus10=1
      call interpol1to3 (lchnk, ncol, daylight, rhum, mplus10, Nnatk, camnull, &  
                      ssa, asym, be, ke, cxstot)
!     BC&OC(n) mode:    ------ fcm not valid here (=0). Use fnbc instead
      mplus10=1
      call interpol4 (lchnk, ncol, daylight, rhum, mplus10, Nnatk, camnull, fnbc, faqm4, &
                     ssa, asym, be, ke, cxstot)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      do k=1,pver
        do icol=1,ncol
          Ctot(icol,k)=0.0_r8
        end do
      enddo  

      do i=0,nmodes    ! mode 0 to 14 
        do k=1,pver
          do icol=1,ncol
!g            if(Nnatk(icol,k,i)>0.0_r8) then
              dCtot(icol,k)=1.e3_r8*be(icol,k,i,8)/(ke(icol,k,i,8)+eps) 
              Ctot(icol,k)=Ctot(icol,k)+dCtot(icol,k)*Nnatk(icol,k,i) 
!g            endif
          end do
        enddo  
      enddo

#ifdef AEROCOM      
!     Mass concentration (ug/m3) and mmr (kg/kg) of aerosol condensed water:
      do k=1,pver
        do icol=1,ncol
          Cwater(icol,k)=Ctot(icol,k)-Ctotdry(icol,k)
          mmr_aerh2o(icol,k)=1.e-9_r8*Cwater(icol,k)/rhoda(icol,k)
        end do
      enddo  
#endif
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!     Optical properties of total aerosol:
      do ib=1,nbands
        do k=1,pver
         do icol=1,ncol
            betot(icol,k,ib)=0.0_r8
            ssatot(icol,k,ib)=0.0_r8
            asymtot(icol,k,ib)=0.0_r8
          end do
        enddo  
      enddo
      do ib=1,nbands
       do i=0,nmodes
        do k=1,pver
          do icol=1,ncol
!g            if(Nnatk(icol,k,i)>0.0_r8) then      
           betot(icol,k,ib)=betot(icol,k,ib)+Nnatk(icol,k,i)*be(icol,k,i,ib)
           ssatot(icol,k,ib)=ssatot(icol,k,ib)+Nnatk(icol,k,i) &
                             *be(icol,k,i,ib)*ssa(icol,k,i,ib)           
           asymtot(icol,k,ib)=asymtot(icol,k,ib)+Nnatk(icol,k,i) &
                         *be(icol,k,i,ib)*ssa(icol,k,i,ib)*asym(icol,k,i,ib)
!g            endif
          end do
        enddo
       enddo
      enddo
#ifdef CMIP6
!ak+   Adding also the volcanic contribution (CMIP6), which is using a CMIP6
!      band numbering identical to the AeroTab numbering (unlike CAM) both
!      for SW and LW. I.e., no remapping is required here.
!   Info from CMIP_CAM6_radiation_v3.nc
! wl1_sun = 0.2, 0.263158, 0.344828, 0.441501, 0.625, 0.77821, 1.24224, 
!    1.2987, 1.62602, 1.94175, 2.15054, 2.5, 3.07692, 3.84615 ;
! wl2_sun = 0.263158, 0.344828, 0.441501, 0.625, 0.77821, 1.24224, 1.2987, 
!    1.62602, 1.94175, 2.15054, 2.5, 3.07692, 3.84615, 12.1951 ;
! wl1_earth = 3.07692, 3.84615, 4.20168, 4.44444, 4.80769, 5.55556, 6.75676, 
!    7.19424, 8.47458, 9.25926, 10.2041, 12.1951, 14.2857, 15.873, 20, 28.5714 ;
! wl2_earth = 3.84615, 4.20168, 4.44444, 4.80769, 5.55556, 6.75676, 7.19424, 
!    8.47458, 9.25926, 10.2041, 12.1951, 14.2857, 15.873, 20, 28.5714, 1000 ;
      do ib=1,nbands
         betot(1:ncol,1:pver,ib) = betot(1:ncol,1:pver,ib) &
             + volc_ext_sun(1:ncol,1:pver,ib)
         ssatot(1:ncol,1:pver,ib) = ssatot(1:ncol,1:pver,ib) &
             + volc_ext_sun(1:ncol,1:pver,ib)*volc_omega_sun(1:ncol,1:pver,ib)
         asymtot(1:ncol,1:pver,ib) = asymtot(1:ncol,1:pver,ib) &
             + volc_ext_sun(1:ncol,1:pver,ib)*volc_omega_sun(1:ncol,1:pver,ib) &
              *volc_g_sun(1:ncol,1:pver,ib)
      enddo
!ak-  and then calculate the total bulk optical parameters
#endif
      do ib=1,nbands
       do k=1,pver
        do icol=1,ncol
          ssatot(icol,k,ib)=ssatot(icol,k,ib)/(betot(icol,k,ib)+eps)
          asymtot(icol,k,ib)=asymtot(icol,k,ib) &
                           /(betot(icol,k,ib)*ssatot(icol,k,ib)+eps)
        end do
       enddo
      enddo

!     initialize modal optical extinctions
      do k=1,pver
        do icol=1,ncol
          bekc0(icol,k)=0.0_r8
          bekc1(icol,k)=0.0_r8
          bekc2(icol,k)=0.0_r8
          bekc3(icol,k)=0.0_r8
          bekc4(icol,k)=0.0_r8
          bekc5(icol,k)=0.0_r8
          bekc6(icol,k)=0.0_r8
          bekc7(icol,k)=0.0_r8
          bekc8(icol,k)=0.0_r8
          bekc9(icol,k)=0.0_r8
          bekc10(icol,k)=0.0_r8
          bekc11(icol,k)=0.0_r8
          bekc12(icol,k)=0.0_r8
          bekc13(icol,k)=0.0_r8
          bekc14(icol,k)=0.0_r8
        end do
      enddo  

!     optical depth (in band 8 = vis.) for each of the 15 modes (for testutskrifter bare...)
      do k=1,pver
        do icol=1,ncol
          bekc0(icol,k) =Nnatk(icol,k,0) *be(icol,k,0,8)
          bekc1(icol,k) =Nnatk(icol,k,1) *be(icol,k,1,8)
          bekc2(icol,k) =Nnatk(icol,k,2) *be(icol,k,2,8)
          bekc3(icol,k) =Nnatk(icol,k,3) *be(icol,k,3,8)
          bekc4(icol,k) =Nnatk(icol,k,4) *be(icol,k,4,8)
          bekc5(icol,k) =Nnatk(icol,k,5) *be(icol,k,5,8)
          bekc6(icol,k) =Nnatk(icol,k,6) *be(icol,k,6,8)
          bekc7(icol,k) =Nnatk(icol,k,7) *be(icol,k,7,8)
          bekc8(icol,k) =Nnatk(icol,k,8) *be(icol,k,8,8)
          bekc9(icol,k) =Nnatk(icol,k,9) *be(icol,k,9,8)
          bekc10(icol,k)=Nnatk(icol,k,10)*be(icol,k,10,8)
          bekc11(icol,k)=Nnatk(icol,k,11)*be(icol,k,11,8)
          bekc12(icol,k)=Nnatk(icol,k,12)*be(icol,k,12,8)
          bekc13(icol,k)=Nnatk(icol,k,13)*be(icol,k,13,8)
          bekc14(icol,k)=Nnatk(icol,k,14)*be(icol,k,14,8)
        end do
      enddo

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!     APPROXIMATE aerosol extinction and absorption at 550nm (0.35-0.64)
      do k=1,pver
        do icol=1,ncol
          betot550(icol,k)=betot(icol,k,8)
          batot550(icol,k)=betot550(icol,k)*(1.0_r8-ssatot(icol,k,8))
        end do
      enddo  
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      do k=1,pver
        do icol=1,ncol
          wak(icol,k) = 0.0_r8
          gak(icol,k) = 0.0_r8
          bak(icol,k) = 0.0_r8
          dayfoc(icol,k) = 0.0_r8    
        enddo
      end do

      do k=1,pver
        do icol=1,ncol
!           dayfoc < 1 when looping only over gridcells with daylight 
          if(daylight(icol)) then
            dayfoc(icol,k) = 1.0_r8
!cak+ with the new bands in CAM, band 8 is still at ca 0.5 um (0.35-0.64)
            wak(icol,k) = ssatot(icol,k,8)
            gak(icol,k) = asymtot(icol,k,8)
            bak(icol,k) = betot(icol,k,8)
!cak-
          endif 
        enddo
      end do

      if(iant==1) then  ! The remainder of the routine run only when all aerosols are included

!       optical parameters in visible light (0.35-0.64um)
        call outfld('WAK     ',wak   ,pcols,lchnk)
        call outfld('GAK     ',gak   ,pcols,lchnk)
        call outfld('BAK     ',bak   ,pcols,lchnk)
        call outfld('DAYFOC  ',dayfoc,pcols,lchnk)
!       total aerosol number concentrations before and after lumping
        do i=0,nmodes    ! mode 0 to 14 
         do k=1,pver
          do icol=1,ncol
            n_aer(icol,k)=n_aer(icol,k)+Nnatk(icol,k,i)
          end do
         enddo  
        enddo

        call outfld('N_AER   ',n_aer ,pcols,lchnk)

!       Initialize fields
        do icol=1,ncol
          akso4c(icol)=0.0_r8
          akbcc(icol)=0.0_r8
          akocc(icol)=0.0_r8
          akcxs(icol)=0.0_r8
          taue550(icol)=0.0_r8
          taua550(icol)=0.0_r8
#ifdef CMIP6
          aodvisvolc(icol)=0.0_r8
          absvisvolc(icol)=0.0_r8
#endif 
!cakx          clearodvis(icol)=0.0_r8
!cakx          clearabsvis(icol)=0.0_r8
#ifdef AEROCOM
          cloudfree(icol)=0.0_r8
#endif
          taukc0(icol)=0.0_r8
          taukc1(icol)=0.0_r8
          taukc2(icol)=0.0_r8
          taukc3(icol)=0.0_r8
          taukc4(icol)=0.0_r8
          taukc5(icol)=0.0_r8
          taukc6(icol)=0.0_r8
          taukc7(icol)=0.0_r8
          taukc8(icol)=0.0_r8
          taukc9(icol)=0.0_r8
          taukc10(icol)=0.0_r8
          taukc11(icol)=0.0_r8
          taukc12(icol)=0.0_r8
          taukc13(icol)=0.0_r8
          taukc14(icol)=0.0_r8
        enddo

        do icol=1,ncol
         do k=1,pver
!         Layer thickness, unit km, and layer airmass, unit kg/m2
          deltah=1.e-4_r8*(pint(icol,k+1)-pint(icol,k))/(rhoda(icol,k)*9.8_r8)  
          deltah_km(icol,k)=deltah
          airmass(icol,k)=1.e3_r8*deltah*rhoda(icol,k)
#ifdef CMIP6
          aodvisvolc(icol)=aodvisvolc(icol)+volc_ext_sun(icol,k,4)*deltah
          absvisvolc(icol)=absvisvolc(icol)+volc_ext_sun(icol,k,4) &
                                *(1.0_r8-volc_omega_sun(icol,k,4))*deltah
#endif
!         Added mass due to nucleation, condensation, coagulation or aqueous chemistry
          akso4c(icol)=akso4c(icol)+deltah* &
!                       (Cnso4(icol,k)+Caitso4(icol,k)+Caso4(icol,k))
                       (Cnso4(icol,k)+Caitso4(icol,k)+Caso4(icol,k)+Cas75(icol,k))
          akbcc(icol) =akbcc(icol)+deltah*  &
                       (Cnbc(icol,k)+Caitbc(icol,k)+Cabc(icol,k) &
                       +Cnbcioc(icol,k)+Caitbcioc(icol,k)+Cabce(icol,k))
          akocc(icol) =akocc(icol)+deltah*  &
                       (Cnoc(icol,k)+Caitoc(icol,k)+Caoc(icol,k) &
                       +Cnocibc(icol,k)+Caitocibc(icol,k))
          akcxs(icol) =akcxs(icol)+cxstot(icol,k)*deltah
!          Optical depths at ca. 550 nm (0.35-0.64um) 
          taue550(icol)=taue550(icol)+betot550(icol,k)*deltah
          taua550(icol)=taua550(icol)+batot550(icol,k)*deltah
          taukc0(icol) =taukc0(icol) +bekc0(icol,k)*deltah
          taukc1(icol) =taukc1(icol) +bekc1(icol,k)*deltah
          taukc2(icol) =taukc2(icol) +bekc2(icol,k)*deltah
          taukc3(icol) =taukc3(icol) +bekc3(icol,k)*deltah
          taukc4(icol) =taukc4(icol) +bekc4(icol,k)*deltah
          taukc5(icol) =taukc5(icol) +bekc5(icol,k)*deltah
          taukc6(icol) =taukc6(icol) +bekc6(icol,k)*deltah
          taukc7(icol) =taukc7(icol) +bekc7(icol,k)*deltah
          taukc8(icol) =taukc8(icol) +bekc8(icol,k)*deltah
          taukc9(icol) =taukc9(icol) +bekc9(icol,k)*deltah
          taukc10(icol)=taukc10(icol)+bekc10(icol,k)*deltah
          taukc11(icol)=taukc11(icol)+bekc11(icol,k)*deltah
          taukc12(icol)=taukc12(icol)+bekc12(icol,k)*deltah
          taukc13(icol)=taukc13(icol)+bekc13(icol,k)*deltah
          taukc14(icol)=taukc14(icol)+bekc14(icol,k)*deltah
         end do  ! k
!        Total and absorption optical depth under cloudfree conditions
!cakx         cloudfree(icol)=1._r8-cltot(icol)
!cakx         clearodvis(icol)=cloudfree(icol)*taue550(icol)
!cakx         clearabsvis(icol)=cloudfree(icol)*taua550(icol)
        end do   ! icol

!      write(*,*) 'pmxsub 8, cltot =', cltot(1)
         
!        Aerosol mass columns in pmxsub.F
!        call outfld('AKSO4C  ',akso4c,pcols,lchnk)  
!        call outfld('AKBCC   ',akbcc ,pcols,lchnk)
!        call outfld('AKOCC   ',akocc,pcols,lchnk)
        call outfld('AKCXS   ',akcxs ,pcols,lchnk)
!       Extinction and absorption for 0.55 um for the total aerosol, and AODs 
#ifdef AEROCOM
        call outfld('BETOT550',betot550,pcols,lchnk)
        call outfld('BATOT550',batot550,pcols,lchnk)
#endif
#ifdef CMIP6
        call outfld('AODVVOLC',aodvisvolc ,pcols,lchnk)
        call outfld('ABSVVOLC',absvisvolc ,pcols,lchnk)
#endif
        call outfld('TAUE550 ',taue550,pcols,lchnk)
        call outfld('TAUA550 ',taua550,pcols,lchnk)
!       Added clear sky variables
!cakx        call outfld('CDODVIS ',clearodvis,pcols,lchnk)
!cakx        call outfld('CABSVIS ',clearabsvis,pcols,lchnk)
!cakx        call outfld('CLDFREE ',cloudfree,pcols,lchnk)
!        call outfld('TAUKC0  ',taukc0 ,pcols,lchnk)
!        call outfld('TAUKC1  ',taukc1 ,pcols,lchnk)
!        call outfld('TAUKC2  ',taukc2 ,pcols,lchnk)
!        call outfld('TAUKC3  ',taukc3 ,pcols,lchnk)
!        call outfld('TAUKC4  ',taukc4 ,pcols,lchnk)
!        call outfld('TAUKC5  ',taukc5 ,pcols,lchnk)
!        call outfld('TAUKC6  ',taukc6 ,pcols,lchnk)
!        call outfld('TAUKC7  ',taukc7 ,pcols,lchnk)
!        call outfld('TAUKC8  ',taukc8 ,pcols,lchnk)
!        call outfld('TAUKC9  ',taukc9 ,pcols,lchnk)
!        call outfld('TAUKC10 ',taukc10,pcols,lchnk)
!        call outfld('TAUKC11 ',taukc11,pcols,lchnk)
!        call outfld('TAUKC12 ',taukc12,pcols,lchnk)
!        call outfld('TAUKC13 ',taukc13,pcols,lchnk)
!        call outfld('TAUKC14 ',taukc14,pcols,lchnk)

!      write(*,*) 'pmxsub 9'

#ifdef AEROCOM  ! AEROCOM***********AEROCOM**************AEROCOM***************below   

!        call outfld('CAM5COND',cam5cond,pcols,lchnk)
!        call outfld('CAM5AQ  ',cam5aq  ,pcols,lchnk)
!        call outfld('CAM5OC  ',cam5oc  ,pcols,lchnk)
!        call outfld('CAM5BC  ',cam5bc  ,pcols,lchnk)
!        call outfld('BEKC5   ',bekc5   ,pcols,lchnk)
!        call outfld('CAM4COND',cam4cond,pcols,lchnk)
!        call outfld('CAM4AQ  ',cam4aq  ,pcols,lchnk)
!        call outfld('CAM4OC  ',cam4oc  ,pcols,lchnk)
!        call outfld('CAM4BC  ',cam4bc  ,pcols,lchnk)
!        call outfld('BEKC4   ',bekc4   ,pcols,lchnk)

!       Initialize fields
        do icol=1,ncol
          daerh2o(icol)=0.0_r8
          vaercols(icol)=0.0_r8
          vaercoll(icol)=0.0_r8
          aaercols(icol)=0.0_r8
          aaercoll(icol)=0.0_r8
          do i=0,nmodes
            dload(icol,i)=0.0_r8
          enddo
        enddo

!       AEROCOM-specific interpolation routine for internally mixed aerosols:

!     BC(ax) mode (hydrophobic, so no rhum needed here):
        call intaeropt0(lchnk, ncol, Nnatk, &
!old          bext1, babs1, bext2, babs2, bebg1, babg1, bebg2, babg2, & 
!old          bebc1, babc1, bebc2, babc2, beoc1, baoc1, beoc2, baoc2, &
!old          bes41, bas41, bes42, bas42, bebglt1, bebggt1, bebclt1, &
!old          bebcgt1, beoclt1, beocgt1, bes4lt1, bes4gt1)  
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebglt1, bebggt1, bebclt1, bebcgt1,                         &
           beoclt1, beocgt1, bes4lt1, bes4gt1,                         &
           backsc550, babg550, babc550, baoc550, basu550)

!     SO4(Ait), BC(Ait) and OC(Ait) modes:
      mplus10=0
        call intaeropt1to3(lchnk, ncol, rhum, mplus10, Nnatk, Cam, &
!old          bext1, babs1, bext2, babs2, bebg1, babg1, bebg2, babg2, & 
!old          bebc1, babc1, bebc2, babc2, beoc1, baoc1, beoc2, baoc2, &
!old          bes41, bas41, bes42, bas42, bebglt1, bebggt1, bebclt1, &
!old          bebcgt1, beoclt1, beocgt1, bes4lt1, bes4gt1)  
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebglt1, bebggt1, bebclt1, bebcgt1,                         &
           beoclt1, beocgt1, bes4lt1, bes4gt1,                         &
           backsc550, babg550, babc550, baoc550, basu550)

!     BC&OC(Ait) mode:    ------ fcm invalid here (=0). Using faitbc instead
      mplus10=0
        call intaeropt4(lchnk, ncol, rhum, mplus10, Nnatk, Cam, faitbc, faqm4, &
!old          bext1, babs1, bext2, babs2, bebg1, babg1, bebg2, babg2, & 
!old          bebc1, babc1, bebc2, babc2, beoc1, baoc1, beoc2, baoc2, &
!old          bes41, bas41, bes42, bas42, bebglt1, bebggt1, bebclt1, &
!old          bebcgt1, beoclt1, beocgt1, bes4lt1, bes4gt1)  
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebglt1, bebggt1, bebclt1, bebcgt1,                         &
           beoclt1, beocgt1, bes4lt1, bes4gt1,                         & 
           backsc550, babg550, babc550, baoc550, basu550)

!     SO4(Ait75) (5), Mineral (6-7) and Sea-salt (8-10) modes:
        call intaeropt5to10(lchnk, ncol, rhum, Nnatk, Cam, fcm, fbcm, faqm, &
!old          bext1, babs1, bext2, babs2, bebg1, babg1, bebg2, babg2, & 
!old          bebc1, babc1, bebc2, babc2, beoc1, baoc1, beoc2, baoc2, &
!old          bes41, bas41, bes42, bas42, bebglt1, bebggt1, bebclt1, &
!old          bebcgt1, beoclt1, beocgt1, bes4lt1, bes4gt1)  
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebglt1, bebggt1, bebclt1, bebcgt1,                         &
           beoclt1, beocgt1, bes4lt1, bes4gt1,                         &
           backsc550, babg550, babc550, baoc550, basu550)

        do k=1,pver
         do icol=1,ncol
!old
!          bext1tot(icol,k)=0.0_r8
!          babs1tot(icol,k)=0.0_r8
!          bext2tot(icol,k)=0.0_r8
!          babs2tot(icol,k)=0.0_r8
!          bebg1tot(icol,k)=0.0_r8
!          babg1tot(icol,k)=0.0_r8
!          bebg2tot(icol,k)=0.0_r8
!          babg2tot(icol,k)=0.0_r8
!          bebc1tot(icol,k)=0.0_r8
!          babc1tot(icol,k)=0.0_r8
!          bebc2tot(icol,k)=0.0_r8
!          babc2tot(icol,k)=0.0_r8
!          beoc1tot(icol,k)=0.0_r8
!          baoc1tot(icol,k)=0.0_r8
!          beoc2tot(icol,k)=0.0_r8
!          baoc2tot(icol,k)=0.0_r8
!          bes41tot(icol,k)=0.0_r8
!          bas41tot(icol,k)=0.0_r8
!          bes42tot(icol,k)=0.0_r8
!          bas42tot(icol,k)=0.0_r8
!old
         bebglt1t(icol,k)=0.0_r8
         bebggt1t(icol,k)=0.0_r8
         bebclt1t(icol,k)=0.0_r8
         bebcgt1t(icol,k)=0.0_r8
         beoclt1t(icol,k)=0.0_r8
         beocgt1t(icol,k)=0.0_r8
         bes4lt1t(icol,k)=0.0_r8
         bes4gt1t(icol,k)=0.0_r8 
         bedustlt1(icol,k)=0.0_r8
         bedustgt1(icol,k)=0.0_r8
         besslt1(icol,k)=0.0_r8
         bessgt1(icol,k)=0.0_r8

         bext440tot(icol,k)=0.0_r8 
         babs440tot(icol,k)=0.0_r8 
         bext500tot(icol,k)=0.0_r8 
         babs500tot(icol,k)=0.0_r8 
         bext550tot(icol,k)=0.0_r8 
         babs550tot(icol,k)=0.0_r8 
         bext670tot(icol,k)=0.0_r8 
         babs670tot(icol,k)=0.0_r8 
         bext870tot(icol,k)=0.0_r8 
         babs870tot(icol,k)=0.0_r8 

         backsc550tot(icol,k)=0.0_r8 
 
         bebg440tot(icol,k)=0.0_r8 
!         babg440tot(icol,k)=0.0_r8 
         bebg500tot(icol,k)=0.0_r8 
!         babg500tot(icol,k)=0.0_r8 
         bebg550tot(icol,k)=0.0_r8 
         babg550tot(icol,k)=0.0_r8 
         bebg670tot(icol,k)=0.0_r8 
!         babg670tot(icol,k)=0.0_r8 
         bebg870tot(icol,k)=0.0_r8 
!         babg870tot(icol,k)=0.0_r8 

         bebc440tot(icol,k)=0.0_r8 
!         babc440tot(icol,k)=0.0_r8 
         bebc500tot(icol,k)=0.0_r8 
!         babc500tot(icol,k)=0.0_r8 
         bebc550tot(icol,k)=0.0_r8 
         babc550tot(icol,k)=0.0_r8 
         bebc670tot(icol,k)=0.0_r8 
!         babc670tot(icol,k)=0.0_r8 
         bebc870tot(icol,k)=0.0_r8 
!         babc870tot(icol,k)=0.0_r8 

         beoc440tot(icol,k)=0.0_r8 
!         baoc440tot(icol,k)=0.0_r8 
         beoc500tot(icol,k)=0.0_r8 
!         baoc500tot(icol,k)=0.0_r8 
         beoc550tot(icol,k)=0.0_r8 
         baoc550tot(icol,k)=0.0_r8 
         beoc670tot(icol,k)=0.0_r8 
!         baoc670tot(icol,k)=0.0_r8 
         beoc870tot(icol,k)=0.0_r8 
!         baoc870tot(icol,k)=0.0_r8 

         besu440tot(icol,k)=0.0_r8 
!         basu440tot(icol,k)=0.0_r8 
         besu500tot(icol,k)=0.0_r8 
!         basu500tot(icol,k)=0.0_r8 
         besu550tot(icol,k)=0.0_r8 
         basu550tot(icol,k)=0.0_r8 
         besu670tot(icol,k)=0.0_r8 
!         basu670tot(icol,k)=0.0_r8 
         besu870tot(icol,k)=0.0_r8 
!         basu870tot(icol,k)=0.0_r8 

         enddo
        enddo

        do i=0,nbmodes
         do k=1,pver
          do icol=1,ncol
!g           if(Nnatk(icol,k,i)>0.0_r8) then      
!old      total internal extinction and absorption for 0.55 um (1) and 0.865 um (2)
!old       bext1tot(icol,k)=bext1tot(icol,k)+Nnatk(icol,k,i)*bext1(icol,k,i)
!old       babs1tot(icol,k)=babs1tot(icol,k)+Nnatk(icol,k,i)*babs1(icol,k,i)
!old       bext2tot(icol,k)=bext2tot(icol,k)+Nnatk(icol,k,i)*bext2(icol,k,i)
!old       babs2tot(icol,k)=babs2tot(icol,k)+Nnatk(icol,k,i)*babs2(icol,k,i)
!      total internal extinction and absorption for 0.44, 0.50, 0.55, 0.68 and 0.87 um
       bext440tot(icol,k)=bext440tot(icol,k)+Nnatk(icol,k,i)*bext440(icol,k,i)
       babs440tot(icol,k)=babs440tot(icol,k)+Nnatk(icol,k,i)*babs440(icol,k,i)
       bext500tot(icol,k)=bext500tot(icol,k)+Nnatk(icol,k,i)*bext500(icol,k,i)
       babs500tot(icol,k)=babs500tot(icol,k)+Nnatk(icol,k,i)*babs500(icol,k,i)
       bext550tot(icol,k)=bext550tot(icol,k)+Nnatk(icol,k,i)*bext550(icol,k,i)
       babs550tot(icol,k)=babs550tot(icol,k)+Nnatk(icol,k,i)*babs550(icol,k,i)
       bext670tot(icol,k)=bext670tot(icol,k)+Nnatk(icol,k,i)*bext670(icol,k,i)
       babs670tot(icol,k)=babs670tot(icol,k)+Nnatk(icol,k,i)*babs670(icol,k,i)
       bext870tot(icol,k)=bext870tot(icol,k)+Nnatk(icol,k,i)*bext870(icol,k,i)
       babs870tot(icol,k)=babs870tot(icol,k)+Nnatk(icol,k,i)*babs870(icol,k,i)
       backsc550tot(icol,k)=backsc550tot(icol,k)+Nnatk(icol,k,i)*backsc550(icol,k,i)
!old      extinction and absorption for 0.55 um (1) and 0.865 um (2)
!old      for the whole background aerosol (icluding SO4,BC&OC for modes 0-5)
!old       bebg1tot(icol,k)=bebg1tot(icol,k)+Nnatk(icol,k,i)*bebg1(icol,k,i)
!old       babg1tot(icol,k)=babg1tot(icol,k)+Nnatk(icol,k,i)*babg1(icol,k,i)
!old       bebg2tot(icol,k)=bebg2tot(icol,k)+Nnatk(icol,k,i)*bebg2(icol,k,i)
!old       babg2tot(icol,k)=babg2tot(icol,k)+Nnatk(icol,k,i)*babg2(icol,k,i)
!      extinction and absorption for 0.44, 0.50, 0.55 (no abs), 0.68 and 0.87 um
!      for the whole background aerosol (icluding SO4,BC&OC for modes 0-5)
       bebg440tot(icol,k)=bebg440tot(icol,k)+Nnatk(icol,k,i)*bebg440(icol,k,i)
!       babg440tot(icol,k)=babg440tot(icol,k)+Nnatk(icol,k,i)*babg440(icol,k,i)
       bebg500tot(icol,k)=bebg500tot(icol,k)+Nnatk(icol,k,i)*bebg500(icol,k,i)
!       babg500tot(icol,k)=babg500tot(icol,k)+Nnatk(icol,k,i)*babg500(icol,k,i)
       bebg550tot(icol,k)=bebg550tot(icol,k)+Nnatk(icol,k,i)*bebg550(icol,k,i)
       babg550tot(icol,k)=babg550tot(icol,k)+Nnatk(icol,k,i)*babg550(icol,k,i)
       bebg670tot(icol,k)=bebg670tot(icol,k)+Nnatk(icol,k,i)*bebg670(icol,k,i)
!       babg670tot(icol,k)=babg670tot(icol,k)+Nnatk(icol,k,i)*babg670(icol,k,i)
       bebg870tot(icol,k)=bebg870tot(icol,k)+Nnatk(icol,k,i)*bebg870(icol,k,i)
!       babg870tot(icol,k)=babg870tot(icol,k)+Nnatk(icol,k,i)*babg870(icol,k,i)
!old      extinction and absorption for 0.55 um (1) and 0.865 um (2)
!old      for each added (internally mixed through Aq./cond./coag.) component (SO4,BC,OC)
!old      Condensed/coagulated SO4 on all modes 1-10, and wet-phase SO4 on modes 4-10:
!old         bes41tot(icol,k)=bes41tot(icol,k)+Nnatk(icol,k,i)*bes41(icol,k,i)
!old         bas41tot(icol,k)=bas41tot(icol,k)+Nnatk(icol,k,i)*bas41(icol,k,i)
!old         bes42tot(icol,k)=bes42tot(icol,k)+Nnatk(icol,k,i)*bes42(icol,k,i)
!old         bas42tot(icol,k)=bas42tot(icol,k)+Nnatk(icol,k,i)*bas42(icol,k,i)
!      extinction and absorption for 0.44, 0.50, 0.55 (no abs), 0.68 and 0.87 um
!      for each added (internally mixed through Aq./cond./coag.) component (SO4,BC,OC).
!      Condensed/coagulated SO4 on all modes 1-10, and wet-phase SO4 on modes 4-10:
       besu440tot(icol,k)=besu440tot(icol,k)+Nnatk(icol,k,i)*besu440(icol,k,i)
!       basu440tot(icol,k)=basu440tot(icol,k)+Nnatk(icol,k,i)*basu440(icol,k,i)
       besu500tot(icol,k)=besu500tot(icol,k)+Nnatk(icol,k,i)*besu500(icol,k,i)
!       basu500tot(icol,k)=basu500tot(icol,k)+Nnatk(icol,k,i)*basu500(icol,k,i)
       besu550tot(icol,k)=besu550tot(icol,k)+Nnatk(icol,k,i)*besu550(icol,k,i)
       basu550tot(icol,k)=basu550tot(icol,k)+Nnatk(icol,k,i)*basu550(icol,k,i)
       besu670tot(icol,k)=besu670tot(icol,k)+Nnatk(icol,k,i)*besu670(icol,k,i)
!       basu670tot(icol,k)=basu670tot(icol,k)+Nnatk(icol,k,i)*basu670(icol,k,i)
       besu870tot(icol,k)=besu870tot(icol,k)+Nnatk(icol,k,i)*besu870(icol,k,i)
!       basu870tot(icol,k)=basu870tot(icol,k)+Nnatk(icol,k,i)*basu870(icol,k,i)
!
!      Coagulated BC and OC on modes 5-10:
       if(i>=5) then
!old         bebc1tot(icol,k)=bebc1tot(icol,k)+Nnatk(icol,k,i)*bebc1(icol,k,i)
!old         babc1tot(icol,k)=babc1tot(icol,k)+Nnatk(icol,k,i)*babc1(icol,k,i)
!old         bebc2tot(icol,k)=bebc2tot(icol,k)+Nnatk(icol,k,i)*bebc2(icol,k,i)
!old         babc2tot(icol,k)=babc2tot(icol,k)+Nnatk(icol,k,i)*babc2(icol,k,i)
!old         beoc1tot(icol,k)=beoc1tot(icol,k)+Nnatk(icol,k,i)*beoc1(icol,k,i)
!old         baoc1tot(icol,k)=baoc1tot(icol,k)+Nnatk(icol,k,i)*baoc1(icol,k,i)
!old         beoc2tot(icol,k)=beoc2tot(icol,k)+Nnatk(icol,k,i)*beoc2(icol,k,i)
!old         baoc2tot(icol,k)=baoc2tot(icol,k)+Nnatk(icol,k,i)*baoc2(icol,k,i)
       bebc440tot(icol,k)=bebc440tot(icol,k)+Nnatk(icol,k,i)*bebc440(icol,k,i)
!       babc440tot(icol,k)=babc440tot(icol,k)+Nnatk(icol,k,i)*babc440(icol,k,i)
       bebc500tot(icol,k)=bebc500tot(icol,k)+Nnatk(icol,k,i)*bebc500(icol,k,i)
!       babc500tot(icol,k)=babc500tot(icol,k)+Nnatk(icol,k,i)*babc500(icol,k,i)
       bebc550tot(icol,k)=bebc550tot(icol,k)+Nnatk(icol,k,i)*bebc550(icol,k,i)
       babc550tot(icol,k)=babc550tot(icol,k)+Nnatk(icol,k,i)*babc550(icol,k,i)
       bebc670tot(icol,k)=bebc670tot(icol,k)+Nnatk(icol,k,i)*bebc670(icol,k,i)
!       babc670tot(icol,k)=babc670tot(icol,k)+Nnatk(icol,k,i)*babc670(icol,k,i)
       bebc870tot(icol,k)=bebc870tot(icol,k)+Nnatk(icol,k,i)*bebc870(icol,k,i)
!       babc870tot(icol,k)=babc870tot(icol,k)+Nnatk(icol,k,i)*babc870(icol,k,i)
       beoc440tot(icol,k)=beoc440tot(icol,k)+Nnatk(icol,k,i)*beoc440(icol,k,i)
!       baoc440tot(icol,k)=baoc440tot(icol,k)+Nnatk(icol,k,i)*baoc440(icol,k,i)
       beoc500tot(icol,k)=beoc500tot(icol,k)+Nnatk(icol,k,i)*beoc500(icol,k,i)
!       baoc500tot(icol,k)=baoc500tot(icol,k)+Nnatk(icol,k,i)*baoc500(icol,k,i)
       beoc550tot(icol,k)=beoc550tot(icol,k)+Nnatk(icol,k,i)*beoc550(icol,k,i)
       baoc550tot(icol,k)=baoc550tot(icol,k)+Nnatk(icol,k,i)*baoc550(icol,k,i)
       beoc670tot(icol,k)=beoc670tot(icol,k)+Nnatk(icol,k,i)*beoc670(icol,k,i)
!       baoc670tot(icol,k)=baoc670tot(icol,k)+Nnatk(icol,k,i)*baoc670(icol,k,i)
       beoc870tot(icol,k)=beoc870tot(icol,k)+Nnatk(icol,k,i)*beoc870(icol,k,i)
!       baoc870tot(icol,k)=baoc870tot(icol,k)+Nnatk(icol,k,i)*baoc870(icol,k,i)
         bebglt1t(icol,k)=bebglt1t(icol,k) &
                       +Nnatk(icol,k,i)*bebglt1(icol,k,i)
         bebggt1t(icol,k)=bebggt1t(icol,k) &
                        +Nnatk(icol,k,i)*bebggt1(icol,k,i)
       endif
       if(i==6.or.i==7) then
         bedustlt1(icol,k)=bedustlt1(icol,k) &
                        +Nnatk(icol,k,i)*bebglt1(icol,k,i)
         bedustgt1(icol,k)=bedustgt1(icol,k) &
                        +Nnatk(icol,k,i)*bebggt1(icol,k,i)
       elseif(i>=8.and.i<=10) then
         besslt1(icol,k)=besslt1(icol,k) &
                        +Nnatk(icol,k,i)*bebglt1(icol,k,i)
         bessgt1(icol,k)=bessgt1(icol,k) &
                        +Nnatk(icol,k,i)*bebggt1(icol,k,i)
       endif
!      Condensed/coagulated SO4 on all modes 1-10, and wet-phase SO4 on modes 4-10:
         bes4lt1t(icol,k)=bes4lt1t(icol,k) &
                        +Nnatk(icol,k,i)*bes4lt1(icol,k,i)
         bes4gt1t(icol,k)=bes4gt1t(icol,k) &
                        +Nnatk(icol,k,i)*bes4gt1(icol,k,i)
!      Coagulated BC and OC on modes 5-10:
       if(i>=5) then
         bebclt1t(icol,k)=bebclt1t(icol,k) &
                        +Nnatk(icol,k,i)*bebclt1(icol,k,i)
         bebcgt1t(icol,k)=bebcgt1t(icol,k) &
                        +Nnatk(icol,k,i)*bebcgt1(icol,k,i)
         beoclt1t(icol,k)=beoclt1t(icol,k) &
                        +Nnatk(icol,k,i)*beoclt1(icol,k,i)
         beocgt1t(icol,k)=beocgt1t(icol,k) &
                        +Nnatk(icol,k,i)*beocgt1(icol,k,i)
       endif
!g           endif
          end do   ! icol
         enddo     ! k
        enddo      ! i

!      extinction/absorptions (km-1) for each background component 
!      in the internal mixture are
        do k=1,pver
          do icol=1,ncol
!old            bint1mi(icol,k)=Nnatk(icol,k,6)*bebg1(icol,k,6) &
!old                         +Nnatk(icol,k,7)*bebg1(icol,k,7)
!old-            baintmi(icol,k)=Nnatk(icol,k,6)*babg1(icol,k,6) &
!old-                         +Nnatk(icol,k,7)*babg1(icol,k,7)
!old            bint2mi(icol,k)=Nnatk(icol,k,6)*bebg2(icol,k,6) &
!old                         +Nnatk(icol,k,7)*bebg2(icol,k,7)
!old added!?                         +Nnatk(icol,k,7)*babg1(icol,k,7)
!old            bint1ss(icol,k)=Nnatk(icol,k,8)*bebg1(icol,k,8) &
!old                         +Nnatk(icol,k,9)*bebg1(icol,k,9)   &
!old                         +Nnatk(icol,k,10)*bebg1(icol,k,10)
!old-            baintss(icol,k)=Nnatk(icol,k,8)*babg1(icol,k,8) &
!old                         +Nnatk(icol,k,9)*babg1(icol,k,9)   &
!old                         +Nnatk(icol,k,10)*babg1(icol,k,10)
!old            bint2ss(icol,k)=Nnatk(icol,k,8)*bebg2(icol,k,8) &
!old                         +Nnatk(icol,k,9)*bebg2(icol,k,9)   &
!old                         +Nnatk(icol,k,10)*bebg2(icol,k,10)
            bint440du(icol,k)=Nnatk(icol,k,6)*bebg440(icol,k,6) &
                             +Nnatk(icol,k,7)*bebg440(icol,k,7)
            bint500du(icol,k)=Nnatk(icol,k,6)*bebg500(icol,k,6) &
                             +Nnatk(icol,k,7)*bebg500(icol,k,7)
            bint550du(icol,k)=Nnatk(icol,k,6)*bebg550(icol,k,6) &
                             +Nnatk(icol,k,7)*bebg550(icol,k,7)
            bint670du(icol,k)=Nnatk(icol,k,6)*bebg670(icol,k,6) &
                             +Nnatk(icol,k,7)*bebg670(icol,k,7)
            bint870du(icol,k)=Nnatk(icol,k,6)*bebg870(icol,k,6) &
                             +Nnatk(icol,k,7)*bebg870(icol,k,7)
            bint440ss(icol,k)=Nnatk(icol,k,8)*bebg440(icol,k,8) &
                             +Nnatk(icol,k,9)*bebg440(icol,k,9) &
                            +Nnatk(icol,k,10)*bebg440(icol,k,10)
            bint500ss(icol,k)=Nnatk(icol,k,8)*bebg500(icol,k,8) &
                             +Nnatk(icol,k,9)*bebg500(icol,k,9) &
                            +Nnatk(icol,k,10)*bebg500(icol,k,10)
            bint550ss(icol,k)=Nnatk(icol,k,8)*bebg550(icol,k,8) &
                             +Nnatk(icol,k,9)*bebg550(icol,k,9) &
                            +Nnatk(icol,k,10)*bebg550(icol,k,10)
            bint670ss(icol,k)=Nnatk(icol,k,8)*bebg670(icol,k,8) &
                             +Nnatk(icol,k,9)*bebg670(icol,k,9) &
                            +Nnatk(icol,k,10)*bebg670(icol,k,10)
            bint870ss(icol,k)=Nnatk(icol,k,8)*bebg870(icol,k,8) &
                             +Nnatk(icol,k,9)*bebg870(icol,k,9) &
                            +Nnatk(icol,k,10)*bebg870(icol,k,10)
            baint550du(icol,k)=Nnatk(icol,k,6)*babg550(icol,k,6) &
                             +Nnatk(icol,k,7)*babg550(icol,k,7)
            baint550ss(icol,k)=Nnatk(icol,k,8)*babg550(icol,k,8) &
                             +Nnatk(icol,k,9)*babg550(icol,k,9) &
                            +Nnatk(icol,k,10)*babg550(icol,k,10)
          end do
        enddo

!     then to the externally mixed SO4(n), BC(n) and OC(n) modes:
      mplus10=1
        call intaeropt1to3(lchnk, ncol, rhum, mplus10, Nnatk, camnull, &
!old          bext1n, babs1n, bext2n, babs2n, bebg1n, babg1n, bebg2n, babg2n, & 
!old          bebc1n, babc1n, bebc2n, babc2n, beoc1n, baoc1n, beoc2n, baoc2n, &
!old          bes41n, bas41n, bes42n, bas42n, bebglt1n, bebggt1n, bebclt1n, &
!old          bebcgt1n, beoclt1n, beocgt1n, bes4lt1n, bes4gt1n)  
           bext440n, bext500n, bext550n, bext670n, bext870n,                &
           bebg440n, bebg500n, bebg550n, bebg670n, bebg870n,                &
           bebc440n, bebc500n, bebc550n, bebc670n, bebc870n,                &
           beoc440n, beoc500n, beoc550n, beoc670n, beoc870n,                &
           besu440n, besu500n, besu550n, besu670n, besu870n,                &
           babs440n, babs500n, babs550n, babs670n, babs870n,                &
           bebglt1n, bebggt1n, bebclt1n, bebcgt1n,                         &
           beoclt1n, beocgt1n, bes4lt1n, bes4gt1n,                         &
           backsc550n, babg550n, babc550n, baoc550n, basu550n)
!     and finally the BC&OC(n) mode:
      mplus10=1
        call intaeropt4(lchnk, ncol, rhum, mplus10, Nnatk, camnull, fnbc, faqm4, &
!old          bext1n, babs1n, bext2n, babs2n, bebg1n, babg1n, bebg2n, babg2n, & 
!old          bebc1n, babc1n, bebc2n, babc2n, beoc1n, baoc1n, beoc2n, baoc2n, &
!old          bes41n, bas41n, bes42n, bas42n, bebglt1n, bebggt1n, bebclt1n, &
!old          bebcgt1n, beoclt1n, beocgt1n, bes4lt1n, bes4gt1n)  
           bext440n, bext500n, bext550n, bext670n, bext870n,                &
           bebg440n, bebg500n, bebg550n, bebg670n, bebg870n,                &
           bebc440n, bebc500n, bebc550n, bebc670n, bebc870n,                &
           beoc440n, beoc500n, beoc550n, beoc670n, beoc870n,                &
           besu440n, besu500n, besu550n, besu670n, besu870n,                &
           babs440n, babs500n, babs550n, babs670n, babs870n,                &
           bebglt1n, bebggt1n, bebclt1n, bebcgt1n,                         &
           beoclt1n, beocgt1n, bes4lt1n, bes4gt1n,                         &
           backsc550n, babg550n, babc550n, baoc550n, basu550n)

      do i=11,14
       do k=1,pver
        do icol=1,ncol
!old          be1x(icol,k,i)=bext1n(icol,k,i-10)
!old          ba1x(icol,k,i)=babs1n(icol,k,i-10)
!old          be2x(icol,k,i)=bext2n(icol,k,i-10)
!old          ba2x(icol,k,i)=babs2n(icol,k,i-10)
          be440x(icol,k,i)=bext440n(icol,k,i-10)
          ba440x(icol,k,i)=babs440n(icol,k,i-10)
          be500x(icol,k,i)=bext500n(icol,k,i-10)
          ba500x(icol,k,i)=babs500n(icol,k,i-10)
          be550x(icol,k,i)=bext550n(icol,k,i-10)
          ba550x(icol,k,i)=babs550n(icol,k,i-10)
          be670x(icol,k,i)=bext670n(icol,k,i-10)
          ba670x(icol,k,i)=babs670n(icol,k,i-10)
          be870x(icol,k,i)=bext870n(icol,k,i-10)
          ba870x(icol,k,i)=babs870n(icol,k,i-10)
          belt1x(icol,k,i)=bebglt1n(icol,k,i-10)
          begt1x(icol,k,i)=bebggt1n(icol,k,i-10)
          backsc550x(icol,k,i)=backsc550n(icol,k,i-10)
        end do
       enddo
      enddo

!     The externally modes' contribution to extinction and absorption:
         do k=1,pver
          do icol=1,ncol
!SO4
!old            bes41xt(icol,k) =Nnatk(icol,k,11)*be1x(icol,k,11)
!old            bas41xt(icol,k) =Nnatk(icol,k,11)*ba1x(icol,k,11)
!old            bes42xt(icol,k) =Nnatk(icol,k,11)*be2x(icol,k,11)
!old            bas42xt(icol,k) =Nnatk(icol,k,11)*ba2x(icol,k,11)
            besu440xt(icol,k) =Nnatk(icol,k,11)*be440x(icol,k,11)
            basu440xt(icol,k) =Nnatk(icol,k,11)*ba440x(icol,k,11)
            besu500xt(icol,k) =Nnatk(icol,k,11)*be500x(icol,k,11)
            basu500xt(icol,k) =Nnatk(icol,k,11)*ba500x(icol,k,11)
            besu550xt(icol,k) =Nnatk(icol,k,11)*be550x(icol,k,11)
            basu550xt(icol,k) =Nnatk(icol,k,11)*ba550x(icol,k,11)
            besu670xt(icol,k) =Nnatk(icol,k,11)*be670x(icol,k,11)
            basu670xt(icol,k) =Nnatk(icol,k,11)*ba670x(icol,k,11)
            besu870xt(icol,k) =Nnatk(icol,k,11)*be870x(icol,k,11)
            basu870xt(icol,k) =Nnatk(icol,k,11)*ba870x(icol,k,11)
            bs4lt1xt(icol,k)=Nnatk(icol,k,11)*belt1x(icol,k,11)
            bs4gt1xt(icol,k)=Nnatk(icol,k,11)*begt1x(icol,k,11)
!BC
            vnbc = fnbc(icol,k)/(fnbc(icol,k) &
                           +(1.0_r8-fnbc(icol,k))*rhopart(l_bc_ni)/rhopart(l_om_ni))
!old            bebc1xt(icol,k) =Nnatk(icol,k,12)*be1x(icol,k,12)  &
!old                       +vnbc*Nnatk(icol,k,14)*be1x(icol,k,14)
!old            babc1xt(icol,k) =Nnatk(icol,k,12)*ba1x(icol,k,12)  &
!old                       +vnbc*Nnatk(icol,k,14)*ba1x(icol,k,14)
!old            bebc2xt(icol,k) =Nnatk(icol,k,12)*be2x(icol,k,12)  &
!old                       +vnbc*Nnatk(icol,k,14)*be2x(icol,k,14)
!old            babc2xt(icol,k) =Nnatk(icol,k,12)*ba2x(icol,k,12)  &
!old                       +vnbc*Nnatk(icol,k,14)*ba2x(icol,k,14)
            bebc440xt(icol,k) =Nnatk(icol,k,12)*be440x(icol,k,12)  &
                         +vnbc*Nnatk(icol,k,14)*be440x(icol,k,14)
            babc440xt(icol,k) =Nnatk(icol,k,12)*ba440x(icol,k,12)  &
                         +vnbc*Nnatk(icol,k,14)*ba440x(icol,k,14)
            bebc500xt(icol,k) =Nnatk(icol,k,12)*be500x(icol,k,12)  &
                         +vnbc*Nnatk(icol,k,14)*be500x(icol,k,14)
            babc500xt(icol,k) =Nnatk(icol,k,12)*ba500x(icol,k,12)  &
                         +vnbc*Nnatk(icol,k,14)*ba500x(icol,k,14)
            bebc550xt(icol,k) =Nnatk(icol,k,12)*be550x(icol,k,12)  &
                         +vnbc*Nnatk(icol,k,14)*be550x(icol,k,14)
            babc550xt(icol,k) =Nnatk(icol,k,12)*ba550x(icol,k,12)  &
                         +vnbc*Nnatk(icol,k,14)*ba550x(icol,k,14)
            bebc670xt(icol,k) =Nnatk(icol,k,12)*be670x(icol,k,12)  &
                         +vnbc*Nnatk(icol,k,14)*be670x(icol,k,14)
            babc670xt(icol,k) =Nnatk(icol,k,12)*ba670x(icol,k,12)  &
                         +vnbc*Nnatk(icol,k,14)*ba670x(icol,k,14)
            bebc870xt(icol,k) =Nnatk(icol,k,12)*be870x(icol,k,12)  &
                         +vnbc*Nnatk(icol,k,14)*be870x(icol,k,14)
            babc870xt(icol,k) =Nnatk(icol,k,12)*ba870x(icol,k,12)  &
                         +vnbc*Nnatk(icol,k,14)*ba870x(icol,k,14)
            bbclt1xt(icol,k)=Nnatk(icol,k,12)*belt1x(icol,k,12) &
                       +vnbc*Nnatk(icol,k,14)*belt1x(icol,k,14)
            bbcgt1xt(icol,k)=Nnatk(icol,k,12)*begt1x(icol,k,12) &
                       +vnbc*Nnatk(icol,k,14)*begt1x(icol,k,14)
!OC
!old            beoc1xt(icol,k)=Nnatk(icol,k,13)*be1x(icol,k,13) &
!old                +(1.0_r8-vnbc)*Nnatk(icol,k,14)*be1x(icol,k,14)
!old            baoc1xt(icol,k)=Nnatk(icol,k,13)*ba1x(icol,k,13) &
!old                +(1.0_r8-vnbc)*Nnatk(icol,k,14)*ba1x(icol,k,14)
!old            beoc2xt(icol,k)=Nnatk(icol,k,13)*be2x(icol,k,13) &
!old                +(1.0_r8-vnbc)*Nnatk(icol,k,14)*be2x(icol,k,14)
!old            baoc2xt(icol,k)=Nnatk(icol,k,13)*ba2x(icol,k,13) &
!old                +(1.0_r8-vnbc)*Nnatk(icol,k,14)*ba2x(icol,k,14)
            beoc440xt(icol,k)=Nnatk(icol,k,13)*be440x(icol,k,13) &
                  +(1.0_r8-vnbc)*Nnatk(icol,k,14)*be440x(icol,k,14)
            baoc440xt(icol,k)=Nnatk(icol,k,13)*ba440x(icol,k,13) &
                  +(1.0_r8-vnbc)*Nnatk(icol,k,14)*ba440x(icol,k,14)
            beoc500xt(icol,k)=Nnatk(icol,k,13)*be500x(icol,k,13) &
                  +(1.0_r8-vnbc)*Nnatk(icol,k,14)*be500x(icol,k,14)
            baoc500xt(icol,k)=Nnatk(icol,k,13)*ba500x(icol,k,13) &
                  +(1.0_r8-vnbc)*Nnatk(icol,k,14)*ba500x(icol,k,14)
            beoc550xt(icol,k)=Nnatk(icol,k,13)*be550x(icol,k,13) &
                  +(1.0_r8-vnbc)*Nnatk(icol,k,14)*be550x(icol,k,14)
            baoc550xt(icol,k)=Nnatk(icol,k,13)*ba550x(icol,k,13) &
                  +(1.0_r8-vnbc)*Nnatk(icol,k,14)*ba550x(icol,k,14)
            beoc670xt(icol,k)=Nnatk(icol,k,13)*be670x(icol,k,13) &
                  +(1.0_r8-vnbc)*Nnatk(icol,k,14)*be670x(icol,k,14)
            baoc670xt(icol,k)=Nnatk(icol,k,13)*ba670x(icol,k,13) &
                  +(1.0_r8-vnbc)*Nnatk(icol,k,14)*ba670x(icol,k,14)
            beoc870xt(icol,k)=Nnatk(icol,k,13)*be870x(icol,k,13) &
                  +(1.0_r8-vnbc)*Nnatk(icol,k,14)*be870x(icol,k,14)
            baoc870xt(icol,k)=Nnatk(icol,k,13)*ba870x(icol,k,13) &
                  +(1.0_r8-vnbc)*Nnatk(icol,k,14)*ba870x(icol,k,14)
            boclt1xt(icol,k)=Nnatk(icol,k,13)*belt1x(icol,k,13) &
                 +(1.0_r8-vnbc)*Nnatk(icol,k,14)*belt1x(icol,k,14)
            bocgt1xt(icol,k)=Nnatk(icol,k,13)*begt1x(icol,k,13) &
                 +(1.0_r8-vnbc)*Nnatk(icol,k,14)*begt1x(icol,k,14)
!
            ec550_aer(icol,k)=bext550tot(icol,k)  &
              +Nnatk(icol,k,11)*be550x(icol,k,11) &
              +Nnatk(icol,k,12)*be550x(icol,k,12) &
              +Nnatk(icol,k,13)*be550x(icol,k,13) &
              +Nnatk(icol,k,14)*be550x(icol,k,14)
            ec550_aer(icol,k)=1.e-3_r8*ec550_aer(icol,k)
            abs550_aer(icol,k)=babs550tot(icol,k)  &
              +Nnatk(icol,k,11)*ba550x(icol,k,11) &
              +Nnatk(icol,k,12)*ba550x(icol,k,12) &
              +Nnatk(icol,k,13)*ba550x(icol,k,13) &
              +Nnatk(icol,k,14)*ba550x(icol,k,14)
            abs550_aer(icol,k)=1.e-3_r8*abs550_aer(icol,k)
            bs550_aer(icol,k)= backsc550tot(icol,k)   &
              +Nnatk(icol,k,11)*backsc550x(icol,k,11) &
              +Nnatk(icol,k,12)*backsc550x(icol,k,12) &
              +Nnatk(icol,k,13)*backsc550x(icol,k,13) &
              +Nnatk(icol,k,14)*backsc550x(icol,k,14)
            bs550_aer(icol,k)=1.e-3_r8*bs550_aer(icol,k)
!
          end do
         enddo
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!       collect AeroCom-fields for optical depth/absorption of each comp,
!old       3D and 2D, at 550 and 865 nm, for all r, r<1um and r>1um     
!       3D and 2D, at 440, 500, 550, 670 and 870 nm, for all d, d<1um and d>1um  
!        initialize 2d-fields
         do icol=1,ncol
           dod440(icol) = 0.0_r8
           abs440(icol) = 0.0_r8
           dod500(icol) = 0.0_r8
           abs500(icol) = 0.0_r8
           dod550(icol) = 0.0_r8
           abs550(icol) = 0.0_r8
           abs550alt(icol) = 0.0_r8
           dod670(icol) = 0.0_r8
           abs670(icol) = 0.0_r8
           dod870(icol) = 0.0_r8
           abs870(icol) = 0.0_r8

           abs550_ss(icol) = 0.0_r8
           abs550_dust(icol) = 0.0_r8
           abs550_so4(icol) = 0.0_r8
           abs550_bc(icol) = 0.0_r8
           abs550_pom(icol) = 0.0_r8
!
           dod440_ss(icol) = 0.0_r8
           dod440_dust(icol) = 0.0_r8
           dod440_so4(icol) = 0.0_r8
           dod440_bc(icol) = 0.0_r8
           dod440_pom(icol) = 0.0_r8
           dod500_ss(icol) = 0.0_r8
           dod500_dust(icol) = 0.0_r8
           dod500_so4(icol) = 0.0_r8
           dod500_bc(icol) = 0.0_r8
           dod500_pom(icol) = 0.0_r8
           dod550_ss(icol) = 0.0_r8
           dod550_dust(icol) = 0.0_r8
           dod550_so4(icol) = 0.0_r8
           dod550_bc(icol) = 0.0_r8
           dod550_pom(icol) = 0.0_r8
           dod670_ss(icol) = 0.0_r8
           dod670_dust(icol) = 0.0_r8
           dod670_so4(icol) = 0.0_r8
           dod670_bc(icol) = 0.0_r8
           dod670_pom(icol) = 0.0_r8
           dod870_ss(icol) = 0.0_r8
           dod870_dust(icol) = 0.0_r8
           dod870_so4(icol) = 0.0_r8
           dod870_bc(icol) = 0.0_r8
           dod870_pom(icol) = 0.0_r8
!old           dod865_ss(icol) = 0.0_r8
!old           dod865_dust(icol) = 0.0_r8
!old           dod865_so4(icol) = 0.0_r8
!old           dod865_bc(icol) = 0.0_r8
!old           dod865_pom(icol) = 0.0_r8
           dod550lt1_ss(icol) = 0.0_r8
           dod550gt1_ss(icol) = 0.0_r8
           dod550lt1_dust(icol) = 0.0_r8
           dod550gt1_dust(icol) = 0.0_r8
           dod550lt1_so4(icol) = 0.0_r8
           dod550gt1_so4(icol) = 0.0_r8
           dod550lt1_bc(icol) = 0.0_r8
           dod550gt1_bc(icol) = 0.0_r8
           dod550lt1_pom(icol) = 0.0_r8
           dod550gt1_pom(icol) = 0.0_r8
! OS Added abs components
!cakx           cloudfree(icol)=1.-cltot(icol)
           clearod550(icol) = 0.0_r8
           clearabs550(icol) = 0.0_r8
           clearabs550alt(icol) = 0.0_r8
!old           clearod865(icol) = 0.0_r8
           clearod870(icol) = 0.0_r8
          do k=1,pver
             abs5503d(icol,k) = 0.0_r8
             abs5503dalt(icol,k) = 0.0_r8
          enddo

         enddo
         do icol=1,ncol
          do k=1,pver
!          Layer thickness, unit km
           deltah=1.e-4_r8*(pint(icol,k+1)-pint(icol,k))/(rhoda(icol,k)*9.8_r8)
!          if(k==pver) write(*,*) 'icol, deltah(pmxsub)=', icol, deltah
!          3D optical depths for monthly averages
!SS
!old           dod5503d_ss(icol,k) = bint1ss(icol,k)*deltah
!old           abs5503d_ss(icol,k) = baintss(icol,k)*deltah
!old           dod8653d_ss(icol,k) = bint2ss(icol,k)*deltah
           dod4403d_ss(icol,k) = bint440ss(icol,k)*deltah
           dod5003d_ss(icol,k) = bint500ss(icol,k)*deltah
           dod5503d_ss(icol,k) = bint550ss(icol,k)*deltah
           abs5503d_ss(icol,k) = baint550ss(icol,k)*deltah
           dod6703d_ss(icol,k) = bint670ss(icol,k)*deltah
           dod8703d_ss(icol,k) = bint870ss(icol,k)*deltah
!DUST
!old           dod5503d_dust(icol,k) = bint1mi(icol,k)*deltah
!old           abs5503d_dust(icol,k) = baintmi(icol,k)*deltah
!old           dod8653d_dust(icol,k) = bint2mi(icol,k)*deltah
           dod4403d_dust(icol,k) = bint440du(icol,k)*deltah
           dod5003d_dust(icol,k) = bint500du(icol,k)*deltah
           dod5503d_dust(icol,k) = bint550du(icol,k)*deltah
           abs5503d_dust(icol,k) = baint550du(icol,k)*deltah
           dod6703d_dust(icol,k) = bint670du(icol,k)*deltah
           dod8703d_dust(icol,k) = bint870du(icol,k)*deltah
!SO4
!old           dod5503d_so4(icol,k) = (bes41tot(icol,k)+bes41xt(icol,k) &       ! condensate + n-mode (11)
!old                                  + Nnatk(icol,k,1)*bebg1(icol,k,1) &       ! background, SO4(Ait) mode (1)
!old                                  + Nnatk(icol,k,5)*bebg1(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
!old           abs5503d_so4(icol,k) = (bas41tot(icol,k)+bas41xt(icol,k) &       ! condensate + n-mode (11)
!old                                  + Nnatk(icol,k,1)*babg1(icol,k,1) &       ! background, SO4(Ait) mode (1)
!old                                  + Nnatk(icol,k,5)*babg1(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
!old           dod8653d_so4(icol,k) = (bes42tot(icol,k)+bes42xt(icol,k) &       ! condensate + n-mode (11)
!old                                  + Nnatk(icol,k,1)*bebg2(icol,k,1) &       ! background, SO4(Ait) mode (1)
!old                                  + Nnatk(icol,k,5)*bebg2(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
           dod4403d_so4(icol,k) = (besu440tot(icol,k)+besu440xt(icol,k) &     ! condensate + n-mode (11)
                                  + Nnatk(icol,k,1)*bebg440(icol,k,1) &       ! background, SO4(Ait) mode (1)
                                  + Nnatk(icol,k,5)*bebg440(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
           dod5003d_so4(icol,k) = (besu500tot(icol,k)+besu500xt(icol,k) &     ! condensate + n-mode (11)
                                  + Nnatk(icol,k,1)*bebg500(icol,k,1) &       ! background, SO4(Ait) mode (1)
                                  + Nnatk(icol,k,5)*bebg500(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
           dod5503d_so4(icol,k) = (besu550tot(icol,k)+besu550xt(icol,k) &     ! condensate + n-mode (11)
                                  + Nnatk(icol,k,1)*bebg550(icol,k,1) &       ! background, SO4(Ait) mode (1)
                                  + Nnatk(icol,k,5)*bebg550(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
           abs5503d_so4(icol,k) = (basu550tot(icol,k)+basu550xt(icol,k) &     ! condensate + n-mode (11)
                                  + Nnatk(icol,k,1)*babg550(icol,k,1) &       ! background, SO4(Ait) mode (1)
                                  + Nnatk(icol,k,5)*babg550(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
           dod6703d_so4(icol,k) = (besu670tot(icol,k)+besu670xt(icol,k) &     ! condensate + n-mode (11)
                                  + Nnatk(icol,k,1)*bebg670(icol,k,1) &       ! background, SO4(Ait) mode (1)
                                  + Nnatk(icol,k,5)*bebg670(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
           dod8703d_so4(icol,k) = (besu870tot(icol,k)+besu870xt(icol,k) &     ! condensate + n-mode (11)
                                  + Nnatk(icol,k,1)*bebg870(icol,k,1) &       ! background, SO4(Ait) mode (1)
                                  + Nnatk(icol,k,5)*bebg870(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
!BC
           vaitbc = faitbc(icol,k)/(faitbc(icol,k) &
                           +(1.0_r8-faitbc(icol,k))*rhopart(l_bc_ni)/rhopart(l_om_ni))
!old           dod5503d_bc(icol,k) = (bebc1tot(icol,k)+bebc1xt(icol,k)  &       ! coagulated + n-mode BC (12)
!old                                  + Nnatk(icol,k,2)*bebg1(icol,k,2) &       ! background, BC(Ait) mode (2)
!old                           + vaitbc*Nnatk(icol,k,4)*bebg1(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
!old                                  + Nnatk(icol,k,0)*bebg1(icol,k,0))*deltah ! background, BC(ax) mode (0)
!old           abs5503d_bc(icol,k) = (babc1tot(icol,k)+babc1xt(icol,k)  &       ! coagulated + n-mode BC (12)
!old                                  + Nnatk(icol,k,2)*babg1(icol,k,2) &       ! background, BC(Ait) mode (2)
!old                           + vaitbc*Nnatk(icol,k,4)*babg1(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
!old                                  + Nnatk(icol,k,0)*babg1(icol,k,0))*deltah ! background, BC(ax) mode (0)
!old           dod8653d_bc(icol,k) = (bebc2tot(icol,k)+bebc2xt(icol,k)  &       ! coagulated + n-mode BC (12)
!old                                  + Nnatk(icol,k,2)*bebg2(icol,k,2) &       ! background, BC(Ait) mode (2)
!old                           + vaitbc*Nnatk(icol,k,4)*bebg2(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
!old                                  + Nnatk(icol,k,0)*bebg2(icol,k,0))*deltah ! background, BC(ax) mode (0)
           dod4403d_bc(icol,k) = (bebc440tot(icol,k)+bebc440xt(icol,k)  &     ! coagulated + n-mode BC (12)
                                  + Nnatk(icol,k,2)*bebg440(icol,k,2) &       ! background, BC(Ait) mode (2)
                           + vaitbc*Nnatk(icol,k,4)*bebg440(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
                                  + Nnatk(icol,k,0)*bebg440(icol,k,0))*deltah ! background, BC(ax) mode (0)
           dod5003d_bc(icol,k) = (bebc500tot(icol,k)+bebc500xt(icol,k)  &     ! coagulated + n-mode BC (12)
                                  + Nnatk(icol,k,2)*bebg500(icol,k,2) &       ! background, BC(Ait) mode (2)
                           + vaitbc*Nnatk(icol,k,4)*bebg500(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
                                  + Nnatk(icol,k,0)*bebg500(icol,k,0))*deltah ! background, BC(ax) mode (0)
           dod5503d_bc(icol,k) = (bebc550tot(icol,k)+bebc550xt(icol,k)  &     ! coagulated + n-mode BC (12)
                                  + Nnatk(icol,k,2)*bebg550(icol,k,2) &       ! background, BC(Ait) mode (2)
                           + vaitbc*Nnatk(icol,k,4)*bebg550(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
                                  + Nnatk(icol,k,0)*bebg550(icol,k,0))*deltah ! background, BC(ax) mode (0)
           abs5503d_bc(icol,k) = (babc550tot(icol,k)+babc550xt(icol,k)  &     ! coagulated + n-mode BC (12)
                                  + Nnatk(icol,k,2)*babg550(icol,k,2) &       ! background, BC(Ait) mode (2)
                           + vaitbc*Nnatk(icol,k,4)*babg550(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
                                  + Nnatk(icol,k,0)*babg550(icol,k,0))*deltah ! background, BC(ax) mode (0)
           dod6703d_bc(icol,k) = (bebc670tot(icol,k)+bebc670xt(icol,k)  &     ! coagulated + n-mode BC (12)
                                  + Nnatk(icol,k,2)*bebg670(icol,k,2) &       ! background, BC(Ait) mode (2)
                           + vaitbc*Nnatk(icol,k,4)*bebg670(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
                                  + Nnatk(icol,k,0)*bebg670(icol,k,0))*deltah ! background, BC(ax) mode (0)
           dod8703d_bc(icol,k) = (bebc870tot(icol,k)+bebc870xt(icol,k)  &     ! coagulated + n-mode BC (12)
                                  + Nnatk(icol,k,2)*bebg870(icol,k,2) &       ! background, BC(Ait) mode (2)
                           + vaitbc*Nnatk(icol,k,4)*bebg870(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
                                  + Nnatk(icol,k,0)*bebg870(icol,k,0))*deltah ! background, BC(ax) mode (0)
!OC
!old           dod5503d_pom(icol,k) = (beoc1tot(icol,k)+beoc1xt(icol,k) &       ! coagulated + n-mode OC (13)
!old                                  + Nnatk(icol,k,3)*bebg1(icol,k,3) &       ! background, OC(Ait) mode (3)
!old                     + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*bebg1(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)
!old           abs5503d_pom(icol,k) = (baoc1tot(icol,k)+baoc1xt(icol,k) &       ! coagulated + n-mode OC (13)
!old                                  + Nnatk(icol,k,3)*babg1(icol,k,3) &       ! background, OC(Ait) mode (3)
!old                     + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*babg1(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)
!old           dod8653d_pom(icol,k) = (beoc2tot(icol,k)+beoc2xt(icol,k) &       ! coagulated + n-mode OC (13)
!old                                  + Nnatk(icol,k,3)*bebg2(icol,k,3) &       ! background, OC(Ait) mode (3)
!old                     + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*bebg2(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)
           dod4403d_pom(icol,k) = (beoc440tot(icol,k)+beoc440xt(icol,k) &     ! coagulated + n-mode OC (13)
                                  + Nnatk(icol,k,3)*bebg440(icol,k,3) &       ! background, OC(Ait) mode (3)
                     + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*bebg440(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)
           dod5003d_pom(icol,k) = (beoc500tot(icol,k)+beoc500xt(icol,k) &     ! coagulated + n-mode OC (13)
                                  + Nnatk(icol,k,3)*bebg500(icol,k,3) &       ! background, OC(Ait) mode (3)
                     + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*bebg500(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)
           dod5503d_pom(icol,k) = (beoc550tot(icol,k)+beoc550xt(icol,k) &     ! coagulated + n-mode OC (13)
                                  + Nnatk(icol,k,3)*bebg550(icol,k,3) &       ! background, OC(Ait) mode (3)
                     + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*bebg550(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)
           abs5503d_pom(icol,k) = (baoc550tot(icol,k)+baoc550xt(icol,k) &     ! coagulated + n-mode OC (13)
                                  + Nnatk(icol,k,3)*babg550(icol,k,3) &       ! background, OC(Ait) mode (3)
                     + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*babg550(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)
           dod6703d_pom(icol,k) = (beoc670tot(icol,k)+beoc670xt(icol,k) &     ! coagulated + n-mode OC (13)
                                  + Nnatk(icol,k,3)*bebg670(icol,k,3) &       ! background, OC(Ait) mode (3)
                     + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*bebg670(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)
           dod8703d_pom(icol,k) = (beoc870tot(icol,k)+beoc870xt(icol,k) &     ! coagulated + n-mode OC (13)
                                  + Nnatk(icol,k,3)*bebg870(icol,k,3) &       ! background, OC(Ait) mode (3)
                     + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*bebg870(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!          Total 3D optical depths/abs. for column integrations
           dod4403d(icol,k) = dod4403d_ss(icol,k)+dod4403d_dust(icol,k) &
                            +dod4403d_so4(icol,k)+dod4403d_bc(icol,k)    &
                            +dod4403d_pom(icol,k)                     
           dod5003d(icol,k) = dod5003d_ss(icol,k)+dod5003d_dust(icol,k) &
                            +dod5003d_so4(icol,k)+dod5003d_bc(icol,k)    &
                            +dod5003d_pom(icol,k)                     
           dod5503d(icol,k) = dod5503d_ss(icol,k)+dod5503d_dust(icol,k) &
                           +dod5503d_so4(icol,k)+dod5503d_bc(icol,k)    &
                           +dod5503d_pom(icol,k)                     
           dod6703d(icol,k) = dod6703d_ss(icol,k)+dod6703d_dust(icol,k) &
                            +dod6703d_so4(icol,k)+dod6703d_bc(icol,k)    &
                            +dod6703d_pom(icol,k)                     
           dod8703d(icol,k) = dod8703d_ss(icol,k)+dod8703d_dust(icol,k) &
                            +dod8703d_so4(icol,k)+dod8703d_bc(icol,k)    &
                            +dod8703d_pom(icol,k)                     
           abs5503d(icol,k) = abs5503d_ss(icol,k)+abs5503d_dust(icol,k) &
                           +abs5503d_so4(icol,k)+abs5503d_bc(icol,k)    &
                           +abs5503d_pom(icol,k)                     
!   (Note: Local abs550alt is up to 6% larger (annually averaged) in typical b.b.
!   regions, compared to abs550. This is most likely most correct, but should be checked!)
           do i=0,10
             abs5503dalt(icol,k) = abs5503dalt(icol,k)+Nnatk(icol,k,i)*babs550(icol,k,i)*deltah
           enddo
           do i=11,14
             abs5503dalt(icol,k) = abs5503dalt(icol,k)+Nnatk(icol,k,i)*babs550n(icol,k,i-10)*deltah
           enddo
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!old            optical depths for r<1um and r>1um
!            optical depths for d<1um and d>1um (r<0.5um and r>0.5um)
!SS
           dod5503dlt1_ss(icol,k) = besslt1(icol,k)*deltah
           dod5503dgt1_ss(icol,k) = bessgt1(icol,k)*deltah
!DUST
           dod5503dlt1_dust(icol,k) = bedustlt1(icol,k)*deltah
           dod5503dgt1_dust(icol,k) = bedustgt1(icol,k)*deltah
!SO4
           dod5503dlt1_so4(icol,k) = (bes4lt1t(icol,k)+bs4lt1xt(icol,k) &   ! condensate + n-mode (11)
                    + Nnatk(icol,k,1)*bebglt1(icol,k,1)                 &   ! background, SO4(Ait) mode (1)
                    + Nnatk(icol,k,5)*bebglt1(icol,k,5))*deltah             ! background, SO4(Ait75) mode (5)
           dod5503dgt1_so4(icol,k) = (bes4gt1t(icol,k)+bs4gt1xt(icol,k) &   ! condensate + n-mode (11)
                    + Nnatk(icol,k,1)*bebggt1(icol,k,1)                 &   ! background, SO4(Ait) mode (1)
                    + Nnatk(icol,k,5)*bebggt1(icol,k,5))*deltah             ! background, SO4(Ait75) mode (5)
!BC
           dod5503dlt1_bc(icol,k) =  (bebclt1t(icol,k)+bbclt1xt(icol,k) &   ! coagulated + n-mode BC (12)
                    + Nnatk(icol,k,2)*bebglt1(icol,k,2)                 &   ! background, BC(Ait) mode (2)
             + vaitbc*Nnatk(icol,k,4)*bebglt1(icol,k,4)                 &   ! background in OC&BC(Ait) mode (4)
                    + Nnatk(icol,k,0)*bebglt1(icol,k,0))*deltah             ! background, BC(ax) mode (0)
           dod5503dgt1_bc(icol,k) =  (bebcgt1t(icol,k)+bbcgt1xt(icol,k) &   ! coagulated + n-mode BC (12)
                    + Nnatk(icol,k,2)*bebggt1(icol,k,2)                 &   ! background, BC(Ait) mode (2)
             + vaitbc*Nnatk(icol,k,4)*bebggt1(icol,k,4)                 &   ! background in OC&BC(Ait) mode (4)
                    + Nnatk(icol,k,0)*bebggt1(icol,k,0))*deltah             ! background, BC(ax) mode (0)
!OC
           dod5503dlt1_pom(icol,k) = (beoclt1t(icol,k)+boclt1xt(icol,k) &   ! coagulated + n-mode OC (13)
                    + Nnatk(icol,k,3)*bebglt1(icol,k,3)                 &   ! background, OC(Ait) mode (3)
       + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*bebglt1(icol,k,4))*deltah             ! background in OC&BC(Ait) mode (4)
           dod5503dgt1_pom(icol,k) = (beocgt1t(icol,k)+bocgt1xt(icol,k) &   ! coagulated + n-mode OC (13)
                    + Nnatk(icol,k,3)*bebggt1(icol,k,3)                 &   ! background, OC(Ait) mode (3)
       + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*bebggt1(icol,k,4))*deltah             ! background in OC&BC(Ait) mode (4)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!          Column integrated optical depths/abs., total and for each constituent
           dod440(icol)         = dod440(icol)+dod4403d(icol,k)
           abs440(icol)         = abs440(icol)+abs4403d(icol,k)
           dod500(icol)         = dod500(icol)+dod5003d(icol,k)
           abs500(icol)         = abs500(icol)+abs5003d(icol,k)
           dod550(icol)         = dod550(icol)+dod5503d(icol,k)
           abs550(icol)         = abs550(icol)+abs5503d(icol,k)
           abs550alt(icol)      = abs550alt(icol)+abs5503dalt(icol,k)
           dod670(icol)         = dod670(icol)+dod6703d(icol,k)
           abs670(icol)         = abs670(icol)+abs6703d(icol,k)
           dod870(icol)         = dod870(icol)+dod8703d(icol,k)
           abs870(icol)         = abs870(icol)+abs8703d(icol,k)
! Added abs components
           abs550_ss(icol)      = abs550_ss(icol)+abs5503d_ss(icol,k)
           abs550_dust(icol)    = abs550_dust(icol)+abs5503d_dust(icol,k)
           abs550_so4(icol)     = abs550_so4(icol)+abs5503d_so4(icol,k)
           abs550_bc(icol)      = abs550_bc(icol)+abs5503d_bc(icol,k)
           abs550_pom(icol)     = abs550_pom(icol)+abs5503d_pom(icol,k)
!
           dod440_ss(icol)      = dod440_ss(icol)+dod4403d_ss(icol,k)
           dod440_dust(icol)    = dod440_dust(icol)+dod4403d_dust(icol,k)
           dod440_so4(icol)     = dod440_so4(icol)+dod4403d_so4(icol,k)
           dod440_bc(icol)      = dod440_bc(icol)+dod4403d_bc(icol,k)
           dod440_pom(icol)     = dod440_pom(icol)+dod4403d_pom(icol,k)
           dod500_ss(icol)      = dod500_ss(icol)+dod5003d_ss(icol,k)
           dod500_dust(icol)    = dod500_dust(icol)+dod5003d_dust(icol,k)
           dod500_so4(icol)     = dod500_so4(icol)+dod5003d_so4(icol,k)
           dod500_bc(icol)      = dod500_bc(icol)+dod5003d_bc(icol,k)
           dod500_pom(icol)     = dod500_pom(icol)+dod5003d_pom(icol,k)
           dod550_ss(icol)      = dod550_ss(icol)+dod5503d_ss(icol,k)
           dod550_dust(icol)    = dod550_dust(icol)+dod5503d_dust(icol,k)
           dod550_so4(icol)     = dod550_so4(icol)+dod5503d_so4(icol,k)
           dod550_bc(icol)      = dod550_bc(icol)+dod5503d_bc(icol,k)
           dod550_pom(icol)     = dod550_pom(icol)+dod5503d_pom(icol,k)
           dod670_ss(icol)      = dod670_ss(icol)+dod6703d_ss(icol,k)
           dod670_dust(icol)    = dod670_dust(icol)+dod6703d_dust(icol,k)
           dod670_so4(icol)     = dod670_so4(icol)+dod6703d_so4(icol,k)
           dod670_bc(icol)      = dod670_bc(icol)+dod6703d_bc(icol,k)
           dod670_pom(icol)     = dod670_pom(icol)+dod6703d_pom(icol,k)
           dod870_ss(icol)      = dod870_ss(icol)+dod8703d_ss(icol,k)
           dod870_dust(icol)    = dod870_dust(icol)+dod8703d_dust(icol,k)
           dod870_so4(icol)     = dod870_so4(icol)+dod8703d_so4(icol,k)
           dod870_bc(icol)      = dod870_bc(icol)+dod8703d_bc(icol,k)
           dod870_pom(icol)     = dod870_pom(icol)+dod8703d_pom(icol,k)
!old           dod865_ss(icol)      = dod865_ss(icol)+dod8653d_ss(icol,k)
!old           dod865_dust(icol)    = dod865_dust(icol)+dod8653d_dust(icol,k)
!old           dod865_so4(icol)     = dod865_so4(icol)+dod8653d_so4(icol,k)
!old           dod865_bc(icol)      = dod865_bc(icol)+dod8653d_bc(icol,k)
!old           dod865_pom(icol)     = dod865_pom(icol)+dod8653d_pom(icol,k)
           dod550lt1_ss(icol)   = dod550lt1_ss(icol)+dod5503dlt1_ss(icol,k)
           dod550gt1_ss(icol)   = dod550gt1_ss(icol)+dod5503dgt1_ss(icol,k)
           dod550lt1_dust(icol) = dod550lt1_dust(icol)+dod5503dlt1_dust(icol,k)
           dod550gt1_dust(icol) = dod550gt1_dust(icol)+dod5503dgt1_dust(icol,k)
           dod550lt1_so4(icol)  = dod550lt1_so4(icol)+dod5503dlt1_so4(icol,k) 
           dod550gt1_so4(icol)  = dod550gt1_so4(icol)+dod5503dgt1_so4(icol,k) 
           dod550lt1_bc(icol)   = dod550lt1_bc(icol)+dod5503dlt1_bc(icol,k)  
           dod550gt1_bc(icol)   = dod550gt1_bc(icol)+dod5503dgt1_bc(icol,k)  
           dod550lt1_pom(icol)  = dod550lt1_pom(icol)+dod5503dlt1_pom(icol,k) 
           dod550gt1_pom(icol)  = dod550gt1_pom(icol)+dod5503dgt1_pom(icol,k)
! Cloudfree components 
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
          enddo ! k

! Cloudfree conditions
         clearod550(icol)=cloudfree(icol)*dod550(icol)
!old         clearod865(icol)=cloudfree(icol)*(dod865_so4(icol)+dod865_bc(icol)+ &
!old            dod865_pom(icol)+dod865_ss(icol)+dod865_dust(icol))
         clearod870(icol)=cloudfree(icol)*dod870(icol)
         clearabs550(icol)=cloudfree(icol)*abs550(icol)
         clearabs550alt(icol)=cloudfree(icol)*abs550alt(icol)
         enddo  ! icol

!       Dry parameters of each aerosol component 
!       BC(ax) mode
        call intdrypar0(lchnk, ncol, Nnatk,                          & 
           cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125, & 
           cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125, &
           cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol,&
           cknorm,cknlt05,ckngt125)
!       SO4(Ait,n), BC(Ait,n) and OC(Ait,n) modes
        call intdrypar1to3(lchnk, ncol, Nnatk, Cam,                  & 
           cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125, & 
           cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125, &
           cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol,&
           aaerosn,aaeroln,vaerosn,vaeroln,cknorm,cknlt05,ckngt125)
!       BC&OC(Ait) mode   ------ fcm not valid here (=0). Use faitbc instead
        call intdrypar4(lchnk, ncol, Nnatk, Cam, faitbc, faqm4,      &
           cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125, &
           cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125, &
           cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol)
!       BC&OC(Ait and n) modes   ------ fcm not valid here (=0). Use fait and fnbc instead
        call intdrypar4r(lchnk, ncol, Nnatk, faitbc, fnbc,           &
           aaerosn,aaeroln,vaerosn,vaeroln,cknorm,cknlt05,ckngt125)
!       SO4(Ait75) (5), mineral (6-7) and Sea-salt (8-10) modes:
        call intdrypar5to10(lchnk, ncol, Nnatk, Cam, fcm, fbcm, faqm,& 
           cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125, & 
           cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125, &
           cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol,&
           cknorm,cknlt05,ckngt125)

        do k=1,pver
           do icol=1,ncol
!           mineral and sea-salt background concentrations, internally mixed
            c_mi(icol,k)    = Nnatk(icol,k,6)*cintbg(icol,k,6)   &
                             +Nnatk(icol,k,7)*cintbg(icol,k,7)
            c_mi05(icol,k)  = Nnatk(icol,k,6)*cintbg05(icol,k,6) &
                             +Nnatk(icol,k,7)*cintbg05(icol,k,7)
            c_mi125(icol,k) = Nnatk(icol,k,6)*cintbg125(icol,k,6)& 
                             +Nnatk(icol,k,7)*cintbg125(icol,k,7) 
            c_ss(icol,k)    = Nnatk(icol,k,8)*cintbg(icol,k,8)   &
                             +Nnatk(icol,k,9)*cintbg(icol,k,9)    &
                             +Nnatk(icol,k,10)*cintbg(icol,k,10)
            c_ss05(icol,k)  = Nnatk(icol,k,8)*cintbg05(icol,k,8) &
                             +Nnatk(icol,k,9)*cintbg05(icol,k,9)  &
                             +Nnatk(icol,k,10)*cintbg05(icol,k,10)
            c_ss125(icol,k) = Nnatk(icol,k,8)*cintbg125(icol,k,8)&
                             +Nnatk(icol,k,9)*cintbg125(icol,k,9) &
                             +Nnatk(icol,k,10)*cintbg125(icol,k,10)
!           internally mixed bc and oc (from coagulation) and so4 concentrations 
!           (sa=so4(aq) and sc=so4(cond+coag), separated because of different density: 
!           necessary for calculation of volume fractions!), and total aerosol surface 
!           areas and volumes. 
            c_bc(icol,k)=0.0_r8
            c_bc05(icol,k)=0.0_r8
            c_bc125(icol,k)=0.0_r8
            c_oc(icol,k)=0.0_r8
            c_oc05(icol,k)=0.0_r8
            c_oc125(icol,k)=0.0_r8
            c_sa(icol,k)=0.0_r8
            c_sa05(icol,k)=0.0_r8
            c_sa125(icol,k)=0.0_r8
            c_sc(icol,k)=0.0_r8
            c_sc05(icol,k)=0.0_r8
            c_sc125(icol,k)=0.0_r8
            aaeros_tot(icol,k)=0.0_r8
            aaerol_tot(icol,k)=0.0_r8
            vaeros_tot(icol,k)=0.0_r8
            vaerol_tot(icol,k)=0.0_r8
            do i=0,nbmodes
               if(i>=5) then
             c_bc(icol,k)    = c_bc(icol,k) &
                            +Nnatk(icol,k,i)*cintbc(icol,k,i)
             c_bc05(icol,k)  = c_bc05(icol,k) &
                            +Nnatk(icol,k,i)*cintbc05(icol,k,i)
             c_bc125(icol,k) = c_bc125(icol,k) &
                            +Nnatk(icol,k,i)*cintbc125(icol,k,i)
             c_oc(icol,k)    = c_oc(icol,k) &
                            +Nnatk(icol,k,i)*cintoc(icol,k,i)
             c_oc05(icol,k)  = c_oc05(icol,k) &
                            +Nnatk(icol,k,i)*cintoc05(icol,k,i)
             c_oc125(icol,k) = c_oc125(icol,k) &
                            +Nnatk(icol,k,i)*cintoc125(icol,k,i)
               endif
             c_sa(icol,k)    = c_sa(icol,k) &
                            +Nnatk(icol,k,i)*cintsa(icol,k,i)
             c_sa05(icol,k)  = c_sa05(icol,k) &
                            +Nnatk(icol,k,i)*cintsa05(icol,k,i)
             c_sa125(icol,k) = c_sa125(icol,k) &
                            +Nnatk(icol,k,i)*cintsa125(icol,k,i)
             c_sc(icol,k)    = c_sc(icol,k) &
                            +Nnatk(icol,k,i)*cintsc(icol,k,i)
             c_sc05(icol,k)  = c_sc05(icol,k) &
                            +Nnatk(icol,k,i)*cintsc05(icol,k,i)
             c_sc125(icol,k) = c_sc125(icol,k) &
                            +Nnatk(icol,k,i)*cintsc125(icol,k,i)
             aaeros_tot(icol,k) = aaeros_tot(icol,k) &
                            +Nnatk(icol,k,i)*aaeros(icol,k,i)
             aaerol_tot(icol,k) = aaerol_tot(icol,k) &
                            +Nnatk(icol,k,i)*aaerol(icol,k,i)
             vaeros_tot(icol,k) =vaeros_tot(icol,k) &
                            +Nnatk(icol,k,i)*vaeros(icol,k,i)
             vaerol_tot(icol,k) = vaerol_tot(icol,k) &
                            +Nnatk(icol,k,i)*vaerol(icol,k,i)
            enddo
!           add dry aerosol area and volume of externally mixed modes
            do i=nbmp1,nmodes
             aaeros_tot(icol,k) = aaeros_tot(icol,k) &
                            +Nnatk(icol,k,i)*aaerosn(icol,k,i)
             aaerol_tot(icol,k) = aaerol_tot(icol,k) &
                            +Nnatk(icol,k,i)*aaeroln(icol,k,i)
             vaeros_tot(icol,k) =vaeros_tot(icol,k) &
                            +Nnatk(icol,k,i)*vaerosn(icol,k,i)
             vaerol_tot(icol,k) = vaerol_tot(icol,k) &
                            +Nnatk(icol,k,i)*vaeroln(icol,k,i)
            end do
!c_er3d           
!           Effective radii for particles smaller and greater than 0.5um, 
!           and for all radii, in each layer
            erlt053d(icol,k)=3.0_r8*vaeros_tot(icol,k) &
                             /(aaeros_tot(icol,k)+eps)
            ergt053d(icol,k)=3.0_r8*vaerol_tot(icol,k) &
                             /(aaerol_tot(icol,k)+eps)
            er3d(icol,k)=3.0_r8*(vaeros_tot(icol,k)+vaerol_tot(icol,k)) &
                          /(aaeros_tot(icol,k)+aaerol_tot(icol,k)+eps)
!c_er3d
!           column integrated dry aerosol surface areas and volumes
!           for r<0.5um and r>0.5um (s and l, respectively).
            aaercols(icol)=aaercols(icol)+aaeros_tot(icol,k)
            aaercoll(icol)=aaercoll(icol)+aaerol_tot(icol,k)
            vaercols(icol)=vaercols(icol)+vaeros_tot(icol,k)
            vaercoll(icol)=vaercoll(icol)+vaerol_tot(icol,k)
!2oct07
            rhobcbyoc=rhopart(l_bc_ni)/rhopart(l_om_ni)
!2oct07
!           then add background and externally mixed BC, OC and SO4 to mass concentrations
              c_bc(icol,k)   = c_bc(icol,k)                                     &
                            +Nnatk(icol,k,2)*cintbg(icol,k,2)                   &
!2oct07                            +Nnatk(icol,k,4)*cintbg(icol,k,4)*faitbc(icol,k)    &
                            +Nnatk(icol,k,4)*cintbg(icol,k,4)*faitbc(icol,k)*rhobcbyoc &
                            +Nnatk(icol,k,0)*cintbg(icol,k,0)                   &
                            +Nnatk(icol,k,12)*cknorm(icol,k,12)                 &
!2oct07                            +Nnatk(icol,k,14)*cknorm(icol,k,14)*fnbc(icol,k)
                            +Nnatk(icol,k,14)*cknorm(icol,k,14)*fnbc(icol,k)*rhobcbyoc
             c_bc05(icol,k)  = c_bc05(icol,k)                                   &
                            +Nnatk(icol,k,2)*cintbg05(icol,k,2)                 &
!2oct07                            +Nnatk(icol,k,4)*cintbg05(icol,k,4)*faitbc(icol,k)  &
                            +Nnatk(icol,k,4)*cintbg05(icol,k,4)*faitbc(icol,k)*rhobcbyoc &
                            +Nnatk(icol,k,0)*cintbg05(icol,k,0)                 &
                            +Nnatk(icol,k,12)*cknlt05(icol,k,12)                &
!2oct07                            +Nnatk(icol,k,14)*cknlt05(icol,k,14)*fnbc(icol,k)
                            +Nnatk(icol,k,14)*cknlt05(icol,k,14)*fnbc(icol,k)*rhobcbyoc
             c_bc125(icol,k) = c_bc125(icol,k)                                  &
                            +Nnatk(icol,k,2)*cintbg125(icol,k,2)                &
!2oct07                            +Nnatk(icol,k,4)*cintbg125(icol,k,4)*faitbc(icol,k) &
                            +Nnatk(icol,k,4)*cintbg125(icol,k,4)*faitbc(icol,k)*rhobcbyoc &
                            +Nnatk(icol,k,0)*cintbg125(icol,k,0)                &
                            +Nnatk(icol,k,12)*ckngt125(icol,k,12)               &
!2oct07                            +Nnatk(icol,k,14)*ckngt125(icol,k,14)*fnbc(icol,k)
                            +Nnatk(icol,k,14)*ckngt125(icol,k,14)*fnbc(icol,k)*rhobcbyoc
             c_oc(icol,k)    = c_oc(icol,k)                                           &
                            +Nnatk(icol,k,3)*cintbg(icol,k,3)                         &
                            +Nnatk(icol,k,4)*cintbg(icol,k,4)*(1.0_r8-faitbc(icol,k))    &
                            +Nnatk(icol,k,13)*cknorm(icol,k,13)                       &
                            +Nnatk(icol,k,14)*cknorm(icol,k,14)*(1.0_r8-fnbc(icol,k))
             c_oc05(icol,k)  = c_oc05(icol,k)                                         &
                            +Nnatk(icol,k,3)*cintbg05(icol,k,3)                       &
                            +Nnatk(icol,k,4)*cintbg05(icol,k,4)*(1.0_r8-faitbc(icol,k))  &
                            +Nnatk(icol,k,13)*cknlt05(icol,k,13)                      &
                            +Nnatk(icol,k,14)*cknlt05(icol,k,14)*(1.0_r8-fnbc(icol,k))
             c_oc125(icol,k) = c_oc125(icol,k)                                        &
                            +Nnatk(icol,k,3)*cintbg125(icol,k,3)                      &
                            +Nnatk(icol,k,4)*cintbg125(icol,k,4)*(1.0_r8-faitbc(icol,k)) &
                            +Nnatk(icol,k,13)*ckngt125(icol,k,13)                     &
                            +Nnatk(icol,k,14)*ckngt125(icol,k,14)*(1.0_r8-fnbc(icol,k))
             c_s4(icol,k)    = c_sa(icol,k)+c_sc(icol,k)          &
                            +Nnatk(icol,k,1)*cintbg(icol,k,1)     &
                            +Nnatk(icol,k,5)*cintbg(icol,k,5)     &
                            +Nnatk(icol,k,11)*cknorm(icol,k,11)
             c_s405(icol,k)  = c_sa05(icol,k)+c_sc05(icol,k)      &
                            +Nnatk(icol,k,1)*cintbg05(icol,k,1)   &
                            +Nnatk(icol,k,5)*cintbg05(icol,k,5)   &
                            +Nnatk(icol,k,11)*cknlt05(icol,k,11)
             c_s4125(icol,k) = c_sa125(icol,k)+c_sc125(icol,k)    &
                            +Nnatk(icol,k,1)*cintbg125(icol,k,1)  &
                            +Nnatk(icol,k,5)*cintbg125(icol,k,5)  &
                            +Nnatk(icol,k,11)*ckngt125(icol,k,11)
!            convert from S to SO4 concentrations      
             c_s4(icol,k)=c_s4(icol,k)/3._r8
             c_s405(icol,k)=c_s405(icol,k)/3._r8
             c_s4125(icol,k)=c_s4125(icol,k)/3._r8
           end do ! icol
        enddo     ! k

!       Effective, column integrated, radii for particles
!       smaller and greater than 0.5um, and for all radii
        do icol=1,ncol
            derlt05(icol)=3.0_r8*vaercols(icol)/(aaercols(icol)+eps)
            dergt05(icol)=3.0_r8*vaercoll(icol)/(aaercoll(icol)+eps)
            der(icol)=3.0_r8*(vaercols(icol)+vaercoll(icol)) &
                       /(aaercols(icol)+aaercoll(icol)+eps)
        enddo

        do icol=1,ncol
         do k=1,pver
!         Layer thickness, unit km
          deltah=1.e-4_r8*(pint(icol,k+1)-pint(icol,k))/(rhoda(icol,k)*9.8_r8)
!         Modal and total mass concentrations for clean and dry aerosol, 
!         i.e. not including coag./cond./Aq. BC,OC,SO4 or condensed water. 
!         Units: ug/m3 for concentrations and mg/m2 (--> kg/m2 later) for mass loading.  
          do i=0,nmodes
            ck(icol,k,i)=cknorm(icol,k,i)*Nnatk(icol,k,i)
            dload3d(icol,k,i)=ck(icol,k,i)*deltah
            dload(icol,i)=dload(icol,i)+dload3d(icol,k,i)
!moved up            n_aer(icol,k)=n_aer(icol,k)+Nnatk(icol,k,i)
          enddo
          nnat_0(icol,k) =Nnatk(icol,k,0)
          nnat_1(icol,k) =Nnatk(icol,k,1)
          nnat_2(icol,k) =Nnatk(icol,k,2)
!=0          nnat_3(icol,k) =Nnatk(icol,k,3)
          nnat_4(icol,k) =Nnatk(icol,k,4)
          nnat_5(icol,k) =Nnatk(icol,k,5)
          nnat_6(icol,k) =Nnatk(icol,k,6)
          nnat_7(icol,k) =Nnatk(icol,k,7)
          nnat_8(icol,k) =Nnatk(icol,k,8)
          nnat_9(icol,k) =Nnatk(icol,k,9)
          nnat_10(icol,k)=Nnatk(icol,k,10)
          nnat_11(icol,k)=Nnatk(icol,k,11)
          nnat_12(icol,k)=Nnatk(icol,k,12)
          nnat_13(icol,k)=Nnatk(icol,k,13)
          nnat_14(icol,k)=Nnatk(icol,k,14)
!         mineral and sea-salt mass concentrations
          cmin(icol,k)=ck(icol,k,6)+ck(icol,k,7)               
          cseas(icol,k)=ck(icol,k,8)+ck(icol,k,9)+ck(icol,k,10)
!          Aerocom: Condensed water loading (mg_m2)
          daerh2o(icol)=daerh2o(icol)+Cwater(icol,k)*deltah
         end do  ! k
         dload_mi(icol)=dload(icol,6)+dload(icol,7)
         dload_ss(icol)=dload(icol,8)+dload(icol,9)+dload(icol,10)
        end do   ! icol

!       Internally and externally mixed dry concentrations (ug/m3) of  
!       SO4, BC and OC, for all r, r<0.5um and r>1.25um...
!-        call outfld('C_BCPM  ',c_bc   ,pcols,lchnk)
!-        call outfld('C_BC05  ',c_bc05 ,pcols,lchnk)
!-        call outfld('C_BC125 ',c_bc125,pcols,lchnk)
!-        call outfld('C_OCPM  ',c_oc   ,pcols,lchnk)
!-        call outfld('C_OC05  ',c_oc05 ,pcols,lchnk)
!-        call outfld('C_OC125 ',c_oc125,pcols,lchnk)
!-        call outfld('C_S4PM  ',c_s4   ,pcols,lchnk)
!-        call outfld('C_S405  ',c_s405 ,pcols,lchnk)
!-        call outfld('C_S4125 ',c_s4125,pcols,lchnk)
!       ... and of background components for all r, r<0.5um and r>1.25um
!-        call outfld('C_MIPM  ',c_mi   ,pcols,lchnk)
!-        call outfld('C_MI05  ',c_mi05 ,pcols,lchnk)
!-        call outfld('C_MI125 ',c_mi125,pcols,lchnk)
!-        call outfld('C_SSPM  ',c_ss   ,pcols,lchnk)
!-        call outfld('C_SS05  ',c_ss05 ,pcols,lchnk)
!-        call outfld('C_SS125 ',c_ss125,pcols,lchnk)
!       total (all r) dry concentrations (ug/m3) and loadings (mg/m2) of 
!       clean background aerosols, and loadings of externally mixed SO4,
!       BC and OC
        call outfld('DLOAD_MI',dload_mi,pcols,lchnk)
        call outfld('DLOAD_SS',dload_ss,pcols,lchnk)
!       condensed water mmr (kg/kg)
        call outfld('MMR_AH2O',mmr_aerh2o,pcols,lchnk)
!       condensed water loading (mg/m2)
        call outfld('DAERH2O ',daerh2o ,pcols,lchnk)
!       number concentrations (1/cm3)
        call outfld('NNAT_0  ',nnat_0 ,pcols,lchnk)
        call outfld('NNAT_1  ',nnat_1 ,pcols,lchnk)
        call outfld('NNAT_2  ',nnat_2 ,pcols,lchnk)
!=0        call outfld('NNAT_3  ',nnat_3 ,pcols,lchnk)
        call outfld('NNAT_4  ',nnat_4 ,pcols,lchnk)
        call outfld('NNAT_5  ',nnat_5 ,pcols,lchnk)
        call outfld('NNAT_6  ',nnat_6 ,pcols,lchnk)
        call outfld('NNAT_7  ',nnat_7 ,pcols,lchnk)
        call outfld('NNAT_8  ',nnat_8 ,pcols,lchnk)
        call outfld('NNAT_9  ',nnat_9 ,pcols,lchnk)
        call outfld('NNAT_10 ',nnat_10,pcols,lchnk)
        call outfld('NNAT_11 ',nnat_11,pcols,lchnk)
        call outfld('NNAT_12 ',nnat_12,pcols,lchnk)
        call outfld('NNAT_13 ',nnat_13,pcols,lchnk)
        call outfld('NNAT_14 ',nnat_14,pcols,lchnk)
        call outfld('AIRMASS ',airmass,pcols,lchnk)
!c_er3d 
!       effective dry radii (um) in each layer
!        call outfld('ERLT053D',erlt053d,pcols,lchnk)
!        call outfld('ERGT053D',ergt053d,pcols,lchnk)
!        call outfld('ER3D    ',er3d    ,pcols,lchnk)
!c_er3d           
!       column integrated effective dry radii (um)
        call outfld('DERLT05 ',derlt05,pcols,lchnk)
        call outfld('DERGT05 ',dergt05,pcols,lchnk)
        call outfld('DER     ',der    ,pcols,lchnk)
!
!       extinction, absorption (m-1) and backscatter coefficients (m-1 sr-1)
        call outfld('EC550AER',ec550_aer,pcols,lchnk)
        call outfld('ABS550_A',abs550_aer,pcols,lchnk)
        call outfld('BS550AER',bs550_aer,pcols,lchnk)
!
!       optical depths and absorption as requested by AeroCom
!old       notation: 1=550, 2=865, 3=3D, D=DOD, A=ABS, LT=LT1, GT=GT1
!       notation: 3=3D, D=DOD, A=ABS, LT=d<1um, GT=d>1um
        call outfld('DOD440  ',dod440  ,pcols,lchnk)
!        call outfld('ABS440  ',abs440  ,pcols,lchnk)
        call outfld('DOD500  ',dod500  ,pcols,lchnk)
!        call outfld('ABS500  ',abs500  ,pcols,lchnk)
        call outfld('DOD550  ',dod550  ,pcols,lchnk)
        call outfld('ABS550  ',abs550  ,pcols,lchnk)
        call outfld('ABS550AL',abs550alt,pcols,lchnk)
        call outfld('DOD670  ',dod670  ,pcols,lchnk)
!        call outfld('ABS670  ',abs670  ,pcols,lchnk)
        call outfld('DOD870  ',dod870  ,pcols,lchnk)
!        call outfld('ABS870  ',abs870  ,pcols,lchnk)
! Added clear sky variables
        call outfld('CDOD550 ',clearod550  ,pcols,lchnk)
        call outfld('CABS550 ',clearabs550  ,pcols,lchnk)
        call outfld('CABS550A',clearabs550alt,pcols,lchnk)
!old        call outfld('CDOD865 ',clearod865  ,pcols,lchnk)
        call outfld('CDOD870 ',clearod870  ,pcols,lchnk)
! OS Added abs 2d
!old        call outfld('A1_SS   ',abs550_ss  ,pcols,lchnk)
!old        call outfld('A1_DUST ',abs550_dust,pcols,lchnk)
!old        call outfld('A1_SO4  ',abs550_so4 ,pcols,lchnk)
!old        call outfld('A1_BC   ',abs550_bc  ,pcols,lchnk)
!old        call outfld('A1_POM  ',abs550_pom ,pcols,lchnk)
        call outfld('A550_SS ',abs550_ss  ,pcols,lchnk)
        call outfld('A550_DU ',abs550_dust,pcols,lchnk)
        call outfld('A550_SO4',abs550_so4 ,pcols,lchnk)
        call outfld('A550_BC ',abs550_bc  ,pcols,lchnk)
        call outfld('A550_POM',abs550_pom ,pcols,lchnk)
!
!old        call outfld('D1_SS   ',dod550_ss  ,pcols,lchnk)
!old        call outfld('D1_DUST ',dod550_dust,pcols,lchnk)
!old        call outfld('D1_SO4  ',dod550_so4 ,pcols,lchnk)
!old        call outfld('D1_BC   ',dod550_bc  ,pcols,lchnk)
!old        call outfld('D1_POM  ',dod550_pom ,pcols,lchnk)
!old        call outfld('D2_SS   ',dod865_ss  ,pcols,lchnk)
!old        call outfld('D2_DUST ',dod865_dust,pcols,lchnk)
!old        call outfld('D2_SO4  ',dod865_so4 ,pcols,lchnk)
!old        call outfld('D2_BC   ',dod865_bc  ,pcols,lchnk)
!old        call outfld('D2_POM  ',dod865_pom ,pcols,lchnk)
        call outfld('D440_SS ',dod440_ss  ,pcols,lchnk)
        call outfld('D440_DU ',dod440_dust,pcols,lchnk)
        call outfld('D440_SO4',dod440_so4 ,pcols,lchnk)
        call outfld('D440_BC ',dod440_bc  ,pcols,lchnk)
        call outfld('D440_POM',dod440_pom ,pcols,lchnk)
        call outfld('D500_SS ',dod500_ss  ,pcols,lchnk)
        call outfld('D500_DU ',dod500_dust,pcols,lchnk)
        call outfld('D500_SO4',dod500_so4 ,pcols,lchnk)
        call outfld('D500_BC ',dod500_bc  ,pcols,lchnk)
        call outfld('D500_POM',dod500_pom ,pcols,lchnk)
        call outfld('D550_SS ',dod550_ss  ,pcols,lchnk)
        call outfld('D550_DU ',dod550_dust,pcols,lchnk)
        call outfld('D550_SO4',dod550_so4 ,pcols,lchnk)
        call outfld('D550_BC ',dod550_bc  ,pcols,lchnk)
        call outfld('D550_POM',dod550_pom ,pcols,lchnk)
        call outfld('D670_SS ',dod670_ss  ,pcols,lchnk)
        call outfld('D670_DU ',dod670_dust,pcols,lchnk)
        call outfld('D670_SO4',dod670_so4 ,pcols,lchnk)
        call outfld('D670_BC ',dod670_bc  ,pcols,lchnk)
        call outfld('D670_POM',dod670_pom ,pcols,lchnk)
        call outfld('D870_SS ',dod870_ss  ,pcols,lchnk)
        call outfld('D870_DU ',dod870_dust,pcols,lchnk)
        call outfld('D870_SO4',dod870_so4 ,pcols,lchnk)
        call outfld('D870_BC ',dod870_bc  ,pcols,lchnk)
        call outfld('D870_POM',dod870_pom ,pcols,lchnk)
        call outfld('DLT_SS  ',dod550lt1_ss,pcols,lchnk)
        call outfld('DGT_SS  ',dod550gt1_ss,pcols,lchnk)
        call outfld('DLT_DUST',dod550lt1_dust,pcols,lchnk)
        call outfld('DGT_DUST',dod550gt1_dust,pcols,lchnk)
        call outfld('DLT_SO4 ',dod550lt1_so4,pcols,lchnk)
        call outfld('DGT_SO4 ',dod550gt1_so4,pcols,lchnk)
        call outfld('DLT_BC  ',dod550lt1_bc,pcols,lchnk)
        call outfld('DGT_BC  ',dod550gt1_bc,pcols,lchnk)
        call outfld('DLT_POM ',dod550lt1_pom,pcols,lchnk)
        call outfld('DGT_POM ',dod550gt1_pom,pcols,lchnk)
!        call outfld('DOD5503D',dod5503d,pcols,lchnk)
        call outfld('ABS5503D',abs5503d,pcols,lchnk)
!old        call outfld('D13_SS  ',dod5503d_ss ,pcols,lchnk)
!old        call outfld('A13_SS  ',abs5503d_ss ,pcols,lchnk)
!old        call outfld('D13_DUST',dod5503d_dust,pcols,lchnk)
!old        call outfld('A13_DUST',abs5503d_dust,pcols,lchnk)
!old        call outfld('D13_SO4 ',dod5503d_so4,pcols,lchnk) 
!old        call outfld('A13_SO4 ',abs5503d_so4,pcols,lchnk)
!old        call outfld('D13_BC  ',dod5503d_bc ,pcols,lchnk)
!old        call outfld('A13_BC  ',abs5503d_bc ,pcols,lchnk)
!old        call outfld('D13_POM ',dod5503d_pom,pcols,lchnk)
!old        call outfld('A13_POM ',abs5503d_pom,pcols,lchnk)
!-        call outfld('D443_SS ',dod4403d_ss  ,pcols,lchnk)
!-        call outfld('D443_DU ',dod4403d_dust,pcols,lchnk)
!-        call outfld('D443_SO4',dod4403d_so4 ,pcols,lchnk)
!-        call outfld('D443_BC ',dod4403d_bc  ,pcols,lchnk)
!-        call outfld('D443_POM',dod4403d_pom ,pcols,lchnk)
!-        call outfld('D503_SS ',dod5003d_ss  ,pcols,lchnk)
!-        call outfld('D503_DU ',dod5003d_dust,pcols,lchnk)
!-        call outfld('D503_SO4',dod5003d_so4 ,pcols,lchnk)
!-        call outfld('D503_BC ',dod5003d_bc  ,pcols,lchnk)
!-        call outfld('D503_POM',dod5003d_pom ,pcols,lchnk)
        call outfld('D553_SS ',dod5503d_ss  ,pcols,lchnk)
        call outfld('D553_DU ',dod5503d_dust,pcols,lchnk)
        call outfld('D553_SO4',dod5503d_so4 ,pcols,lchnk)
        call outfld('D553_BC ',dod5503d_bc  ,pcols,lchnk)
        call outfld('D553_POM',dod5503d_pom ,pcols,lchnk)
!-        call outfld('D673_SS ',dod6703d_ss  ,pcols,lchnk)
!-        call outfld('D673_DU ',dod6703d_dust,pcols,lchnk)
!-        call outfld('D673_SO4',dod6703d_so4 ,pcols,lchnk)
!-        call outfld('D673_BC ',dod6703d_bc  ,pcols,lchnk)
!-        call outfld('D673_POM',dod6703d_pom ,pcols,lchnk)
!-        call outfld('D873_SS ',dod8703d_ss  ,pcols,lchnk)
!-        call outfld('D873_DU ',dod8703d_dust,pcols,lchnk)
!-        call outfld('D873_SO4',dod8703d_so4 ,pcols,lchnk)
!-        call outfld('D873_BC ',dod8703d_bc  ,pcols,lchnk)
!-        call outfld('D873_POM',dod8703d_pom ,pcols,lchnk)

 
!     Extra AeroCom diagnostics requiring tablel ook-ups with RH=0%. Note:
!     intaeropt0 gives the same results as with ambient RH and is not called. 

      do k=1,pver
        do icol=1,ncol
          rhum(icol,k) = 0.01_r8
          ec550dry_aer(icol,k)=0.0_r8 
          abs550dry_aer(icol,k)=0.0_r8 
        end do
      end do
      do icol=1,ncol
        dod550dry(icol)=0.0_r8 
        abs550dry(icol)=0.0_r8 
      end do

!     SO4(Ait), BC(Ait) and OC(Ait) modes:
      mplus10=0
        call intaeropt1to3(lchnk, ncol, rhum, mplus10, Nnatk, Cam, &
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebglt1, bebggt1, bebclt1, bebcgt1,                         &
           beoclt1, beocgt1, bes4lt1, bes4gt1,                         &
           backsc550, babg550, babc550, baoc550, basu550)

!     BC&OC(Ait) mode:    ------ fcm invalid here (=0). Using faitbc instead
      mplus10=0
        call intaeropt4(lchnk, ncol, rhum, mplus10, Nnatk, Cam, faitbc, faqm4, &
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebglt1, bebggt1, bebclt1, bebcgt1,                         &
           beoclt1, beocgt1, bes4lt1, bes4gt1,                         &
           backsc550, babg550, babc550, baoc550, basu550)
  
!     SO4(Ait75) (5), Mineral (6-7) and Sea-salt (8-10) modes:
        call intaeropt5to10(lchnk, ncol, rhum, Nnatk, Cam, fcm, fbcm, faqm, &
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebglt1, bebggt1, bebclt1, bebcgt1,                         &
           beoclt1, beocgt1, bes4lt1, bes4gt1,                         &
           backsc550, babg550, babc550, baoc550, basu550)

!     then to the externally mixed SO4(n), BC(n) and OC(n) modes:
      mplus10=1
        call intaeropt1to3(lchnk, ncol, rhum, mplus10, Nnatk, camnull, &
           bext440n, bext500n, bext550n, bext670n, bext870n,                &
           bebg440n, bebg500n, bebg550n, bebg670n, bebg870n,                &
           bebc440n, bebc500n, bebc550n, bebc670n, bebc870n,                &
           beoc440n, beoc500n, beoc550n, beoc670n, beoc870n,                &
           besu440n, besu500n, besu550n, besu670n, besu870n,                &
           babs440n, babs500n, babs550n, babs670n, babs870n,                &
           bebglt1n, bebggt1n, bebclt1n, bebcgt1n,                         &
           beoclt1n, beocgt1n, bes4lt1n, bes4gt1n,                         &
           backsc550n, babg550n, babc550n, baoc550n, basu550n)

!     and finally the BC&OC(n) mode:
      mplus10=1
        call intaeropt4(lchnk, ncol, rhum, mplus10, Nnatk, camnull, fnbc, faqm4, &
           bext440n, bext500n, bext550n, bext670n, bext870n,                &
           bebg440n, bebg500n, bebg550n, bebg670n, bebg870n,                &
           bebc440n, bebc500n, bebc550n, bebc670n, bebc870n,                &
           beoc440n, beoc500n, beoc550n, beoc670n, beoc870n,                &
           besu440n, besu500n, besu550n, besu670n, besu870n,                &
           babs440n, babs500n, babs550n, babs670n, babs870n,                &
           bebglt1n, bebggt1n, bebclt1n, bebcgt1n,                         &
           beoclt1n, beocgt1n, bes4lt1n, bes4gt1n,                         &
           backsc550n, babg550n, babc550n, baoc550n, basu550n)

!       Dry extinction and absorption optical depth:
        do k=1,pver
          do icol=1,ncol

            do i=0,10
              ec550dry_aer(icol,k)  = ec550dry_aer(icol,k)+Nnatk(icol,k,i)*bext550(icol,k,i)
              abs550dry_aer(icol,k) = abs550dry_aer(icol,k)+Nnatk(icol,k,i)*babs550(icol,k,i)
            enddo
            do i=11,14
              ec550dry_aer(icol,k)  = ec550dry_aer(icol,k)+Nnatk(icol,k,i)*bext550n(icol,k,i-10)
              abs550dry_aer(icol,k) = abs550dry_aer(icol,k)+Nnatk(icol,k,i)*babs550n(icol,k,i-10)
            enddo

            deltah=1.e-4_r8*(pint(icol,k+1)-pint(icol,k))/(rhoda(icol,k)*9.8_r8)
            dod550dry(icol) = dod550dry(icol)+ec550dry_aer(icol,k)*deltah
            abs550dry(icol) = abs550dry(icol)+abs550dry_aer(icol,k)*deltah

            ec550dry_aer(icol,k)  = 1.e-3_r8*ec550dry_aer(icol,k)
            abs550dry_aer(icol,k) = 1.e-3_r8*abs550dry_aer(icol,k)

          enddo
        enddo 

        call outfld('ECDRYAER',ec550dry_aer,pcols,lchnk)
        call outfld('ABSDRYAE',abs550dry_aer,pcols,lchnk)
        call outfld('OD550DRY',dod550dry,pcols,lchnk)
        call outfld('AB550DRY',abs550dry,pcols,lchnk)

#endif  ! ***********AEROCOM***********AEROCOM**************AEROCOM***************above

      endif ! iant==1

!      write(*,*) 'pmxsub 10'

      return
end subroutine pmxsub
