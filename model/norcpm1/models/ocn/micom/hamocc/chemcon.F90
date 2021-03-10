      SUBROUTINE CHEMCON(modal,kpie,kpje,kpke,psao,ptho,       &
     &           pddpo,pdlxp,pdlyp,ptiestu,kplmonth,omask)

!**********************************************************************
!
!**** *CHEMCON* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - rename: tracer(i,j,k)=ocetra(i,j,k,itracer)
!     - interfacing with ocean model.
!     - use of zonal mean temperature/salinity.
!     akw3(kpje,kpke),akb3(kpje,kpke),ak13(kpje,kpke),ak23(kpje,kpke)
!     aksp(kpje,kpke) again 3d fields! (EMR suggestion)
!
!     Patrick Wetzel,    *MPI-Met, HH*    15.04.02
!     - rename to CHEMCON
!     - chemc(i,j,7) to chemcm(i,j,8,12)
!     
!     
!     Purpose
!     -------
!     .
!
!     Method
!     -------
!     .
!
!     *CALL*       *CHEMCON(modal,kpie,kpje,kpke,psao,ptho,
!                          pddpo,pdlxp,pdlyp,ptiestu,kplmonth)*
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *modal*   - .
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *psao*    - salinity [psu.].
!     *REAL*    *ptho*    - potential temperature [deg C].
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd dimension) [m].
!     *REAL*    *ptiestu* - level depths [m].
!     *INTEGER* *kplmonth*  - month in ocean model
!
!     Externals
!     ---------
!     .
!**********************************************************************

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      use mo_param1_bgc 

      USE mo_control_bgc
      USE mod_xc, only: mnproc

      implicit none

      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: ptho (kpie,kpje,kpke)
      REAL :: psao (kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL :: ptiestu(kpie,kpje,kpke+1)
      REAL :: omask(kpie,kpje)

      REAL :: ZERO,TENM7,SMICR,THOUSI,PERC,FOURTH,THIRD
      REAL :: HALF,ONE,TWO,TEN
      REAL :: devkb,devk1t,devk2t,devkbt,devkst,devks
      REAL :: rgas,bor1,bor2,oxyco
      REAL :: c00,c01,c02,c03,c04,c05,c10,c11,c12,c13,c20,c21,c22,c23
      REAL :: cb0,cb1,cb2,cb3,cw0,cw1,cw2
      REAL :: ox0,ox1,ox2,ox3,ox4,ox5,ox6
      REAL :: an0,an1,an2,an3,an4,an5,an6
      REAL :: akcc1, akcc2, akcc3, akcc4
      REAL :: salchl, temzer, arafra, calfra,sucall,aracal
      REAL :: devk1, devk2, t,q,s,cl
      REAL :: cek0, ckb, ck1, ck2, ckw, oxy, ani
      REAL :: ak1, ak2, akb, akw, ak0, aksp0
      REAL :: rrr,rho2,bor,p,cp,tc
      REAL :: a1,a2,a3,b1,b2,b3,atn2o,rs
        
      INTEGER :: i,j,k,kplmonth,kpie,kpje,kpke,modal,kchem,kmon

!      WRITE(io_stdo_bgc,*)                                             &
!     &    ' CHEMCON is called with modal,kpie,kpje,kpke,kplmonth:  ',  &
!     &                        modal,kpie,kpje,kpke,kplmonth

! That happens when called from INI_BGC
      if (kplmonth .eq. 0) then
         kplmonth = 1
      endif
!
!
!     -----------------------------------------------------------------
!*         1. SET HALF PRECISION CONSTANTS
!             --- ---- --------- ---------
!
      ZERO=0.
      TENM7=10.**(-7.0)
      SMICR=1.E-6
      THOUSI=1./1000.
      PERC=0.01
      FOURTH=0.25
      THIRD=1./3.
      HALF=0.5
      ONE=1.
      TWO=2.
      TEN=10.

!     -----------------------------------------------------------------
!*         3. SET CONVERSION FACTOR SALINITY -> CHLORINITY
!             ------ ------- -- ---- ------- -- ----------
!             (AFTER WOOSTER ET AL., 1969)
!
      SALCHL=1./1.80655
!
!     -----------------------------------------------------------------
!*         4. SET ZERO DEG CENTIGRADE AT KELVIN SCALE
!             --- ---- --- ---------- -- ------ -----
!
       TEMZER=273.16
!
!     -----------------------------------------------------------------
!*         5. SET MEAN TOTAL [CA++] IN SEAWATER (MOLES/KG)
!             (SEE BROECKER A. PENG, 1982, P. 26)
!             ([CA++](MOLES/KG)=1.026E-2*(S/35.) AFTER
!             CULKIN(1965), CF. BROECKER ET AL. 1982)
!             ------------- --- -------- -- --- -----
!
      CALCON=1.03E-2
!
!     -----------------------------------------------------------------
!*         6. SET COEFFICIENTS FOR APPARENT SOLUBILITY EQUILIBRIUM
!             OF CALCITE (INGLE, 1800, EQ. 6)
!             -- ------- ------- ----- --- ----------- -----------
!
       AKCC1=-34.452
       AKCC2=-39.866
       AKCC3=110.21
       AKCC4=-7.5752E-6
!     -----------------------------------------------------------------
!*        6A. SET FRACTION OF ARAGONITE IN BIOGENIC CACO3 PARTICLES
!             --- -------- --------- -- -------- ----- ---------
!
       ARAFRA=0.
!
!     -----------------------------------------------------------------
!*        6B. FRACTION OF CALCITE IN BIOGENIC CACO3 PARTICLES
!             -------- ------- -- -------- ----- ---------
!
       CALFRA=1.-ARAFRA
       SUCALL=ARAFRA+CALFRA

!
!     -----------------------------------------------------------------
!*         7. FACTOR TO GET APPARENT SOLUBILITY PRODUCT FOR
!             ARAGONIT BY MULTIPLICATION WITH APPARENT SOLUBILITY
!             PRODUCT FOR CALCIT (CF. BERNER, 1976,
!             OR BROECKER ET AL., 1982)
!             -- -------- -- ---- ----- -- ------ --- ------- -----
!
      ARACAL=ARAFRA*1.45+CALFRA
!     WRITE(io_stdo_bgc,*) 'ARACAL=',ARACAL
!
!     -----------------------------------------------------------------
!*         8. SET COEFFICIENTS FOR SEAWATER PRESSURE CORRECTION OF

!             (AFTER CULBERSON AND PYTKOWICZ, 1968, CF. BROECKER
!             ET AL., 1982)
!             ------- -------- --- ---------- ----- --- -------- -
!
      DEVK1=24.2
      DEVK2=16.4
      DEVKB=27.5
      DEVK1T=0.085
      DEVK2T=0.04
      DEVKBT=0.095
!
!     -----------------------------------------------------------------
!*         9. SET COEFFICIENTS FOR PRESSURE CORRECTION OF SOLUBILITY
!             PRODUCT OF CACO3 CORRESPONDING TO ARAGONITE/CALCITE
!             RATIO AFTER EDMOND AND GIESKES (1970),
!             P. 1285
!             -- ---- -- --------- ----- ------ --- ------- -------
!
      DEVKST=0.23
      DEVKS=32.8*ARAFRA+35.4*CALFRA
!     WRITE(io_stdo_bgc,*) '***DEVKS=',DEVKS,' DEVKST=',DEVKST,TEN
!
!     -----------------------------------------------------------------
!*        11. SET UNIVERSAL GAS CONSTANT
!             --- --------- --- --------
!
         RGAS=83.143
!     -----------------------------------------------------------------
!*        12. SET BORON CONCENTRATION IN SEA WATER
!*            IN G/KG PER O/OO CL ACCORDING
!             TO RILEY AND SKIRROW, 1965 (P.250)
!             -- ----- --- -------- ---- ------- -
!
      BOR1=0.00023
!
!     -----------------------------------------------------------------
!*        13. SET INVERSE OF ATOMIC WEIGHT OF BORON [G**-1]
!             (USED TO CONVERT SPECIFIC TOTAL BORAT INTO CONCENTRATIONS)
!             ----- -- ------- -------- ----- ----- ---- ---------------
!
      BOR2=1./10.82
!
!     -----------------------------------------------------------------
!*        14. SET INVERS OF NORMAL MOLAL VOLUME OF AN IDEAL GAS
!             [CM**3]
!             ---------- -- ------ ----- ------ -- -- ----- ---
!
      OXYCO=1./22414.4
!
!     -----------------------------------------------------------------
!*        15. SET VOLUMETRIC SOLUBILITY CONSTANTS FOR CO2 IN ML/L
!             WEISS, R. F. (1974)
!             CARBON DIOXIDE IN WATER AND SEAWATER: THE SOLUBILITY OF A
!             NON IDEAL GAS. MARINE CHEMISTRY, VOL. 2, 203-215.
!     -----------------------------------------------------------------

      C00=-58.0931          ! C null null
      C01=90.5069
      C02=22.2940
      C03=0.027766
      C04=-0.025888
      C05=0.0050578
!
!     -----------------------------------------------------------------
!*        16. SET COEFF. FOR 1. DISSOC. OF CARBONIC ACID
!             (EDMOND AND GIESKES, 1970)
!             ------- --- -------- ------- -------- ----
!
      C10=812.27
      C11=3.356
      C12=-0.00171
      C13= 0.000091

!     -----------------------------------------------------------------
!*        17. SET COEFF. FOR 2. DISSOC. OF CARBONIC ACID
!             (EDMOND AND GIESKES, 1970)
!             ------- --- -------- ------- -------- ----
!
      C20=1450.87
      C21=4.604
      C22=-0.00385
      C23= 0.000182
!
!     -----------------------------------------------------------------
!*        18. SET COEFF. FOR 1. DISSOC. OF BORIC ACID
!             (EDMOND AND GIESKES, 1970)
!             ------- --- -------- ------- ----- ----
!
      CB0=2291.90
      CB1=0.01756
      CB2=-3.385
      CB3=-0.32051
!
!     -----------------------------------------------------------------
!*        19. SET COEFF. FOR DISSOC. OF WATER
!             (DICKSON AND RILEY, 1979, EQ. 7, COEFFICIENT
!             CW2 CORRECTED FROM 0.9415 TO 0.09415 AFTER
!             PERS. COMMUN. TO B. BACASTOW, 1988)
!             ----- ------- -- -- --------- ------ -------
!
      CW0=3441.
      CW1=2.241
      CW2=-0.09415
!
!     -----------------------------------------------------------------
!*        20. SET VOLUMETRIC SOLUBILITY CONSTANTS FOR O2 IN ML/L
!             (WEISS, 1970)
!             ------- ------ --------- --------- --- -- -- ----
!
      OX0=-173.4292
      OX1=249.6339
      OX2=143.3483
      OX3=-21.8492
      OX4=-0.033096
      OX5=0.014259
      OX6=-0.0017

!     -----------------------------------------------------------------
!*            SET VOLUMETRIC SOLUBILITY CONSTANTS FOR N2 IN ML/L
!             WEISS, R. F. (1970) THE SOLUBILITY OF NITROGEN
!             OXYGEN AND ARGON IN WATER AND SEAWATER.
!             DEEP-SEA RESEARCH, VOL. 17, 721-735.
!     -----------------------------------------------------------------

       AN0=-172.4965
       AN1=248.4262
       AN2=143.0738
       AN3=-21.7120
       AN4=-0.049781
       AN5=0.025018
       AN6=-0.0034861

!      Constants for laughing gas solubility 
!      (WEISS, 1974, MARINE CHEMISTRY)
!      --------------------------------------  
       a1=-62.7062
       a2=97.3066
       a3=24.1406
       b1=-0.058420
       b2=0.033193
       b3=-0.0051313
       atn2o=3.e-7
       
!
!
!     -----------------------------------------------------------------
!*        21. CHEMICAL CONSTANTS - SURFACE LAYER
!             -------- --------- - ------- -----
!
      IF(modal.LE.0) THEN
      DO 1 j=1,kpje
      DO 1 i=1,kpie
      IF(omask(i,j).GT.0.5) THEN
!      IF(pddpo(i,j,1).GT.0.5) THEN
!
!*        21.1 SET ABSOLUTE TEMPERATURE
!              ------------------------
      T=ptho(i,j,1)+TEMZER
      Q=T*PERC
      S=MAX(25.,psao(i,j,1))

!      Laughing gas solubility (WEISS, 1974)
!      --------------------------------------  
       rs=a1+a2*(100./t)+a3*log(t/100)                  &
     &    +s*( b1 +b2*(t/100.) + b3*(t/100.)**2)

       satn2o(i,j)=atn2o*exp(rs)

!
!*        21.2 CHLORINITY (WOOSTER ET AL., 1969)
!              ---------------------------------
!
      CL=S*SALCHL
!
!*        21.3 LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1974)
!              -------------------------------------------------

      CEK0=C00+C01/Q+C02*ALOG(Q)+S*(C03+C04*Q+C05*Q**2)

!
!*        21.4 PK1, PK2 OF CARB. ACID, PKB OF BORIC ACID 
!              -----------------------------------------
!*             AFTER EDMOND AND GIESKES (1970)
!              -------------------------------

      CKB=CB0/T+CB1*T+CB2+CB3*CL**THIRD
      CK1=C10/T+C11+C12*S*ALOG(T)+C13*S**2
      CK2=C20/T+C21+C22*S*ALOG(T)+C23*S**2

!
!*        21.5 PKW (H2O) (DICKSON AND RILEY, 1979)
!              ------------------------------------

      CKW=CW0/T+CW1+CW2*SQRT(S)

!
!*        19.4 PKW (H2O) LIT. ?
!*****CKW COULD ADDITIONALLY BE EXPRESSED SALIN. DEPENDENT *********
!

!
!*        21.6 LN(K0) OF SOLUBILITY OF O2 (EQ. 4, WEISS, 1970)
!              -----------------------------------------------

      OXY=OX0+OX1/Q+OX2*ALOG(Q)+OX3*Q+S*(OX4+OX5*Q+OX6*Q**2)

!*       SOLUBILITY OF N2 
!        WEISS, R. F. (1970), DEEP-SEA RESEARCH, VOL. 17, 721-735.
!              -----------------------------------------------

      ANI=AN0+AN1/Q+AN2*ALOG(Q)+AN3*Q+S*(AN4+AN5*Q+AN6*Q**2)

!
!*        21.7 K1, K2 OF CARB. ACID, KB OF BORIC ACID (EDMOND AND GIESKES,1970)
!              ----------------------------------------------------------------
      AK1=TEN**(-CK1)
      AK2=TEN**(-CK2)
      AKB=TEN**(-CKB)
!
!*        21.8 IONIC PRODUCT OF WATER KW (H2O) (DICKSON AND RILEY, 1979)
!              ----------------------------------------------------------------
      AKW=TEN**(-CKW)
      AKW3(I,J,1)=AKW
!
!*       21.9 CO2 SOLUBILITY IN SEAWATER (WEISS, 1974, CF. EQ. 12)
!              ----------------------------------------------------------------
      AK0=EXP(CEK0)*SMICR
!
!*       21.10 DENSITY OF SEAWATER AND TOTAL BORATE IN MOLES/L
      RRR=RHO2(S,ptho(i,j,1),ZERO) *THOUSI
      BOR=BOR1*RRR*CL*BOR2
!
!*       21.11 SET CHEMICAL CHEMICAL CONSTANTS
      CHEMCM(i,j,5,kplmonth)=AK0
      CHEMCM(i,j,4,kplmonth)=ak1
      CHEMCM(i,j,3,kplmonth)=ak2
      CHEMCM(i,j,1,kplmonth)=akb
      CHEMCM(i,j,2,kplmonth)=AKW
      CHEMCM(i,j,6,kplmonth)=BOR
!
!*       21.12 O2/N2 SOLUBILITY IN SEAWATER (WEISS, 1970)
!              ----------------------------------------------------------------
      CHEMCM(i,j,7,kplmonth)=EXP(OXY)*OXYCO
      CHEMCM(i,j,8,kplmonth)=EXP(ANI)*OXYCO

!      CHEMC(I,J,7,LAUMON)=EXP(OXY)*OXYCO/196800.
!      CHEMC(I,J,8,LAUMON)=EXP(ani)*OXYCO/802000.

      ENDIF
      
1     CONTINUE

      ENDIF
!Paddy:

      IF(modal.EQ.-13) THEN
         IF(mnproc.eq.1) THEN
           WRITE(io_stdo_bgc,*) 'CHEMCON: initializing all month '
         ENDIF
         DO j=1,kpje
         DO i=1,kpie
         DO kchem=1,8
         DO kmon=1,12
            CHEMCM(i,j,kchem,kmon) = CHEMCM(i,j,kchem,kplmonth)
         ENDDO
         ENDDO
	 ENDDO
	 ENDDO	 
      ENDIF 


!      IF ( kchck .NE. 0 ) THEN
!         DO kchem=1,8
!            WRITE(io_stdo_bgc,*) 'chemcm:k=',kchem
!            CALL EXTR(kpie,kpje,chemcm(1,1,kchem,kplmonth),pddpo(1,1,1),   &
!     &          rmasko,io_stdo_bgc)
!         ENDDO

!       i = kpie/2 - p_ioff
!         IF(i>=1 .AND. i<=kpie) THEN
!           WRITE(io_stdo_bgc,*)                                        &
!     &    'Matrix with chem. constants (CHEMCM) at i=kpie/2; k=1->8:'
!           WRITE(io_stdo_bgc,6003)                                    &
!     &       (j,(chemcm(i,j,k,kplmonth),k=1,8),j=1,kpje)
!         ENDIF
!      ENDIF

!6003  FORMAT(I3.0,'(',I5.0,'m):',7e12.5)

!
!     -----------------------------------------------------------------
!*        22. CHEMICAL CONSTANTS - DEEP OCEAN
!             ----------------------------------------------------------------

      DO 2 k=1,kpke

      DO 2 i=1,kpie
      DO 2 j=1,kpje

!
!*        22.1 APPROX. SEAWATER PRESSURE AT U-POINT DEPTH (BAR)
!              ----------------------------------------------------------------
      P=1.025E-1*ptiestu(i,j,k)
!      WRITE(io_stdo_bgc,*) 'CHEMCON: P=1.025E-1*ptiestu(k)', P,ptiestu(k),k
!
!*        22.1.1 Zonal mean temperature/salinity
!                -------------------------------
!
!      tzsum = 0.0
!      szsum = 0.0
!      vzsum = 0.0
!      DO i=1,kpie
!         IF(omask(i,j).GT.0.5) THEN
!            tzsum = tzsum+ptho(i,j,k)*pdlxp(i,j)*pdlyp(i,j)*pddpo(i,j,k)
!            szsum = szsum+psao(i,j,k)*pdlxp(i,j)*pdlyp(i,j)*pddpo(i,j,k)
!            vzsum = vzsum+            pdlxp(i,j)*pdlyp(i,j)*pddpo(i,j,k)
!         ENDIF
!      ENDDO
!      IF(vzsum.GT.0.5) THEN
!         tzmean = tzsum/vzsum
!         szmean = szsum/vzsum
!      ELSE
!         tzmean = 0.0
!         szmean = 0.0
!      ENDIF
!
!
!*        22.2 SET LIMITS FOR SEAWATER TEMP. AND SALINITY
!              ----------------------------------------------------------------
!              (THIS IS DONE TO AVOID COMPUTATIONAL CRASH AT DRY
!               POINTS DURING CALCULATION OF CHEMICAL CONSTANTS)
!
!*        22.3 SET [H+] (FIRST GUESS)
!              ----------------------------------------------------------------

!
!*        22.4 SET ABSOLUTE TEMPERATURE
!              ----------------------------------------------------------------
      t=ptho(i,j,k)+273.16
      q=t*perc
      s=max(34.,psao(i,j,k))

!
!*        22.5 CHLORINITY (WOOSTER ET AL., 1969)
!              ----------------------------------------------------------------
      CL=S*SALCHL
!
!*        22.6 LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1974)
!              ----------------------------------------------------------------
      CEK0=C00+C01/Q+C02*ALOG(Q)+S*(C03+C04*Q+C05*Q**2)
!
!*        22.7 PK1, PK2 OF CARBONIC ACID, PKB OF BORIC ACID 
!              ----------------------------------------------------------------
!              AFTER EDMOND AND GIESKES (1970)

      CKB=CB0/T+CB1*T+CB2+CB3*CL**THIRD
!
      CK1=C10/T+C11+C12*S*ALOG(T)+C13*S**2
      CK2=C20/T+C21+C22*S*ALOG(T)+C23*S**2

!*        22.8 LN(K0) OF SOLUBILITY OF O2 (EQ. 4, WEISS, 1970)
!              ----------------------------------------------------------------
      OXY=OX0+OX1/Q+OX2*ALOG(Q)+OX3*Q+S*(OX4+OX5*Q+OX6*Q**2)

      satoxy(i,j,k)=exp(oxy)*oxyco
!
!*        22.9 K1, K2 OF CARBONIC ACID, KB OF BORIC ACID, KW (H2O) (LIT.?)
      AK1=TEN**(-CK1)
      AK2=TEN**(-CK2)
      AKB=TEN**(-CKB)
!
!*       22.10 APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE IN SEAWATER
!              ----------------------------------------------------------------
!              (S=27-43, T=2-25 DEG C) AT P=0 (ATMOSPH. PRESSURE)
!              (INGLE, 1800, EQ. 6)

       AKSP0=1.E-7*(AKCC1+AKCC2*S**(THIRD)+AKCC3*ALOG10(S)+AKCC4*T**2)
       CKW=CW0/T+CW1+CW2*SQRT(S)
       AKW3(I,J,K)=TEN**(-CKW)
!
!*       22.11 FORMULA FOR CP AFTER EDMOND AND GIESKES (1970)
!              ----------------------------------------------------------------

!              (REFERENCE TO CULBERSON AND PYTKOQICZ (1968) AS MADE
!              IN BROECKER ET AL. (1982) IS INCORRECT; HERE RGAS IS
!              TAKEN TENFOLD TO CORRECT FOR THE NOTATION OF P IN
!              DBAR INSTEAD OF BAR AND THE EXPRESSION FOR CP IS
!              MULTIPLIED BY LN(10.) TO ALLOW USE OF EXP-FUNCTION
!              WITH BASIS E IN THE FORMULA FOR AKSPP (CF. EDMOND
!              AND GIESKES (1970), P. 1285 AND P. 1286 (THE SMALL
!              FORMULA ON P. 1286 IS RIGHT AND CONSISTENT WITH THE
!              SIGN IN PARTIAL MOLAR VOLUME CHANGE AS SHOWN ON
!              P. 1285))
!      WRITE(io_stdo_bgc,*)  'CHEMCON: CP=P/(RGAS*T)', CP,P,RGAS,T 
      CP=P/(RGAS*T)
!
!*       22.12 KB OF BORIC ACID, K1,K2 OF CARBONIC ACID PRESSURE
!              CORRECTION AFTER CULBERSON AND PYTKOWICZ (1968)
!              (CF. BROECKER ET AL., 1982)

      TC=ptho(i,j,k)

!      WRITE(io_stdo_bgc,*)  ' CHEMCON: modal&(CP*(DEVKB-DEVKBT*TC)) ',modal,CP,DEVKB,DEVKBT,TC
      
      AKB3(I,J,K)=AKB*EXP(CP*(DEVKB-DEVKBT*TC))
      AK13(I,J,K)=AK1*EXP(CP*(DEVK1-DEVK1T*TC))
      AK23(I,J,K)=AK2*EXP(CP*(DEVK2-DEVK2T*TC))
!
!        22.13 APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE (OR ARAGONITE)
!              ----------------------------------------------------------------
!              AS FUNCTION OF PRESSURE FOLLWING EDMOND AND GIESKES (1970)
!              (P. 1285) AND BERNER (1976)
      AKSP(I,J,K)=ARACAL*AKSP0*EXP(CP*(DEVKS-DEVKST*TC))



!
!*       22.14 DENSITY OF SEAWATER AND TOTAL BORATE CONCENTR. [MOLES/L]
!              ----------------------------------------------------------------

      rrr=rho2(s,tc,p)*thousi
      bor=bor1*rrr*cl*bor2
      rrrcl=salchl*1.025*bor1*bor2
!
!     -----------------------------------------------------------------
!*        23. INITIATE [H+] AND [CO3--]
!             -------- ---- --- -------
!
!*       23.1  FIRST GUESSES FOR [CO3--] AND [H+]
!              ----------------------------------------------------------------

!      CARALK=OCETRA(i,j,k,IALKALI)-BOR/(ONE+TENM7/AKB3(I,J,K))

!
!*       20.15 DENSITY OF SEAWATER AND TOTAL BORATE IN MOLES/L
!              ----------------------------------------------------------------
      RRR=RHO2(S,TC,P)*THOUSI
      BOR=BOR1*RRR*CL*BOR2

!
!*       20.16 [CO3--], FIRST ESTIMATE
!              ----------------------------------------------------------------

2     CONTINUE     

!     -----------------------------------------------------------------
!*        21. ITERATION TO INITIATE [H+] AND [CO3--]
!             --------- -- -------- ---- --- -------
!      BT=BOR
!      BREMS=1.5
!      DO 43 KI=1,30
!         zer(ki)=0.
!         DO 43 k=1,kpke
!         DO 43 j=1,kpje
!         DO 43 i=1,kpie
!
!            IF(pddpo(i,j,k).GT.0.5) THEN
!
!            ak1=ak13(i,j,k)
!            ak2=ak23(i,j,k)
!            H=HI(i,j,k)
!            R=CO3(i,j,k)
!            ALKA=OCETRA(i,j,k,IALKALI)
!            C=OCETRA(i,j,k,ISCO212)
!            T1=h/ak13(i,j,k)
!            T2=h/ak23(i,j,k)
!            AKW=AKW3(J,K)
!            BT=rrrcl*psao(i,j,k)
!            AKB=AKB3(J,K)
!            alk=ocetra(i,j,k,ialkali)
!            A=!*(2.+t2)/(1.+t2+t2*t1)  +AKW/H-H+BT/(1.+H/AKB)-ALK
!            A=!*(2.+t2)/(1.+t2+t2*t1)  +AKW/H-H+BT/(1.+H/AKB)-ALK
!            zer(ki)=zer(ki)+a**2
!            DADH=!*(1./(AK2*(1.+T2+T2*T1))-(2.+T2)*(1./AK2+2.*T1/AK2)/
!     1          (1.+T2+T2*T1)**2)
!     1          -AKW/H**2-1.-(BT/AKB)/(1.+H/AKB)**2
!            dddhhh=a/dadh
!            reduk=MAX(1.,1.2*abs(dddhhh/h))
!            H=H-dddhhh/reduk
!  
!            HI(i,j,k)=HI(i,j,k)-dddhhh
!            co3(i,j,k)
!     1      =c/(1.+hi(i,j,k)*(1.+hi(i,j,k)/ak13(i,j,k))/ak23(i,j,k))
!
!           ENDIF
!
!43      CONTINUE

!      WRITE(io_stdo_bgc,*)  ' '
!      WRITE(io_stdo_bgc,*)  'CHEMCON: convergence ', zer

      RETURN
      END
