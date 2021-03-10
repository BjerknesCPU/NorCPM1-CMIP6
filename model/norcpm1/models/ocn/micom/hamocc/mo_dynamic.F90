      MODULE mo_dynamic


!***********************************************************************
!
!**** *MODULE mo_dynamic* - Variables for dynamic MLD sampling.
!
!     Patrick Wetzel    *MPI-Met, HH*    28.05.04
!  
!     Purpose
!     -------
!     - declaration and memory allocation
!
!**********************************************************************
      implicit none

        INTEGER, PARAMETER ::                                  &
     &          kdphyto  =1,                                   &
     &          kdgrazer =2,                                   &  
     &          kddoc    =3,                                   &   
     &          kddic    =4,                                   &  
     &          kdphosph =5,                                   &  
     &          kdoxygen =6,                                   &  
     &          kdiron   =7,                                   &
     &          kdano3   =8,                                   &
     &          kdalkali =9,                                   & 
     &          kdsilica =10,                                  &
     &          kdtemp   =11,                                  &
     &          kdsal    =12,                                  &
     &          nbgcdyn  =12

      
        INTEGER, PARAMETER ::                                &
     &          kdadv  =1,				     &
     &          kddif  =2,				     &  
     &          kdpre  =3,				     &  
     &          kdgmp  =4,				     &   
     &          kdbio  =5,				     &  
     &          kdtot  =5

        INTEGER, PARAMETER ::                                &
     &          kdyndiff  =1,				     &
     &          kdynsave  =2

      INTEGER :: nc_dyn_id
      
      INTEGER, DIMENSION (:,:),   ALLOCATABLE :: nbgc_mld
      
      REAL, DIMENSION (:,:),   ALLOCATABLE ::     bgc_zmld
      REAL, DIMENSION (:,:),   ALLOCATABLE ::     bgc_nmld
      
      REAL, DIMENSION (:,:,:,:),   ALLOCATABLE :: bgcdyn
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: bgcdyntmp

      CONTAINS

      SUBROUTINE ALLOC_MEM_DYNAMIC(kpie,kpje,kpke)

      use mod_xc
      use mo_control_bgc
      use mo_bgcmean
      use mo_param1_bgc 
      
      INTEGER :: kpie,kpje,kpke
      INTEGER :: errstat
      
        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcdyn ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',nbgcdyn
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',kdtot
        ENDIF

        ALLOCATE (bgcdyn(kpie,kpje,nbgcdyn,kdtot),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory bgcdyn'
      
        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcdyntmp .'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',nbgcdyn
        ENDIF

        ALLOCATE (bgcdyntmp(kpie,kpje,nbgcdyn),stat=errstat)   
        if(errstat.ne.0) stop 'not enough memory bgcdyntmp'
      
        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgc_zmld ..'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (bgc_zmld(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory bgc_zmld'
        
        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgc_nmld ..'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (bgc_nmld(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory bgc_nmld'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable nbgc_mld ..'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (nbgc_mld(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory nbgc_mld'

      END SUBROUTINE ALLOC_MEM_DYNAMIC

      END MODULE mo_dynamic
