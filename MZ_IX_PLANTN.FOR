C======================================================================
C  MZ_IX_PLANTN, Subroutine
C
C  Handles plant N dynamics
C----------------------------------------------------------------------
C  Revision history
C
C                 Written
C  02/08/1993 PWW Header revision and minor changes    
C  06/01/1994 WTB Modified 
C  06/28/1994 JR/BDB Changed water content dependent factor 
C  03/29/2001 WBD Converted to modular format      
C  12/01/2001 WDB Major revision for 2002 release   
C  03/12/2003 CHP Changed senescence variable to composite (SENESCE)
C                   as defined in ModuleDefs.for
!  03/02/2005 CHP/JIL Add check for negative RNH4U & RNO3U.
!  11/07/2005 CHP Replaced FAC with SOILPROP variable KG2PPM
!  01/31/2006 CHP/JIL Fixed computation of RNLOSS to match root loss
!                     calculation in MZ_GROSUB (from 5% to .5%)
!  03/30/2006 CHP Added composition of senesced matter for SOM modules
!  10/09/2007 JIL Adapted and modified for IXIM model
!  09/09/2012 JIL Modified to handle all plant N dynamics
C----------------------------------------------------------------------
C  Called : MZ_IX_GROSUB 
C
C  Calls  : None
C----------------------------------------------------------------------
!      SUBROUTINE MZ_IX_NUPTAK (DYNAMIC,                              !Control
      SUBROUTINE MZ_IX_PLANTN (DYNAMIC,                               !Control
     &      CARBO,DLAYR,DTT,EAROWT,EARWT,ECNP,GLA,ISTAGE,KG2PPM,LFN,  !Input
     &      LFWTD,LCNP,LFWT,LEAFWT,LEAFWTG,LFCN,LL,NEXP,NH4,NLAYR,NX, !Input
     &      NO3,PGROEA,PGROLF,PGRORT,PGROST,PLIGRT,PLTPOP,PTF,RCNP,   !Input
     &      RLV,RTWT,SAT,SCNP,SENESCE,SHF,STMWT,STOVWT,SW,TCNP,YRDOY, !Input
     &      XSTAGE,ISWNIT,GRAINN,GROGRN,LA,LEAFWTS,LNEXP,NFAC,SEN,    !Input
     &      SENGRT,SLAX,SSLAG,SWFAC,TANCE,TEMPM,TURFAC,VANC,VMNC,     !Input
     &      EARN,LEAFN,RANC,ROOTN,STEMN,STOVN,TANC,TRNU,UNH4,UNO3,    !Output
     &      EANC,LFAN,LFANG,LFANS,NLF,NLFG,NLFS,SUMLFWTG,SLN,SUMNLFG, !Output
     &      SSLA,SSLAS)                                               !Output

      USE ModuleDefs
      IMPLICIT  NONE
      SAVE
C----------------------------------------------------------------------
C                         Define Variables
C----------------------------------------------------------------------

      REAL        ANDEM              
      REAL        CARBO 
      REAL        DLAYR(NL)   
      REAL        DLEAFN            !Daily leaf N
      REAL        DLFN
      REAL        DEARN
      REAL        DROOTN
      REAL        DSTEMN      
      REAL        DSTOVN      
      REAL        DTT
      INTEGER     DYNAMIC
      REAL        EANC
      REAL        EARN
      REAL        EAROWT
      REAL        EARWT
      REAL        ECNP
      REAL        ENDEM
      REAL        GNP
      REAL        KG2PPM(NL)     
      REAL        FACTOR      
      REAL        FNH4        
      REAL        FNO3        
      REAL        GLA(50)
      REAL        GRAINN
      REAL        GROGRN
      INTEGER     I
      INTEGER     ISTAGE
      CHARACTER*1 ISWNIT
      INTEGER     L           
      INTEGER     L1          
      REAL        LA(50)
      REAL        LANC
      REAL        LCNP
      REAL        LEAFGN
      REAL        LEAFN
      REAL        LEAFSN
      REAL        LEAFWT(50)
      REAL        LEAFWTG(50)
      REAL        LEAFWTS(50)
      REAL        LFAN(50)
      REAL        LFANG(50)
      REAL        LFANS(50)
      REAL        LFCN(50)
      REAL        LFDMEXP 
      REAL        LFDMPOOL
      INTEGER     LFN
      REAL        LFND(50)
      REAL        LFNEXP
      REAL        LFNLAB
      REAL        LFNPOOL     !Leaf N being remobilized, g N/leaf
      REAL        LFNRT
      REAL        LFWT
      REAL        LFWTD(50)
      REAL        LL(NL)      
      REAL        LNDEM
      INTEGER     LNEXP
      REAL        NDEM        
      REAL        NEXP
      REAL        NFAC
      REAL        NH4(NL)     
      INTEGER     NLAYR       
      REAL        NLF(50)
      REAL        NLFG(50)
      REAL        NLFS(50)
      REAL        NO3(NL)     
      REAL        NOUT
      REAL        NPOOL
      REAL        NPOOL1
      REAL        NPOOL2
      REAL        NPOOL3
      REAL        NRED
      REAL        NSDR        
      REAL        NSINK
      REAL        NUF         
      REAL        NRF
      REAL        NX
      REAL        PGROEA
      REAL        PGROLF
      REAL        PGRORT      
      REAL        PGROST
      REAL        PLTPOP      
      REAL        PTF         
      REAL        RANC        
      REAL        RCNP        
      REAL        RFAC        
      REAL        RLV(NL)     
      REAL        RMNC
      REAL        RNDEM       
      REAL        RNH4U(NL)   
      REAL        RNLAB
      REAL        RNLOSS      
      REAL        RNO3U(NL)   
      REAL        RNOUT
      REAL        ROOTN       
      REAL        RTWT
      REAL        SANC
      REAL        SEN(50)
      REAL        SENGRT(50)
      REAL        SFAC
      REAL        SHF(NL)     
      REAL        SHOOTN
      REAL        SHOOTWT
      REAL        SLAX
      REAL        SLN(50)
      REAL        SMDFR       
      REAL        SNAV
      REAL        SNDEM
      REAL        SNH4(NL)    
      REAL        SNO3(NL)    
      REAL        STEMN
      REAL        STOVN       
      REAL        STOVWT      
      REAL        SAT(NL)     
      REAL        SCNP
      REAL        SSLA(50)
      REAL        SSLAG(50)
      REAL        SSLAS(50)
      REAL        STMWT
      REAL        SUMLFND
      REAL        SUMLFWT
      REAL        SUMLFWTG
      REAL        SUMNLF
      REAL        SUMNLFG
      REAL        SW(NL)      
      REAL        SWFAC
      REAL        TANC        
      REAL        TANCE
      REAL        TCNP        
      REAL        TEMPM
      REAL        TFAC
      REAL        TNDEM       
      REAL        TNLAB
      REAL        TOTNA
      REAL        TRGTSN
      REAL        TRLV        
      REAL        TRNLOS      
      REAL        TRNU
      REAL        TURFAC        
      REAL        UNH4(NL)    
      REAL        UNO3(NL)    
      REAL        VANC  
      REAL        VMNC
      REAL        XMIN        
      REAL        XNDEM       
      REAL        XNF
      REAL        XNU
      REAL        XSTAGE      
      REAL        XWT
      INTEGER     YRDOY


      TYPE (ResidueType) SENESCE

!----------------------------------------------------------------------
!     CHP 3/30/2006
!     Proportion of lignin in roots
      REAL PLIGRT

C----------------------------------------------------------------------
C     Initialize variables
C----------------------------------------------------------------------
!----------------------------------------------------------------------
!                     DYNAMIC = RUNINIT
!----------------------------------------------------------------------

      IF(DYNAMIC.EQ.RUNINIT.OR.DYNAMIC.EQ.SEASINIT) THEN

      ANDEM    = 0.0
      DLEAFN   = 0.0
      DEARN    = 0.0
	DLFN     = 0.0
      DROOTN   = 0.0
      DSTEMN   = 0.0
      DSTOVN   = 0.0
      EANC     = 0.0
      EARN     = 0.0
      ENDEM    = 0.0
      FACTOR   = 1.0
      FNH4     = 0.0
      FNO3     = 0.0
      GNP      = 0.0
      I        = 0
      L        = 0
      L1       = 0
      LANC     = 0.0
      LEAFGN   = 0.0
      LEAFN    = 0.0
      LEAFSN   = 0.0
      LFDMEXP  = 0.0
      LFDMPOOL = 0.0
      LFNLAB   = 0.0
      LFNPOOL  = 0.0
      LFNRT    = 0.0
      LNDEM    = 0.0
      NDEM     = 0.0
      NOUT     = 0.0
      NPOOL    = 0.0
      NPOOL1   = 0.0
      NPOOL2   = 0.0
      NPOOL3   = 0.0
      NRED     = 0.0
      NRF      = 0.0
      
      NSINK    = 0.0
      NUF      = 0.0
      RFAC     = 0.0
      RMNC     = 0.0
      RNDEM    = 0.0
      RNLAB    = 0.0
      RNLOSS   = 0.0
      RNOUT    = 0.0
!      ROOTN    = 0.0   Initialized in GROSUB
      SANC     = 0.0
      SFAC     = 0.0
      SHOOTN   = 0.0
      SHOOTWT  = 0.0
      SMDFR    = 0.0
      SNAV     = 0.0
      SNDEM    = 0.0
      STEMN    = 0.0
      STOVN    = 0.0
      SUMLFND  = 0.0
      SUMLFWT  = 0.0
      SUMLFWTG = 0.0
      SUMNLF   = 0.0
      SUMNLFG  = 0.0
      TFAC     = 0.0
      TNDEM    = 0.0
      TNLAB    = 0.0
      TOTNA    = 0.0
      TRGTSN   = 0.0
      TRLV     = 0.0
      TRNLOS   = 0.0
      TRNU     = 0.0
      XMIN     = 0.0
      XNDEM    = 0.0
      XNF      = 0.0
      XNU      = 0.0
      XWT      = 0.0
      
	DO I=1,50
	  LFAN(I)     = 0.0
	  LFANG(I)    = 0.0
	  LFANS(I)    = 0.0
	  NLF(I)      = 0.0
	  NLFG(I)     = 0.0
	  NLFS(I)     = 0.0
	  SLN(I)      = 0.0
      END DO
      LFAN(1)  = TANCE
      LFANG(1) = TANCE
      NLF(1)   = LEAFWT(1)*LFAN(1)
      NLFG(1)  = LEAFWTG(1)*LFANG(1)
        RNH4U = 0.0
        RNO3U = 0.0
        SNH4  = 0.0
        SNO3  = 0.0
        UNH4  = 0.0
        UNO3  = 0.0

!-----------------------------------------------------------------------
!                     DYNAMIC = INTEGR
!-----------------------------------------------------------------------
      ELSEIF (DYNAMIC .EQ. INTEGR) THEN   

!-----------------------------------------------------------------------
!     Update actual N concentration in organs g N/plant / g DW/plant
!-----------------------------------------------------------------------

      IF (RTWT.GT.0.0) RANC   = ROOTN / RTWT       ! g N/plant / g DW/plant WHY NOT?? JIL
      IF (STMWT.GT.0.0) SANC  = STEMN / STMWT
      IF (LFWT .GT.0.0) LANC  = LEAFN / LFWT
      IF (EAROWT.GT.0.0) THEN
         EANC = EARN / EAROWT
      ELSE
         EANC = 0.0
      ENDIF

!-----------------------------------------------------------------------
!     Update leaf N concentration, g N / g DW
!-----------------------------------------------------------------------

      DO I = 1, LFN                   !Update leaf N concentrations
        IF (GLA(I) .LE. 0.0) THEN
          NLFS(I)= NLFG(I)+NLFS(I)          
          NLF(I) = NLFS(I)                  
          NLFG(I)= 0.0                      
          LEAFWTS(I) = LEAFWTS(I)+LEAFWTG(I)
          LEAFWT(I) = LEAFWTS(I)            
          LEAFWTG(I) = 0.0                  
          GLA(I) = 0.0                      
          SEN(I) = LA(I)                    
        ENDIF

        NLF(I) = NLFG(I)+NLFS(I)
        IF (LEAFWTG(I) .GT. 0.0) THEN
          LFANG(I) = NLFG(I) / LEAFWTG(I)
        ELSE
          LFANG(I) = 0.0
        ENDIF
        IF (LEAFWTS(I) .GT. 0.0) THEN
          LFANS(I) = NLFS(I) / LEAFWTS(I)
        ELSE
          LFANS(I) = 0.0
        ENDIF
        IF (LEAFWT(I) .GT. 0.0) THEN
          LFAN(I) = NLF(I) / LEAFWT(I)
        ELSE
          LFAN(I) = 0.0
        ENDIF

!---------------------------------------------------------------------------
!     Move leaf N from green to senesced portion due to senescence, g N/leaf
!---------------------------------------------------------------------------
        IF (GLA(I) .GT. 0.0) THEN
!          SLN(I) = NLFG(I)/GLA(I)*10000.
          SLN(I) = NLFG(I)/(GLA(I)+SENGRT(I))   !This SLN only used to move N linked to senescence
!          SLN(I) = AMIN1(SLN(I),3.0)           !SLN: g/cm2
          NOUT = SENGRT(I)*SLN(I)
          NLFG(I) = NLFG(I)-NOUT
          NLFS(I) = NLFS(I)+NOUT
        ELSE
          SLN(I) = 0.0
          NLFS(I)= NLFG(I)+NLFS(I)          !leaf senesce if LFANG(I)<0.8%          
          NLF(I) = NLFS(I)                                                          
          NLFG(I)= 0.0                                                              
        ENDIF

      END DO

!---------------------------------------------------------------------------------
!     Export N and associated DM from expanded green leaves to corresponding pools 
!---------------------------------------------------------------------------------

      LFNPOOL  = 0.0
      LFDMPOOL = 0.0
!      SUMLFWTG = 0.0
      SUMNLFG  = 0.0
      NPOOL3 = 0.0
      DO I = 1, LFN
        NLF(I) = NLFG(I)+NLFS(I)
        IF (LEAFWTG(I) .GT. 0.0) THEN
          LFANG(I) = NLFG(I) / LEAFWTG(I)
        ELSE
          LFANG(I) = 0.0
        ENDIF
        IF (LEAFWTS(I) .GT. 0.0) THEN
          LFANS(I) = NLFS(I) / LEAFWTS(I)
        ELSE
          LFANS(I) = 0.0
        ENDIF
        IF (LEAFWT(I) .GT. 0.0) THEN
          LFAN(I) = NLF(I) / LEAFWT(I)
        ELSE
          LFAN(I) = 0.0
        ENDIF

          LFNEXP  = 0.0
          LFDMEXP = 0.0
          IF(GLA(I) .GT. 0.01) THEN

            IF (I .LE. LNEXP .AND. LFANG(I) .GT. 0.008) THEN    !N export from expanded leaf
                                                     !Reduction in LFANG(I)=0.0000134*DTT
              LFNEXP = LEAFWTG(I)*0.0000134*DTT*NEXP !Target N export from leaf, g N/lf
              IF (NLFG(I)-LFNEXP.LT.0.008*LEAFWTG(I)) THEN
                LFNEXP = NLFG(I)-0.008*LEAFWTG(I)    !Nmin instead of 0.08%??
              ENDIF
              LFDMEXP  = 6.25*LFNEXP                 !Target DM export as protein, g DM/leaf
              IF (LFDMEXP .GT. LEAFWTG(I)) THEN
                LFDMEXP=LEAFWTG(I)
                LFNEXP=LFDMEXP/6.25
              ENDIF
              NLFG(I)=NLFG(I)-LFNEXP
              LEAFWTG(I)=LEAFWTG(I)-LFDMEXP
              NLF(I)=NLFG(I)+NLFS(I)
              LEAFWT(I)=LEAFWTG(I)+LEAFWTS(I)

              LFNPOOL  = LFNPOOL+LFNEXP            !N pool available for mobilization, g N/leaf
              LFDMPOOL = LFDMPOOL+LFDMEXP          !DM pool available for mobilization, g DM/leaf

              IF (LEAFWTG(I) .GT. 0.0) THEN
                LFANG(I)=NLFG(I)/LEAFWTG(I)
              ELSE
                LFANG(I) = 0.0
              ENDIF
              IF (LEAFWT(I) .GT. 0.0) THEN
                LFAN(I)=NLF(I)/LEAFWT(I)
              ELSE
                LFAN(I) = 0.0
              ENDIF
              
            ENDIF  
            IF (LFANG(I) .GT. 0.008) THEN 
              NPOOL3 = NPOOL3+(LFANG(I)-0.008)*LEAFWTG(I)  ! g N in leaves available for mobilization to grain
              SUMNLFG = SUMNLFG+NLFG(I)
!              SUMLFWTG = SUMLFWTG+LEAFWTG(I)
            ENDIF
          ENDIF
      END DO

!---------------------------------------------------------------------------------
!     N mobilization to grain from vegetative organs 
!---------------------------------------------------------------------------------

            IF (ISWNIT.NE.'N' .AND. (ISTAGE.EQ.5.OR.ISTAGE.EQ.6)) THEN !
                                                                       !
           ! Grain N allowed to vary between .01 and .018.             !
           ! High temp., low soil water, and high N increase grain N   !
                                                                       !
              SFAC   = 1.125 - 0.1250*TURFAC                           !
              TFAC   = 0.690 + 0.0125*TEMPM                            !
              GNP    = (0.004+ 0.0130*NFAC)*AMAX1(SFAC,TFAC)           !N concentration in new grain mass
              NSINK  = GROGRN*GNP                                      !Demand grain N, g N/plant
              
              IF (LFNPOOL .GE. NSINK) THEN
                LFNPOOL = LFNPOOL - NSINK                              !
              ELSE
                GRAINN = GRAINN + LFNPOOL
                NSINK = NSINK - LFNPOOL
                LFNPOOL = 0.0
                IF (NSINK .NE. 0.0) THEN                                !
                  RMNC   = 0.75*RCNP                                    !

                  IF (RANC .LT. RMNC) THEN                              !
                    RANC = RMNC                                         !
                  ENDIF                                                 !
                                                      
!                  VANC   = STOVN/STOVWT                                !
!                  VANC   = (SUMNLFG+STEMN+EARN)/(SUMLFWTG+STMWT+EAROWT)  !Current N concentration in shoot
                  VANC   = (STEMN+EARN)/(STMWT+EAROWT)
                  IF (VANC .LT. VMNC) THEN                              !
                    VANC = VMNC                                         !
                  ENDIF                                                 !
                                                                        !
!                  NPOOL1 = STOVWT*(VANC-VMNC)                          !
                  NPOOL1 = (STMWT+EAROWT)*(VANC-VMNC)                   ! g N/plant
                  NPOOL2 = RTWT  *(RANC-RMNC)                           !
                  XNF    = 0.15  + 0.25*NFAC                            !Labile fraction in pools
                  TNLAB  = XNF   * NPOOL1                               !
                  RNLAB  = XNF   * NPOOL2                               !
                  NPOOL  = TNLAB + RNLAB                                !Supply of N, g/plant, available for mobilization
                  IF (NPOOL .LT. NSINK) THEN
                    LFNLAB = XNF * NPOOL3
                    IF (LFNLAB .GT. NSINK-NPOOL) THEN
                      LFNLAB = NSINK-NPOOL
                    ENDIF
!                    NPOOL = NPOOL + LFNLAB
                  ELSE
                    LFNLAB = 0.0
                  ENDIF
                  NPOOL = NPOOL + LFNLAB
!                IF (ISTAGE .EQ. 5 .AND. ICSDUR .EQ. 1) THEN            !

! JIL 09/07/2007 Put a max to GPP based on N supply                    !
!                 GPP  = AMIN1 (GPP*NDEF3,(NPOOL/(0.062*0.0095))) !This should go to KNUMBER JIL
! JIL N stress already taken into account in GPP calculation           !
!                  GPP  = AMIN1 (GPP,(NPOOL/(0.062*0.0095)))            !
                          ! Corrected GPSM calculations .. PWW 2-2-94  !
                                                                       !
!                 GPSM = GPP * EARS                                    !
!                ENDIF                                                  !

                  NSDR   = NPOOL/NSINK                                  !Nitrogen supply to demand ratio
 
                  IF (NSDR .LT. 1.0) THEN                               !
                    NSINK = NSINK*NSDR                                  !Reduce demand if not enough N
                  ENDIF                                                 !
!                  IF (NSINK .GT. TNLAB) THEN                            !Demand > available shoot N??
                  IF (NSINK .GT. TNLAB+RNLAB) THEN                      !Demand > available shoot+root N??
!                    STOVN = STOVN - TNLAB                                !Mobilize N from shoot
!                    SUMNLFG  = 0.0
                    DO I = 1, LFN                                       !Mobilize N from leaves
                      IF (GLA(I).GT.0.01.AND.LFANG(I).GT.0.008) THEN    !poportional to current N content
                        NLFG(I) = NLFG(I)-NLFG(I)*LFNLAB/SUMNLFG
                        NLF(I) =NLFG(I)+NLFS(I)
                      ENDIF
                    END DO
                    STEMN = STEMN - STEMN*TNLAB/(STEMN+EARN)            !And from shoot
                    EARN  = EARN  - EARN *TNLAB/(STEMN+EARN)
!                    RNOUT = NSINK - TNLAB                              !
                    RNOUT = NSINK - TNLAB - LFNLAB
                    ROOTN = ROOTN - RNOUT                               !And root
                    RANC  = ROOTN / RTWT                                !
                  ELSEIF (NSINK .GT. TNLAB) THEN
                    STEMN = STEMN - STEMN*TNLAB/(STEMN+EARN)            !Mobilize N from shoot
                    EARN  = EARN  - EARN *TNLAB/(STEMN+EARN)
                    RNOUT = NSINK - TNLAB                               !
                    ROOTN = ROOTN - RNOUT                               !And root
                    RANC  = ROOTN / RTWT                                !
                  ELSE
!                    STOVN = STOVN - NSINK                                !Mobilize N only from shoot
                    STEMN = STEMN - STEMN*NSINK/(STEMN+EARN)
                    EARN  = EARN  - EARN *NSINK/(STEMN+EARN)
!                    VANC  = STOVN / STOVWT                               !
                    VANC   = (STEMN+EARN)/(STMWT+EAROWT)
                  ENDIF                                                 !
                                                                        !
                ENDIF                                                   !
              ENDIF                                                     !
              GRAINN = GRAINN + NSINK                                   !Update grain N content, g N/plant
              GROGRN = GROGRN + (NSINK * 6.25)                          !N mobilizad as protein JIL
            ENDIF                       !<------------------------------!

C-----------------------------------------------------------------------
C     Calculate potential soil N supply (kg N/ha)
C-----------------------------------------------------------------------
      
      DO L = 1, NLAYR
         RNO3U(L) = 0.0
         RNH4U(L) = 0.0
!         TRLV     = TRLV    + RLV(L)
         TRLV     = TRLV    + RLV(L) * DLAYR(L)
         !KG2PPM(L) = 10.0/(BD(L)*DLAYR(L))
         SNO3(L) = NO3(L) / KG2PPM(L)  ! Converts (mg N/kg soil) into (kg N/ha)
         SNH4(L) = NH4(L) / KG2PPM(L)  ! JIL 
      END DO

C-----------------------------------------------------------------------
C     Calculate N demand (DNG=g N/plant) for new growth (PDWI=g/plt)
C-----------------------------------------------------------------------

!      IF (PDWI .EQ. 0.0) THEN
!         PDWI = 1.0
!      ENDIF
!      DNG = PDWI*TCNP            !JIL 06/20/2013 DNG currently not used
!      IF (XSTAGE .LE. 1.2) THEN
!         DNG = 0.0
!      ENDIF

C-----------------------------------------------------------------------
C     Calculate total N demand (NDEM) of tops (TNDEM) and roots (RNDEM),
C     all expressed in g N/pl. Convert total N demand to kg N/ha (ANDEM)
C     Grain N has been mobilized already thus, TNDEM is calculated from 
C     expanding leaves, stem and ear                                JIL
C-----------------------------------------------------------------------

!      TNDEM  = STOVWT * (TCNP-TANC) + DNG
      RNDEM  = RTWT   * (RCNP-RANC) + PGRORT*RCNP
      STEMN  = STEMN + LFNPOOL
      STMWT  = STMWT + LFDMPOOL
      IF (STMWT .GT. 0.0) SANC = STEMN / STMWT
	SNDEM  = AMAX1(STMWT  * (SCNP-SANC) + PGROST*SCNP, 0.0)     !g N/plant
!	ENDEM  = AMAX1(EAROWT  * (ECNP-EANC) + PGROEA*ECNP, 0.0)    !g N/plant
!	LNDEM  = AMAX1(LFWT   * (LCNP-LANC) + PGROLF*LCNP, 0.0)     !g N/plant
!	SNDEM  = STMWT  * (SCNP-SANC) + PGROST*SCNP                 !g N/plant
	IF (EAROWT .GT. 0.0) THEN
	  EANC   = EARN / EAROWT
	  ENDEM  = AMAX1(EAROWT  * (ECNP-EANC) + PGROEA*ECNP, 0.0)
	ENDIF
	
      SUMLFND = 0.0
      DO I = 1, LFN
        LFND(I) = 0.0
        IF (GLA(I) .GT. 0.0 .AND. LEAFWTG(I) .GT. 0.00001) THEN
          IF (I .GE. LNEXP) THEN           !Expanding leaf: active N sink
            IF (ISTAGE .LE. 3) THEN
!              LFND(I) = AMAX1(LEAFWTG(I)*(LFCN(I)-LFANG(I)) + 
!     &                        LFWTD(I)*LFCN(I),0.0)
!              LFND(I) = LEAFWTG(I)*(LFCN(I)-LFANG(I))+LFWTD(I)*LFCN(I)
              LFND(I) = AMAX1(LEAFWTG(I)*(LCNP-LFANG(I)) +
     &                        LFWTD(I)*LCNP,0.0)
!            ELSEIF (ISTAGE .GT. 3) THEN
!              LFND(I) = AMAX1(LEAFWTG(I)*(LFCN(I)-LFANG(I)),0.0)
            ENDIF
          ENDIF
        ENDIF
        SUMLFND = SUMLFND + LFND(I)
      ENDDO
!      ONDEM  = LNDEM + SNDEM + ENDEM     !Organ N demand; ONDEM=TNDEM?
!JIL  Under water stress reduce shoot N demand?
      SUMLFND = SUMLFND*SWFAC
      SNDEM   = SNDEM*SWFAC
      ENDEM   = ENDEM*SWFAC
!      SUMLFND = SUMLFND
!      SNDEM   = SNDEM
!      ENDEM   = ENDEM

      SHOOTWT = LFWT + STMWT + EARWT
      SHOOTN  = LEAFN + STEMN + EARN + GRAINN
!      TRGTSN = 3.0 *(1.-EXP(-0.0094*SHOOTWT))  !JIL Target shoot N, gN/pl

 !JIL Target shoot N (gN/pl) calculated from Plenet & Lemaire, 2000
      TRGTSN = 0.0832 * SHOOTWT**0.63            !Using critical N
!      TRGTSN = 6.667 * (1.-EXP(-0.0036*SHOOTWT)) !Using max N
      IF (SHOOTN .GT. TRGTSN) THEN              
!        IF (SHOOTN/TRGTSN .LT. 1.25) THEN
          SUMLFND = SUMLFND*(TRGTSN/SHOOTN)
          SNDEM   = SNDEM*(TRGTSN/SHOOTN)
          ENDEM   = ENDEM*(TRGTSN/SHOOTN)
!        ELSE
!          ANDEM = 0.0
!          GOTO 100
!        ENDIF
      ENDIF     

      TNDEM  = SUMLFND + SNDEM + ENDEM    !g N/plt
      
!      NDEM   = ONDEM + RNDEM
      NDEM   = TNDEM  + RNDEM

C-----------------------------------------------------------------------
C     Convert total N demand from g N/plt to kg N/ha (ANDEM)
C-----------------------------------------------------------------------

      ANDEM  = NDEM   * PLTPOP*10.0

C-----------------------------------------------------------------------
C     Calculate potential N supply in soil layers with roots (TRNU)
C     based on soil NO3 & NH4 availability  JIL
C-----------------------------------------------------------------------

      DO L = 1, NLAYR
         IF (RLV(L) .NE. 0.0) THEN
            L1 = L
C
C        The following code was written by JTR to simplify the code for the
C        generic model and to make the functions more like the water uptake
C        functions.  Done on June 28, 1994.
C
C        New water content dependent factor for uptake and new uptake
C        methods that work like the soil water uptake  -- JTR 6/94
C
            SMDFR    = 1.5-6.0*((SW(L)-LL(L))/(SAT(L)-LL(L))-0.5)**2
            SMDFR    = AMAX1 (SMDFR,0.0)
            SMDFR    = AMIN1 (SMDFR,1.0)
            RFAC     = 1.0-EXP(-8.0*RLV(L))
            FNH4     = SHF(L)*0.075  !JIL Was 0.075
            FNO3     = SHF(L)*0.075
            RNH4U(L) = RFAC*SMDFR*FNH4*(NH4(L)-0.5)*DLAYR(L)
            RNO3U(L) = RFAC*SMDFR*FNO3*NO3(L)*DLAYR(L)
            RNH4U(L) = MAX(RNH4U(L),0.0)
            RNO3U(L) = MAX(RNO3U(L),0.0) 
            TRNU     = TRNU + RNO3U(L) + RNH4U(L)     !kg N/ha
         ENDIF
      END DO

C-----------------------------------------------------------------------
C     Setting a maximum N uptake based on the capacity of the system JIL
C-----------------------------------------------------------------------

	XNU = NX *(1.0-EXP(-0.8*CARBO))*PLTPOP*10.0    !kg N/ha
!	XNU = NX *(1.0-EXP(-0.75*CARBO))*PLTPOP*10.0    !kg N/ha

	IF (TRNU .GT. XNU) THEN
	  NRF = XNU/TRNU

        TRNU = 0.0
	  DO L = 1, L1
	    RNO3U(L) = RNO3U(L) * NRF
	    RNH4U(L) = RNH4U(L) * NRF
          TRNU     = TRNU + RNO3U(L) + RNH4U(L)       !kg N/ha
	  END DO
	ENDIF

C-----------------------------------------------------------------------
C     Calculate factor (NUF) to reduce N uptake to level of demand
C-----------------------------------------------------------------------

100   IF (ANDEM .LE. 0.0) THEN
         TRNU  = 0.0
         NUF   = 0.0
       ELSE
         ANDEM = AMIN1 (ANDEM,TRNU)
         IF (TRNU .EQ. 0.0) RETURN
         NUF   = ANDEM/TRNU
         TRNU  = 0.0
      ENDIF
      IF (LNEXP .LE. 3) THEN
        NUF = 1.0
      ENDIF

C-------------------------------------------------------------------------
C     Calculate N uptake in soil layers with roots based on demand (kg/ha)
C-------------------------------------------------------------------------

      DO L = 1, L1
         UNO3(L) = RNO3U(L)*NUF
         UNH4(L) = RNH4U(L)*NUF
         XMIN    = 0.25/KG2PPM(L)  !At least 0.25 mg N/kg soil of NO3, 
	                             !         0.50 mg N/kg soil of NH4 
	                             !must remain in soil layer. 
							       !Original CERES-Maize values 0.1 and 2.0  JIL

         UNO3(L) = MIN (UNO3(L),SNO3(L) - XMIN)
         UNO3(L) = MAX(UNO3(L),0.0)  
         XMIN    = 0.5/KG2PPM(L)
         UNH4(L) = MIN (UNH4(L),SNH4(L) - XMIN)
         UNH4(L) = MAX(UNH4(L),0.0) 
         TRNU    = TRNU + UNO3(L) + UNH4(L)   !kg[N]/ha
      END DO

      IF (PLTPOP .GT. 0.0) THEN
        TRNU = TRNU/(PLTPOP*10.0)             !g[N]/plant
      ELSE
        TRNU = 0.0
      ENDIF
C-----------------------------------------------------------------------
C     Allocate stover and root N
C-----------------------------------------------------------------------

      IF (NDEM .GT. TRNU) THEN
          XNDEM  = TRNU
          FACTOR = XNDEM / NDEM
          NDEM   = XNDEM
          TNDEM  = TNDEM * FACTOR
          SUMLFND = SUMLFND * FACTOR
          SNDEM   = SNDEM * FACTOR
          ENDEM   = ENDEM * FACTOR
!          ONDEM  = ONDEM * FACTOR
          RNDEM  = RNDEM * FACTOR
      ELSE
          FACTOR = 1.0
      ENDIF
         DLEAFN = 0.0
         DSTEMN = 0.0
         DEARN  = 0.0

      IF (NDEM .LE. 0.0 .OR. TRNU .LE. 0.0) THEN
         DSTOVN = 0.0
         DROOTN = 0.0
       ELSE
         !
         ! Calculate root senescence losses @ 0.5%/day
         !
         TRNLOS = 0.0
         DO L = 1, L1
!CHP/JIL 2/1/06 
            !RNLOSS = 0.0
            !IF (TANC .GT. TCNP) THEN
              !RNLOSS = RANC * RTWT * 0.05 * PLTPOP * RLV(L) / TRLV
               RNLOSS = RANC * RTWT * 0.005 * PLTPOP * 
     &                                RLV(L)*DLAYR(L) / TRLV
          !   g[N]/m2 = g[N]/g[root]       * plants/m2
          !                  * g[root]/plant        * fraction
            !ENDIF

            !Calculate N in senesced roots (kg/ha)
            SENESCE % ResE(L,1) = RNLOSS * 10.0
!                      kg[N]/ha =  g/m2  * 10.

            !Back calculate senesced root mass from N senesced.
            IF (RANC .GT. 0.0) THEN
              SENESCE % ResWt(L) = SENESCE % ResE(L,1) / RANC   
            ELSE                               !kg[dry matter]/ha
              SENESCE % ResWt(L) = SENESCE % ResE(L,1) * 10.0 / 0.40   
!              kg[dry matter]/ha =       kg[N]/ha * kg[C]/kg[N]         
!                                                 / kg[C]/kg[dry matter]
            ENDIF
            !Compute lignin, cellulose and carbohydrate portions
            SENESCE % ResLig(L) = SENESCE % ResWt(L) * PLIGRT

            TRNLOS  = TRNLOS + RNLOSS
            !FON(L)  = FON(L) + RNLOSS     !FON is local variable here
         END DO

         ! Adjust DSTOVN and DROOTN to compensate for N lost to FON
         IF(NDEM.GT.0.0.AND.PLTPOP.GT.0.0) THEN
!           DSTOVN = TNDEM / NDEM*TRNU-      PTF*TRNLOS/(PLTPOP*10.0)
!           DROOTN = RNDEM / NDEM*TRNU-(1.0-PTF)*TRNLOS/(PLTPOP*10.0)
!     JIL 06/25/2007
           IF (TNDEM .GT. 0.0) THEN
             DSTOVN = TNDEM / NDEM*TRNU-      PTF*TRNLOS/PLTPOP
!           DSTOVN = ONDEM / NDEM* TRNU-      PTF*TRNLOS/(PLTPOP)
!           DROOTN = RNDEM / NDEM*TRNU-(1.0-PTF)*TRNLOS/(PLTPOP)
!          g N/pl =   fraction  *g N/pl- fraction*g N/m2/(pl/m2) 
             IF (DSTOVN .GT. 0.0) THEN
!               DLEAFN = (LNDEM / ONDEM*DSTOVN)           !g N/plant
!               DLEAFN = (LNDEM / ONDEM*DSTOVN) + LFNPOOL !g N/plant
               DLEAFN = SUMLFND / TNDEM*DSTOVN
               DSTEMN = SNDEM / TNDEM*DSTOVN
               DEARN  = ENDEM / TNDEM*DSTOVN
             ENDIF
           ELSE
             DSTOVN = 0.0
           ENDIF
           DROOTN = RNDEM / NDEM*TRNU-(1.0-PTF)*TRNLOS/(PLTPOP)
         ENDIF
      ENDIF

!      STOVN = STOVN + DSTOVN                   !g N/pl
      ROOTN = ROOTN + DROOTN                   !g N/pl

      IF (DLEAFN .GT. 0.0) THEN                !Allocating per-leaf N
        DO I = 1,LFN
          DLFN = 0.0
          IF (GLA(I) .GT. 0.0 .AND. I .GE. LNEXP) THEN       !Distribute N uptake among expanding leaves
            IF (LEAFWTG(I) .GT. 0.00001 .AND. ISTAGE .LE. 3) THEN
              IF (SUMLFND .GT. 0.0) THEN
                DLFN = LFND(I)/(SUMLFND/FACTOR) * DLEAFN  !g N/leaf
              ENDIF
            ENDIF
          ENDIF
          NLFG(I) = NLFG(I) + DLFN             !g N/leaf
        ENDDO
      ENDIF

C-----------------------------------------------------------------------
C     Before exit update N variables in all leaves and check leaf N
C-----------------------------------------------------------------------
      STEMN = STEMN + DSTEMN                   !g N/pl
      IF(STMWT.GT. 0.0) SANC = STEMN / STMWT   !g N/g DW
      SUMNLF  = 0.0
      SUMNLFG = 0.0
      SUMLFWTG= 0.0
      SUMLFWT = 0.0
      DO I = 1,LFN
        IF (LEAFWTG(I) .GT. 0.0) THEN
          LFANG(I) = NLFG(I) / LEAFWTG(I)
        ELSE
          LFANG(I) = 0.0
        ENDIF

!JIL N fix for leaves with no N 
        IF(LFANG(I).LT.LCNP.AND.I.GE.LNEXP.AND.LFND(I).GT.0.0)THEN
          IF(SANC/SCNP .GT. 0.75) THEN
            NLFG(I)=NLFG(I)+LEAFWTG(I)*(LCNP-LFANG(I))
            STEMN = STEMN-LEAFWTG(I)*(LCNP-LFANG(I))
            IF(STMWT.GT.0.0) SANC=STEMN/STMWT
          ENDIF
        ENDIF    
        NLF(I) = NLFG(I)+NLFS(I)
        IF (LEAFWTG(I) .GT. 0.0) THEN
          LFANG(I) = NLFG(I) / LEAFWTG(I)
        ELSE
          LFANG(I) = 0.0
        ENDIF
        IF (LEAFWTS(I) .GT. 0.0) THEN
          LFANS(I) = NLFS(I) / LEAFWTS(I)
        ELSE
          LFANS(I) = 0.0
        ENDIF
        IF (LEAFWT(I) .GT. 0.0) THEN
          LFAN(I) = NLF(I) / LEAFWT(I)
        ELSE
          LFAN(I) = 0.0
        ENDIF
        
        IF (GLA(I) .GT. 0.0) THEN
          SLN(I) = NLFG(I)/GLA(I)*10000.       !SLN will be used to calculate tomorrow's leaf Assimilation
          SLN(I) = AMIN1(SLN(I),3.0)           !SLN: g/m2
        ELSE
          SLN(I) = 0.0
        ENDIF
        
        SUMNLF  = SUMNLF + NLF(I)            !Check vs LEAFN
        SUMNLFG = SUMNLFG + NLFG(I)          !Used to mobilization to grain during ISTAGE 5
        SUMLFWT = SUMLFWT + LEAFWT(I)
        SUMLFWTG= SUMLFWTG + LEAFWTG(I)
      ENDDO

      LEAFN = SUMNLF
!      STEMN = STEMN + DSTEMN                   !g N/pl (Moved above)
      EARN  = EARN  + DEARN                    !g N/pl
      STOVN = LEAFN + STEMN + EARN 
      LFWT  = SUMLFWT
      STOVWT= LFWT+STMWT+EAROWT

      IF(STOVWT.GT.0.0) TANC = STOVN / STOVWT  !g N/g DW  JIL
      IF(LFWT .GT. 0.0) LANC = LEAFN / LFWT    !g N/g DW
      IF(STMWT.GT. 0.0) SANC = STEMN / STMWT   !g N/g DW
      IF(EAROWT.GT.0.0) THEN
        EANC = EARN / EAROWT                   !g N/g DW
      ELSE
        EANC = 0.0
      ENDIF

      IF(RTWT.GT.0.1*RTWT .AND. RTWT.GT.0.0) THEN
        RANC  = ROOTN / (RTWT-0.01*RTWT)      !g N/g DW
      ENDIF

	   OPEN(UNIT=68,FILE='PerLeaf.DAT',STATUS='UNKNOWN',
     &        ACCESS='APPEND')
 
!	   WRITE (68,1100) XSTAGE,SCNP,SANC
!     &                   LFCN(9),LFCN(14),EANC,LFAN(14)
	   WRITE (68,1100) YRDOY,SSLA(1),SSLA(2),SSLA(3),
     &   SSLA(4),SSLA(5),SSLA(6),SSLA(7),SSLA(8),
     &   SSLA(9),SSLA(10),SSLA(11),SSLA(12),SSLA(13),
     &   SSLA(14),SSLA(15),SSLA(16),SSLA(17),SSLA(18),
     &   SSLA(19),SSLA(20)
         CLOSE (UNIT=68)
!1100   FORMAT (F6.3,2X,F7.3, 2X,F7.3)
1100  FORMAT (I7,X,20(X,F7.3))
      ENDIF       !Endif for DYNAMIC LOOP

      RETURN
!      END SUBROUTINE MZ_IX_NUPTAK
      END SUBROUTINE MZ_IX_PLANTN


C--------------------------------------------------------------------------------------------------
C                         Define Variables
C--------------------------------------------------------------------------------------------------

! ANDEM       !Crop N demand (kg N/ha)
! ANO3        !Total extractable nitrate N in soil profile (kg N/ha)
! ANH4        !Total extractable ammonium N in soil profile (kg N/ha)
! DLAYR(L)    !Soil thickness in layer L (cm)
! DNG         !N demand of potential new growth of tops (g N/plant)
! DROOTN      !Daily change in plant root nitrogen content (g N/plant)
! DSTOVN      !Daily change in plant stover nitrogen content (g N/plant)
! FACTOR      !Ratio of root N uptake to plant N demand
! FNH4        !Unitless soil ammonium supply index
! FNO3        !Unitless soil nitrate supply index
! FON(20)     !Fresh organic nitrogen in soil layer L due to root senescence, kg N/ha
! GNP         !Nitrogen concentration in new grain growth, gN/g dry matter
! GRAINN      !Grain nitrogen content, g N/plant
! KG2PPM(20)  !Factor that converts mg elemental N/kg soil to kg N/ha for soil layer L
! L           !Index counter
! L1          !Lowest soil layer with roots
! LL(20)      !Lower limit of plant extractable water for soil layer L, cm3/cm3
! NDEM        !Plant nitrogen demand, g N/plant
! NH4(20)     !Ammonium in soil layer L, mg elemental N/kg soil
! NLAYR       !Number of soil layer
! NO3(20)     !Nitrate in soil layer L (mg elemental N/kg soil)
! NSDR        !Plant N supply/demand ratio used to modify grain N content
! NSINK       !Grain demand for N during grain filling (g/plant/day)
! NUF         !Plant N supply to demand used to modify N uptake
! PDWI        !Potential increment in new shoot growth, g/plant
! PGRORT      !Potential increment in new root growth, g/plant
! PLTPOP      !Plant population, plants/m2
! PTF         !Ratio of above ground biomass to total biomass
! RANC        !Root actual N concentration, g N/g root
! RCNP        !Root critical nitrogen concentration, g N/g root dry weight
! RFAC        !Interim variable describing the effects of root length density on potential N uptake from a layer
! RLV(20)     !Root length density for soil layer L, cm root/cm2 soil
! RMNC        !Root minimum nitrogen concentration (g N/g root dry weight)
! RNDEM       !Plant root demand for nitrogen (g/plant)
! RNH4U(20)   !Potential ammonia uptake from layer L, kg N/ha
! RNLOSS      !Loss of N from the plant via root exudation in one layer (g N/m2)
! RNO3U(20)   !Potential nitrate uptake from layer L, kg N/ha
! RNOUT       !Nitrogen content in root available for mobilization to grain, g N/plant
! ROOTN       !Root nitrogen content, g N/plant
! RTWT        !Root weight, g/plant
! SFAC        !Drought stress factor for grain nitrogen concentration
! SHF(20)     !Relative root distribution in soil layer L (0-1)
! SMDFR       !Soil moisture deficit factor affecting N uptake
! SNH4(20)    !Ammonium nitrogen in soil layer L, kg N/ha
! SNO3(20)    !Nitrate content in soil layer L, kg N/ha
! STOVN       !Nitrogen content in stover, g N/plant
! STOVWT      !Stover weight (Stem + leaf), g/plant
! SAT(20)     !Saturated water holding capacity for soil layer L, cm3/cm3999
! SUMNLFG     !Total N content in green portion of leaves, gN/leaf
! SW(20)      !Soil water content in layer L, cm3 water/cm3 soil
! TANC        !Nitrogen concentration in above ground biomass, g N/g dry weight
! TCNP        !Critical nitrogen concentration in tops, g N/g dry weight
! TFAC        !Temperature stress factor for grain nitrogen concentration
! TNDEM       !Plant tops demand for nitrogen (g N/plant)
! TRLV        !Total root length density, cm root/cm2 soil
! TRNLOS      !Total plant N lost by root exudation (g N/m2)
! TRNU        !Total potential root nitrogen uptake, kg N/ha converted to g N/pl
! UNH4(20)    !Plant uptake of ammonium from layer (kg N/ha/day)
! UNO3(20)    !Plant uptake of nitrate from a layer (kg N/ha/day)
! VANC        !Plant vegetative actual N concentration, g N/g plant
! XMIN        !Minimum amount of NO3, NH4 remaining in soil layer (kg N/ha)
! XNDEM       !Temporary variable
! XNF         !Modified nitrogen factor based on critical N concentration in vegetative biomass
! XNU         !Maximum allowable daily N uptake (kg N/ha)
! XSTAGE      !Non-integer growth stage indicator

