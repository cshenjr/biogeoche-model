cOPTION CROSSREF ON
cOPTION SYMTABLE ON
cPAGEWIDTH 132
C
C********************************************************************************
C 
C 
      SUBROUTINE   TUNER
C 
C                                           PATHOGEN MODEL - PATH
C                                           ---------------------
C 
C*******************************************************************************
C 
C 
C     SYSTEMS                                                              UNITS
C     -------                                                              -----
C        1 - Salinity (SAL)                                                  PPT
C        2 - Tracer or Total Residual Chlorine (TRC)                        mg/L
C        3 - Total Coliform Bacteria (TCOLI)                             #/100mL
C        4 - Fecal Coliform Bacteria (FCOLI)                             #/100mL
C        5 - Enterococci (ENTERO)                                        #/100mL
C 
C*******************************************************************************
C
C     CONSTANTS
C     --------- 
C                  NAMES AND DESCRIPTIONS OF CONSTANTS
C   NO   NAME      DESCRIPTION                                             UNITS
C 
C    1   KCBC    BASE TEMPERATURE-DEPENDENT COLIFORM BACTERIA MORTALITY     /DAY
C                RATE AT 20 DEG C 
C    2   KCBT    TEMPERATURE COEFFICIENT FOR BACTERIA MORTALITY RATES
C    3   KCSAL   TOTAL/FECAL COLIFORM DIE-OFF RATE DUE TO SEAWATER          /DAY
C    4   KCLGHT  TOTAL COLIFORM DIE-OFF RATE DUE TO SOLAR RADIATION     /LYS-DAY
C    5   SOLAR   ON/OFF SWITCH FOR SOLAR RADIATION DIE-OFF RATE 
C                =0, NO DIE-OFF DUE TO SOLAR RADIATION
C                =1, BACTERIAL DIE-OFF RATE INCLUDES SOLAR RADIATION
C                    EFFECTS
C    6   KEBC    BASE TEMPERATURE-DEPENDENT ENTEROCOCCI BACTERIA MORTALITY
C                RATE AT 20 DEG C                                           /DAY
C    7   KEBT    TEMPERATURE COEFFICIENT FOR ENTEROCOCCI MORTALITY RATES
C    8   KESAL   ENTEROCOCCUS DIE-OFF RATE DUE TO SEAWATER                  /DAY
C    9   KELGHT  ENTEROCOCCUS DIE-OFF RATE DUE TO SOLAR RADIATION       /LYS-DAY
C   10   KEOPT   EXTINCTION COEFFICIENT OPTION 
C                = 0 KE IS A CONSTANT (SPATIALLY AND TEMPORALLY
C                    INVARIANT
C                = 1 KE IS SPATIALLY VARIABLE BUT CONSTANT IN TIME
C                    (USING 2-D PARAMETER ARRAY)
C                = 2 KE IS SPATIALLY INVARIANT BUT VARIES IN TIME
C                    (USING TIME-VARIABLE FUNCTION)
C                = 3 KE IS SPATIALLY VARIABLE AND CAN VARY IN TIME,
C                    (USING THE PRODUCT OF A 2-D PARAMETER ARRAY
C                     AND ONE TIME-VARIABLE FUNCTION)
C                = 4 KE IS SPATIALLY AND TEMPORALLY VARIABLE
C                    (REQUIRES SEPARATE INPUT FILE)
C   11   KECONST BASE EXTINCTION COEFFICIENT (USED WHEN KEOPT=0)             /M
C   12   TRCOPT  TRACER DECAY RATE OPTION
C                = 0 CONSERVATIVE
C                = 1 FIRST-ORDER DECAY
C                = 2 SECOND-ORDER DECAY
C   13   KTRC    DECAY RATE FOR TRACER OR TRC                               /DAY
C               
C 
C*******************************************************************************
C 
C     2-D PARAMETERS
C     --------------
C                       NAMES AND DESCRIPTIONS OF PARAMETERS
C   NO   NAME      DESCRIPTION                                             UNITS
C
C    1   KEBS    EXTINCTION COEFFICIENTS (USED WHEN KEOPT=1,3)                /M
C               
C  
C
C*******************************************************************************
C
C     TIME-VARIABLE FUNCTIONS
C     -----------------------
C
C                       NAMES AND DESCRIPTIONS OF PARAMETERS
C   NO   NAME      DESCRIPTION                                             UNITS
C    1   ITOTSF  TOTAL DAILY SOLAR RADIATION                        LANGLEYS/DAY
C    2   F       FRACTION OF DAYLIGHT
C    3   KETVF   EXTINCTION COEFFICIENT (USED FOR KEOPT=2,3)                  /M
C
C*******************************************************************************
!
!=======================================================================
! The file was revised to write out state variables in NetCDF format   ! 
!                                 YUN LI, UMCES/HPL Feb-18-2011        !
!=======================================================================
!
      USE netcdf
! 
      SAVE

      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
      PARAMETER    (NGDMP=5)
      CHARACTER    GDNAMES(NGDMP)*8,DDNAMES(5,NOSYS)*8

C        STATE-VARIABLES
      REAL
     .      SAL(NX,NY,NZ)    , CDSAL(NX,NY,NZ)
     .    , TRACR(NX,NY,NZ)  , CDTRACR(NX,NY,NZ)
     .    , TCOLI(NX,NY,NZ)  , CDTCOLI(NX,NY,NZ)
     .    , FCOLI(NX,NY,NZ)  , CDFCOLI(NX,NY,NZ)
     .    , ENTERO(NX,NY,NZ) , CDENTERO(NX,NY,NZ)
      REAL
     .      SAL_DDA(NX,NY,NZ)    , TRACR_DDA(NX,NY,NZ)
     .    , TCOLI_DDA(NX,NY,NZ)  , FCOLI_DDA(NX,NY,NZ)
     .    , ENTERO_DDA(NX,NY,NZ) , ECMSAL_DDA(NX,NY,NZ)
      REAL
     .      SAL_GDA(NX,NY,NZ)    , TRACR_GDA(NX,NY,NZ)
     .    , TCOLI_GDA(NX,NY,NZ)  , FCOLI_GDA(NX,NY,NZ)
     .    , ENTERO_GDA(NX,NY,NZ) , ECMSAL_GDA(NX,NY,NZ)

C        LABELED COMMON FOR SEGMENT KES
      COMMON /WCKE/ IKE,NKE,TKE(366)

      REAL
     .    SKE(NX,NY)       , ATTENL(NX,NY,NZ)
     .   ,RIAVG(NX,NY,NZ)  , KCBCT(NX,NY,NZ)
     .   ,KEBCT(NX,NY,NZ)  , KCSALT(NX,NY,NZ)
     .   ,KESALT(NX,NY,NZ) , SKCOLI(NX,NY,NZ)
     .   ,KOCL(NX,NY,NZ)   , SKENTE(NX,NY,NZ)
     .   ,KCL(NX,NY,NZ)    , KEL(NX,NY,NZ)

      REAL
     .    TEMP1(NX,NY)

      REAL
     .    KCBC,KCBT,KCSAL,KCLGHT,SOLAR,KEBC,KEBT,KESAL,KELGHT
     .   ,ITOTSF,ISURF,TEMP,KTRC,KEOPT,KECONST,KETVF
     
      REAL
     .    KCBTT(350),KEBTT(350)

      INTEGER*2
     .    SYSGDP(40)
      INTEGER
     .    ISTP20(NX,NY,NZ)

      EQUIVALENCE
     .  (CARAY(1,1,1,1),SAL(1,1,1))   ,(CDARAY(1,1,1,1),CDSAL(1,1,1))
     . ,(CARAY(1,1,1,2),TRACR(1,1,1)) ,(CDARAY(1,1,1,2),CDTRACR(1,1,1))
     . ,(CARAY(1,1,1,3),TCOLI(1,1,1)) ,(CDARAY(1,1,1,3),CDTCOLI(1,1,1))
     . ,(CARAY(1,1,1,4),FCOLI(1,1,1)) ,(CDARAY(1,1,1,4),CDFCOLI(1,1,1))
     . ,(CARAY(1,1,1,5),ENTERO(1,1,1)),(CDARAY(1,1,1,5),CDENTERO(1,1,1))

      EQUIVALENCE
     .    (CONST(1),KCBC )  , (CONST(2),KCBT )   ,  (CONST(3),KCSAL)
     .  , (CONST(4),KCLGHT) , (CONST(5),SOLAR)   ,  (CONST(6),KEBC)
     .  , (CONST(7),KEBT)   , (CONST(8),KESAL)   ,  (CONST(9),KELGHT)
     .  , (CONST(10),KEOPT) , (CONST(11),KECONST),  (CONST(12),TRCOPT)
     .  , (CONST(13),KTRC)

      REAL
     .    KEBS(NX,NY)
      EQUIVALENCE
     .    (PARAM2D(1,1,1),KEBS(1,1))
!
!-----------------------------------------------------------------------
! Name of this rca file
!-----------------------------------------------------------------------
!
      RCANA='pathogens.f'
!
!----------------------------------------------------------------------
! Create NetCDF file for restart
!----------------------------------------------------------------------
!
      IF(NETCDFOPT.EQ.1 .AND. rstTflag.EQ.0) THEN
        RSTFILNA = TRIM(OUTFILPX)//'_rst.nc'
        CALL nccrt_rcasys(RSTFILNA)
      ENDIF

      IF(INITB.EQ.1)   GO TO 50 
C        SET-UP AND WRITE INFORMATION NEEDED BY GDP
      DO ISYS=1,40
        IF(ISYS.LE.NGDMP)  THEN
          SYSGDP(ISYS)=0
        ELSE
          SYSGDP(ISYS)=1
        ENDIF
      ENDDO
C        FIRST ... NAMES OF THE GLOBAL DUMP VARIABLES
      GDNAMES( 1) = 'Salinity'
      GDNAMES( 2) = 'TRC     '
      GDNAMES( 3) = 'T Coli  '
      GDNAMES( 4) = 'F Coli  '
      GDNAMES( 5) = 'Entero  '

!
!-------------------------------------------------------
! Optionally rewrite the RCAF10 FILE
!-------------------------------------------------------
!
      IF(NETCDFOPT.EQ.1) THEN
      ELSE
C               REWRITE RCAF10 FILE
       REWIND(10)
       WRITE(10)   NX,NY,NZ,NOSYS,NGDMP
       WRITE(10)   GDNAMES
       WRITE(10)   SYSGDP
       WRITE(10)   FSM
      ENDIF

C        WRITE DDNAMES TO RCAF12
      DDNAMES(1,1) = 'Salinity'
      DDNAMES(2,1) = 'HydroSal'
      DDNAMES(3,1) = 'HydroTmp'
      DDNAMES(4,1) = 'Dummy   '
      DDNAMES(5,1) = 'Dummy   '
      DDNAMES(1,2) = 'TRC     '
      DDNAMES(2,2) = 'Dummy   '
      DDNAMES(3,2) = 'Dummy   '
      DDNAMES(4,2) = 'Dummy   '
      DDNAMES(5,2) = 'Dummy   '
      DDNAMES(1,3) = 'T Coli  '
      DDNAMES(2,3) = 'KCbase  '
      DDNAMES(3,3) = 'KCsalt  '
      DDNAMES(4,3) = 'KClight '
      DDNAMES(5,3) = 'KCtotal '
      DDNAMES(1,4) = 'AvgLght '
      DDNAMES(2,4) = 'Dummy   '
      DDNAMES(3,4) = 'Dummy   '
      DDNAMES(4,4) = 'Dummy   '
      DDNAMES(5,4) = 'Dummy   '
      DDNAMES(1,5) = 'Entero  '
      DDNAMES(2,5) = 'KEbase  '
      DDNAMES(3,5) = 'KEsalt  '
      DDNAMES(4,5) = 'KElight '
      DDNAMES(5,5) = 'KEtotal '
      WRITE(12)  DDNAMES

      IF(SOLAR.EQ.0.)  THEN
         KCLGHT=0.
         KELGHT=0.
      ELSE
C        SET UP EXTINCTION COEFFICIENTS BASED ON INPUT OPTION SELECTED
       IF(KEOPT.EQ.0.)  THEN
         DO IY=1,NY
          DO IX=1,NX
           IF(FSM(IX,IY).EQ.1.) THEN
            SKE(IX,IY) = KECONST
           ELSE
            SKE(IX,IY) = 0.
           ENDIF
          ENDDO
         ENDDO
       ELSEIF(KEOPT.EQ.1.)  THEN
         DO IY=1,NY
          DO IX=1,NX
           IF(FSM(IX,IY).EQ.1.) THEN
            SKE(IX,IY) = KEBS(IX,IY)
           ELSE
            SKE(IX,IY) = 0.
           ENDIF
          ENDDO
         ENDDO
       ELSEIF(KEOPT.EQ.4.) THEN
C        READ TIME-VARIABLE EXTINCTION COEFFICIENTS
         OPEN(28,FILE=KINFILNA(1),FORM='FORMATTED',ERR=905)
C        READ TIME INTERVALS TO UPDATE EXTINCTION COEFFICIENTS (DAYS)
         READ(28,1200,ERR=910)  NKE,(TKE(I),I=1,NKE)
 1200    FORMAT(I5/(8F10.0))
         READ(28,1300,ERR=915)  ((SKE(IX,IY),IX=1,NX),IY=1,NY)
 1300    FORMAT(8F10.0)
         WRITE(OUT,1155) TKE(IKE)
 1155    FORMAT(//30X,'SEGMENT KEs (/METER) AT TIME = ',F5.0/
     .   6X,'X-->',25X,2(24X,'1'),24X,'2'/4X,'Y',5X,2(24X,'5',24X,'0')/)
         DO IY=1,NY
         WRITE(OUT,1156)  IY,(SKE(IX,IY),IX=1,NX)
 1156    FORMAT(I5,5X,20F5.1/(6X,20F5.1))
         ENDDO
         IKE=2
       ENDIF
      ENDIF

      DO 35 ITEMP=1,350
        TEMP = (FLOAT(ITEMP-1)/10.+0.05) - 20.
        KCBTT(ITEMP) = KCBT**TEMP
        KEBTT(ITEMP) = KEBT**TEMP
   35 CONTINUE

C        INITIALIZE ARRAY FOR GLOBAL DUMP AVERAGING, IF REQUIRED
      IF(IGDOPT.EQ.1)  THEN
c$doacross local(iz,iy,ix) , share(sal_gda,dye1_gda,bac1_gda,
c$&        bac2_gda,bac3_gda,ecmsal_gda)
        DO 42 IZ=1,NZ
        DO 42 IY=1,NY
        DO 42 IX=1,NX
          SAL_GDA(IX,IY,IZ) = 0.
          TRACR_GDA(IX,IY,IZ) = 0.
          TCOLI_GDA(IX,IY,IZ) = 0.
          FCOLI_GDA(IX,IY,IZ) = 0.
          ENTERO_GDA(IX,IY,IZ) = 0.
          ECMSAL_GDA(IX,IY,IZ) = 0.
   42   CONTINUE
      IAVGGDCNTR = 0
      ENDIF

C        INITIALIZE ARRAY FOR DETAILED DUMP AVERAGING, IF REQUIRED
      IF(IDDOPT.EQ.1)  THEN
c$doacross local(iz,iy,ix) , share(sal_dda,dye1_dda,bac1_dda,
c$&        bac2_dda,bac3_dda,ecmsal_dda)
        DO 44 IZ=1,NZ
        DO 44 IY=1,NY
        DO 44 IX=1,NX
          SAL_DDA(IX,IY,IZ) = 0.
          TRACR_DDA(IX,IY,IZ) = 0.
          TCOLI_DDA(IX,IY,IZ) = 0.
          FCOLI_DDA(IX,IY,IZ) = 0.
          ENTERO_DDA(IX,IY,IZ) = 0.
          ECMSAL_DDA(IX,IY,IZ) = 0.
   44   CONTINUE
      IAVGDDCNTR = 0
      ENDIF


C931206  INITIALIZE
      DO 48 IZ=1,NZ
       DO 48 IY=1,NY
        DO 48 IX=1,NX
         CARAY(IX,IY,IZ,1) = HYDSAL(IX,IY,IZ)
C        CARAY(IX,IY,IZ,2) = 0.0
C        CARAY(IX,IY,IZ,3) = 0.0
   48 CONTINUE

C 
C                  S Y S T E M   L O O P S
C 

   50 CONTINUE
      IF(SOLAR.EQ.1.)  THEN
       IF(NXCALL13T.LE.TIME)   THEN
         CALL RCA13(NXCALL13T)
       ENDIF
       ITOTSF = (TIME-NXFUNT(1))*MFUNC(1) + BFUNC(1)
       F = (TIME-NXFUNT(2))*MFUNC(2) + BFUNC(2)
       KETVF = (TIME-NXFUNT(3))*MFUNC(3) + BFUNC(3)
       TSUNRISE = 0.5-F/2
       TSUNSET = 0.5+F/2
       TDAY = MOD(TIME,1.0)
       IF(TDAY.LT.TSUNRISE.OR.TDAY.GT.TSUNSET) THEN
        ISURF = 0.0
       ELSE
        ISURF = ITOTSF/F/0.635*SIN(3.14159*(TDAY-TSUNRISE)/F)
       ENDIF
C        UPDATE SEGMENT KEs AS APPROPRIATE
       IF(KEOPT.EQ.2.)  THEN
         DO IY=1,NY
          DO IX=1,NX
           IF(FSM(IX,IY).EQ.1.)  THEN
            SKE(IX,IY) = KETVF
           ELSE
            SKE(IX,IY) = 0.0
           ENDIF
          ENDDO
         ENDDO
       ELSEIF(KEOPT.EQ.3.)  THEN
         DO IY=1,NY
          DO IX=1,NX
           IF(FSM(IX,IY).EQ.1.)  SKE(IX,IY) = KETVF*KEBS(IX,IY)
          ENDDO
         ENDDO
       ENDIF
       IF(KEOPT.EQ.4.)  THEN
        IF(TKE(IKE).LE.TIME)  THEN
          READ(28,1300,ERR=915)  ((SKE(IX,IY),IX=1,NX),IY=1,NY)
          WRITE(OUT,1155) TKE(ITEMP)
          DO IY=1,NY
           WRITE(OUT,1156)  IY,(SKE(IX,IY),IX=1,NX)
          ENDDO
          IKE=IKE+1
        ENDIF
       ENDIF
      ENDIF

      DO 60 IZ=1,NZ
      DO 60 IY=1,NY
      DO 60 IX=1,NX
       IF(FSM(IX,IY).NE.1.)  GO TO 60
        IF(HYDTEMP(IX,IY,IZ).GE.0.0 .AND. HYDTEMP(IX,IY,IZ).LE.35.) THEN
         ISTP20(IX,IY,IZ) = INT(10.*HYDTEMP(IX,IY,IZ)+1)
        ELSEIF(HYDTEMP(IX,IY,IZ).LT.0) THEN
         ISTP20(IX,IY,IZ) = 1
        ELSEIF(HYDTEMP(IX,IY,IZ).GT.35.) THEN
         ISTP20(IX,IY,IZ) = 350
        ENDIF
   60 CONTINUE

C        LOOP FOR DETAILED DUMP AVERAGING, IF REQUIRED
      IF(IDDOPT.EQ.1)  THEN
c$doacross local(iz,iy,ix) , share(sal_dda,ecmsal_dda,dye1_dda,bac1_dda,
c$&        bac2_dda,bac3_dda)
        DO 70 IZ=1,NZ
         DO 70 IY=1,NY
          DO 70 IX=1,NX
            SAL_DDA(IX,IY,IZ)  = SAL_DDA(IX,IY,IZ) + SAL(IX,IY,IZ)
            ECMSAL_DDA(IX,IY,IZ)=ECMSAL_DDA(IX,IY,IZ)+HYDSAL(IX,IY,IZ)
            TRACR_DDA(IX,IY,IZ) = TRACR_DDA(IX,IY,IZ) + TRACR(IX,IY,IZ)
            TCOLI_DDA(IX,IY,IZ) = TCOLI_DDA(IX,IY,IZ) + TCOLI(IX,IY,IZ)
            FCOLI_DDA(IX,IY,IZ) = FCOLI_DDA(IX,IY,IZ) + FCOLI(IX,IY,IZ)
            ENTERO_DDA(IX,IY,IZ) = ENTERO_DDA(IX,IY,IZ)+ENTERO(IX,IY,IZ)
   70   CONTINUE
        IAVGDDCNTR = IAVGDDCNTR + 1
      ENDIF

      IF(SOLAR.EQ.1.)  THEN
C     GET LIGHT ATTENUATION COEFFICIENTS
       DO 75 IY=1,NY
        DO 75 IX=1,NX
         IF(FSM(IX,IY).LE.0.)  GO TO 75
         DO 74 IZ=1,NZ
          IF(IZ.EQ.1)   THEN
             ITOT=ISURF
          ELSE
             ITOT=ATTENL(IX,IY,IZ-1)
          ENDIF
          ATTENL(IX,IY,IZ) = ITOT*EXP(-SKE(IX,IY)*HBAR(IX,IY)*DZ(IZ))
  74    CONTINUE
  75   CONTINUE
      ENDIF

C        SYSTEM LOOPS
C        SALINITY and TRC/TRACER

      DO 80 IZ=1,NZ
       DO 80 IY=1,NY
        DO 80 IX=1,NX
         CDSAL(IX,IY,IZ) = 0.0
         CDTRACR(IX,IY,IZ) = 0.0
   80 CONTINUE
      IF(TRCOPT.EQ.0) GO TO 95
      IF(TRCOPT.EQ.1) THEN
        DO 85 IZ=1,NZ
         DO 85 IY=1,NY
          DO 85 IX=1,NX
           IF(FSM(IX,IY).LE.0.)  GO TO 85
           CDTRACR(IX,IY,IZ) = -KTRC*TRACR(IX,IY,IZ)
   85   CONTINUE
      ELSEIF(TRCOPT.EQ.2) THEN
        DO 90 IZ=1,NZ
         DO 90 IY=1,NY
          DO 90 IX=1,NX
           IF(FSM(IX,IY).LE.0.)  GO TO 90
           CDTRACR(IX,IY,IZ) = -KTRC*TRACR(IX,IY,IZ)*TRACR(IX,IY,IZ)
   90   CONTINUE
      ENDIF

C        THREE TYPEs OF BACTERIA
   95 DO 100 IZ=1,NZ
       DO 100 IY=1,NY
        DO 100 IX=1,NX
         IF(FSM(IX,IY).LE.0.)  GO TO 100
C        BASE BACTERIAL DECAY RATES
        KCBCT(IX,IY,IZ) = KCBC*KCBTT(ISTP20(IX,IY,IZ))
        KEBCT(IX,IY,IZ) = KEBC*KEBTT(ISTP20(IX,IY,IZ))
C        SALINITY INDUCED BACTERIAL DECAY RATES
C        Note: 100 percent seawater = 35.ppt salinity (2.857143=100/35)
        KCSALT(IX,IY,IZ) = KCSAL*2.857143*HYDSAL(IX,IY,IZ)
     .                     *KCBTT(ISTP20(IX,IY,IZ))
        KESALT(IX,IY,IZ) = KESAL*2.857143*HYDSAL(IX,IY,IZ)
     .                     *KEBTT(ISTP20(IX,IY,IZ))
C        LIGHT INDUCED BACTERIAL DECAY RATES
         IF(SOLAR.EQ.0.)  THEN
           KCL(IX,IY,IZ)=0.0
           KEL(IX,IY,IZ)=0.0
         ELSE
C        LIGHT ATTENUATION
          IF(IZ.EQ.1)   THEN
           ITOT=ISURF
          ELSE
           ITOT=ATTENL(IX,IY,IZ-1)
          ENDIF
C        COMPUTE AVGERAGE LIGHT
          TEMP1(IX,IY) = SKE(IX,IY)*HBAR(IX,IY)*DZ(IZ)
          RIAVG(IX,IY,IZ) =ITOT/TEMP1(IX,IY)*(1.0 - EXP(-TEMP1(IX,IY)))
          KCL(IX,IY,IZ) = KCLGHT*RIAVG(IX,IY,IZ)
          KEL(IX,IY,IZ) = KELGHT*RIAVG(IX,IY,IZ)
         ENDIF
        SKCOLI(IX,IY,IZ) = KCBCT(IX,IY,IZ) + KCSALT(IX,IY,IZ)
     .                   + KCL(IX,IY,IZ)
        SKENTE(IX,IY,IZ) = KEBCT(IX,IY,IZ) + KESALT(IX,IY,IZ)
     .                   + KEL(IX,IY,IZ)

        CDTCOLI(IX,IY,IZ) = -SKCOLI(IX,IY,IZ)*TCOLI(IX,IY,IZ)
        CDFCOLI(IX,IY,IZ) = -SKCOLI(IX,IY,IZ)*FCOLI(IX,IY,IZ)
        CDENTERO(IX,IY,IZ) = -SKENTE(IX,IY,IZ)*ENTERO(IX,IY,IZ)
  100 CONTINUE

C        PUT ECOM COMPUTED SALINITY AT BOUNDARY
      DO 115 I=1,NOBCALL
        IX = IBCALL(1,I)
        IY = IBCALL(2,I)
        IZ = IBCALL(3,I)
        CARAY(IX,IY,IZ,1) = HYDSAL(IX,IY,IZ)
  115 CONTINUE

      IF(IDISK.EQ.2 .OR. IDISK.EQ.3)  THEN
       DO 131 IDMP=1,NDMPS
         IX = IFDMPS(IDMP,1)
         IY = IFDMPS(IDMP,2)
         IZ = IFDMPS(IDMP,3)
         IF(IDDOPT.EQ.0)  THEN
          CALL RCAWBUF(1,SAL(IX,IY,IZ),HYDSAL(IX,IY,IZ),
     .         HYDTEMP(IX,IY,IZ),DUMMY,DUMMY)
         ELSE
          CALL RCAWBUF(1,SAL_DDA(IX,IY,IZ)/IAVGDDCNTR,
     .         ECMSAL_DDA(IX,IY,IZ)/IAVGDDCNTR,HYDTEMP(IX,IY,IZ),
     .         DUMMY,DUMMY)
         ENDIF
  131  CONTINUE

       DO 132 IDMP=1,NDMPS
         IX = IFDMPS(IDMP,1)
         IY = IFDMPS(IDMP,2)
         IZ = IFDMPS(IDMP,3)
         IF(IDDOPT.EQ.0)  THEN
          CALL RCAWBUF(2,TRACR(IX,IY,IZ),CDTRACR(IX,IY,IZ),
     .         DUMMY,DUMMY,DUMMY)
         ELSE
          CALL RCAWBUF(2,TRACR_DDA(IX,IY,IZ)/IAVGDDCNTR,
     .         CDTRACR(IX,IY,IZ),DUMMY,DUMMY,DUMMY)
         ENDIF
  132  CONTINUE

        DO 145 IDMP=1,NDMPS
         IX = IFDMPS(IDMP,1)
         IY = IFDMPS(IDMP,2)
         IZ = IFDMPS(IDMP,3)
         IF(IDDOPT.EQ.0)  THEN
          CALL RCAWBUF(3,TCOLI(IX,IY,IZ), KCBCT(IX,IY,IZ),
     .         KCSALT(IX,IY,IZ),KCL(IX,IY,IZ),SKCOLI(IX,IY,IZ))
         ELSE
          CALL RCAWBUF(3,TCOLI_DDA(IX,IY,IZ)/IAVGDDCNTR,KCBCT(IX,IY,IZ),
     .         KCSALT(IX,IY,IZ),KCL(IX,IY,IZ),SKCOLI(IX,IY,IZ))
         ENDIF
  145   CONTINUE

        DO 146 IDMP=1,NDMPS
         IX = IFDMPS(IDMP,1)
         IY = IFDMPS(IDMP,2)
         IZ = IFDMPS(IDMP,3)
         IF(IDDOPT.EQ.0)  THEN
          CALL RCAWBUF(4,FCOLI(IX,IY,IZ),RIAVG(IX,IY,IZ),
     .         DUMMY,DUMMY,DUMMY)
         ELSE
          CALL RCAWBUF(4,FCOLI_DDA(IX,IY,IZ)/IAVGDDCNTR,
     .         RIAVG(IX,IY,IZ),DUMMY,DUMMY,DUMMY)
         ENDIF
  146   CONTINUE

        DO 147 IDMP=1,NDMPS
         IX = IFDMPS(IDMP,1)
         IY = IFDMPS(IDMP,2)
         IZ = IFDMPS(IDMP,3)
         IF(IDDOPT.EQ.0)  THEN
          CALL RCAWBUF(5,ENTERO(IX,IY,IZ),KEBCT(IX,IY,IZ),
     .         KESALT(IX,IY,IZ),KEL(IX,IY,IZ),SKENTE(IX,IY,IZ))
         ELSE
          CALL RCAWBUF(5,ENTERO_DDA(IX,IY,IZ)/IAVGDDCNTR,
     .         KEBCT(IX,IY,IZ),KESALT(IX,IY,IZ),KEL(IX,IY,IZ),
     .         SKENTE(IX,IY,IZ))
         ENDIF
  147   CONTINUE
      ENDIF

C        CLOSE BUFFER

      IF(IDISK.EQ.2 .OR. IDISK.EQ.3)  THEN
         CALL RCAWRIT

C        RE-INITIALIZE ARRAY FOR DETAILED DUMP AVERAGING, IF REQUIRED
        IF(IDDOPT.EQ.1)  THEN
c$doacross local(iz,iy,ix) , share(sal_dda,ecmsal_dda,dye1_dda,bac1_dda,
c$&        bac2_dda,bac3_dda)
          DO 210 IZ=1,NZ
           DO 210 IY=1,NY
            DO 210 IX=1,NX
              SAL_DDA(IX,IY,IZ) = 0.
              ECMSAL_DDA(IX,IY,IZ) = 0.
              TRACR_DDA(IX,IY,IZ) = 0.
              TCOLI_DDA(IX,IY,IZ) = 0.
              FCOLI_DDA(IX,IY,IZ) = 0.
              ENTERO_DDA(IX,IY,IZ) = 0.
  210       CONTINUE
          IAVGDDCNTR = 0
        ENDIF
      ENDIF

C        PERFORM GLOBAL DUMP AVERAGING, IF REQUIRED
      IF(IGDOPT.EQ.1)  THEN
c$doacross local(iz,iy,ix) , share(sal_gda,ecmsal_gda,dye1_gda,bac1_gda,
c$&        bac2_gda,bac3_gda)
        DO 220 IZ=1,NZ
         DO 220 IY=1,NY
          DO 220 IX=1,NX
            SAL_GDA(IX,IY,IZ) = SAL_GDA(IX,IY,IZ) + SAL(IX,IY,IZ)
            ECMSAL_GDA(IX,IY,IZ)=ECMSAL_GDA(IX,IY,IZ)+HYDSAL(IX,IY,IZ)
            TRACR_GDA(IX,IY,IZ) = TRACR_GDA(IX,IY,IZ) + TRACR(IX,IY,IZ)
            TCOLI_GDA(IX,IY,IZ) = TCOLI_GDA(IX,IY,IZ) + TCOLI(IX,IY,IZ)
            FCOLI_GDA(IX,IY,IZ) = FCOLI_GDA(IX,IY,IZ) + FCOLI(IX,IY,IZ)
            ENTERO_GDA(IX,IY,IZ) = ENTERO_GDA(IX,IY,IZ)+ENTERO(IX,IY,IZ)
  220   CONTINUE
  225   CONTINUE
      ENDIF

      IAVGGDCNTR = IAVGGDCNTR + 1


C        CONVERT TO MASS UNITS

c$doacross local(iz,iy,ix,isys)
      DO 230 ISYS=1,NOSYS
       DO 230 IZ=1,NZ
        DO 230 IY=1,NY
         DO 230 IX=1,NX
          CDARAY(IX,IY,IZ,ISYS) = BVOL(IX,IY,IZ)*CDARAY(IX,IY,IZ,ISYS)
  230 CONTINUE


  260 CONTINUE

  
C         CHECK IF TIME TO DUMP TO DISK
      IF(IDISK.EQ.0)   RETURN

C        GLOBAL DUMPS

      IF(IDISK.EQ.1 .OR. IDISK.EQ.3)  THEN
       IF(IGDOPT.EQ.0)  THEN
!
!------------------------------------------------------------
! Write out state variable into NetCDF files 
!------------------------------------------------------------
!
        IF(NETCDFOPT.EQ.1) THEN
         WRITE(OUT,'(A,i5.5,A,f15.5,A,i5.5,A,A)'),
     .           'WRITE: TOTAL IREC = ',IREC,' Time = ',time,
     .           ' NetCDF REC = ', ncIREC, ' in ',
     .           TRIM(ADJUSTL(OUTFILNA))
         status=nf90_open(TRIM(ADJUSTL(OUTFILNA)),nf90_write,ncID)
         CALL nccheck_status(status,OUTFILNA,RCANA)
         CALL ncwrt_t1dvar(ncID,idpath_TIME,TIME)
         CALL ncwrt_g2dvar(ncID,idpath_ETA ,ETA )
         CALL ncwrt_g3dvar(ncID,idpath_SAL ,SAL )
         CALL ncwrt_g3dvar(ncID,idTRACR    ,TRACR)
         CALL ncwrt_g3dvar(ncID,idTCOLI    ,TCOLI)
         CALL ncwrt_g3dvar(ncID,idFCOLI    ,FCOLI)
         CALL ncwrt_g3dvar(ncID,idENTERO   ,ENTERO)
         status=nf90_close(ncID)
         CALL nccheck_status(status,OUTFILNA,RCANA)
         ncIREC=ncIREC+1
        ELSE
         WRITE(10)  TIME
         WRITE(11)   SAL
         WRITE(11)   TRACR
         WRITE(11)   TCOLI
         WRITE(11)   FCOLI
         WRITE(11)   ENTERO
        ENDIF
       ELSE
c$doacross local(iz,iy,ix) , share(sal_gda,ecmsal_gda,dye1_gda,bac1_gda,
c$&        bac2_gda,bac3_gda)
         DO 360 IZ=1,NZ
          DO 360 IY=1,NY
           DO 360 IX=1,NX
            SAL_GDA(IX,IY,IZ) = SAL_GDA(IX,IY,IZ)/FLOAT(IAVGGDCNTR)
            ECMSAL_GDA(IX,IY,IZ)= ECMSAL_GDA(IX,IY,IZ)/FLOAT(IAVGGDCNTR)
            TRACR_GDA(IX,IY,IZ) = TRACR_GDA(IX,IY,IZ)/FLOAT(IAVGGDCNTR)
            TCOLI_GDA(IX,IY,IZ) = TCOLI_GDA(IX,IY,IZ)/FLOAT(IAVGGDCNTR)
            FCOLI_GDA(IX,IY,IZ) = FCOLI_GDA(IX,IY,IZ)/FLOAT(IAVGGDCNTR)
            ENTERO_GDA(IX,IY,IZ)= ENTERO_GDA(IX,IY,IZ)/FLOAT(IAVGGDCNTR)
  360    CONTINUE
!
!------------------------------------------------------------
! Write out global dumped state variable into NetCDF files 
!------------------------------------------------------------
!
         IF(NETCDFOPT.EQ.1) THEN
          WRITE(OUT,'(A,i5.5,A,f15.5,A,i5.5,A,A)'),
     .              'WRITE: TOTAL IREC = ',IREC,' Time = ',time,
     .              ' NetCDF REC = ', ncIREC, ' in ',
     .              TRIM(ADJUSTL(OUTFILNA))
          status=nf90_open(TRIM(ADJUSTL(OUTFILNA)),nf90_write,ncID)
          CALL nccheck_status(status,OUTFILNA,RCANA)
          IF(INITB.EQ.0)
     .    CALL ncwrt_t1dvar(ncID,idpath_TIME    ,TIME)
          IF(INITB.GE.1)
     .    CALL ncwrt_t1dvar(ncID,idpath_TIME    ,
     .                      TIME-(FLOAT(IPRNTGSECS)/86400.)/2.)
          CALL ncwrt_g3dvar(ncID,idpath_SAL_GDA ,SAL_GDA)
          CALL ncwrt_g3dvar(ncID,idTRACR_GDA    ,TRACR_GDA)
          CALL ncwrt_g3dvar(ncID,idTCOLI_GDA    ,TCOLI_GDA)
          CALL ncwrt_g3dvar(ncID,idFCOLI_GDA    ,FCOLI_GDA)
          CALL ncwrt_g3dvar(ncID,idENTERO_GDA   ,ENTERO_GDA)
          status=nf90_close(ncID)
          CALL nccheck_status(status,OUTFILNA,RCANA)
          ncIREC=ncIREC+1
         ELSE
          IF(INITB.EQ.0) WRITE(10) TIME
          IF(INITB.GE.1) WRITE(10) TIME-(FLOAT(IPRNTGSECS)/86400.)/2.
          WRITE(11)   SAL_GDA
          WRITE(11)   TRACR_GDA
          WRITE(11)   TCOLI_GDA
          WRITE(11)   FCOLI_GDA
          WRITE(11)   ENTERO_GDA
         ENDIF
       ENDIF
         

       IF(IGDOPT.EQ.1)  THEN
c$doacross local(iz,iy,ix) , share(sal_gda,ecmsal_gda,dye1_gda,bac1_gda,
c$&        bac2_gda,bac3_gda)
         DO 380 IZ=1,NZ
          DO 380 IY=1,NY
           DO 380 IX=1,NX
            SAL_GDA(IX,IY,IZ) = 0.
            ECMSAL_GDA(IX,IY,IZ) = 0.
            TRACR_GDA(IX,IY,IZ) = 0.
            TCOLI_GDA(IX,IY,IZ) = 0.
            FCOLI_GDA(IX,IY,IZ) = 0.
            ENTERO_GDA(IX,IY,IZ) = 0.
  380    CONTINUE
       ENDIF

       IAVGGDCNTR = 0

       ENDIF

C        INITIAL CONDITION FILE
!
!------------------------------------------------------------
!  Optionally write restart file 
!------------------------------------------------------------
!
      IF(NETCDFOPT.EQ.1) THEN
        IF(ITIMESECS.GE.NXPRTR) THEN
        status=nf90_open(TRIM(ADJUSTL(RSTFILNA)),nf90_write,ncID_rst)
        CALL nccheck_status(status,RSTFILNA,RCANA)
        CALL ncwrt_rstvar(ncID_rst,ideutro_TIME,idTvar,TIME,CARAY)
        status=nf90_close(ncID_rst)
        CALL nccheck_status(status,RSTFILNA,RCANA)
        WRITE(OUT,'(A,i5.5,A,f15.5,A,i5.5,A,A)'),
     .           ' WRITE: RESTART path   IREC = ',IREC,
     .           ' Time = ',time,
     .           ' NetCDF REC = ', rstTflag,
     .           ' in ', TRIM(ADJUSTL(RSTFILNA))
        IF(rstTflag.EQ.0) rstTflag=1
        rstTflag=rstTflag+1
        IF(rstTflag.GT.2) rstTflag=1
        NXPRTR = NXPRTR + IPRNTRSECS
        ENDIF
      ELSE
        REWIND 15
        WRITE(15)  CARAY
      ENDIF

      RETURN

  900 WRITE(OUT,9990)
 9990 FORMAT(///5X,'INPUT ERROR WHILE READING IGDOPT,IDDOPT'//)
      CALL EXIT

C       ERROR OPENING KE FILE
  905 WRITE(OUT,9105)  KINFILNA(1)  
 9105 FORMAT(' Error encountered opening extinction coefficient file'/
     . 'file name =',A40)
      CALL EXIT
C       ERROR ON KE TIME BREAK READ
  910 WRITE(OUT,9110)  
 9110 FORMAT(' Error encountered in reading time breaks for the extincti
     .on coefficients')
      CALL EXIT
C       ERROR ON SEGMENT KE READ
  915 WRITE(OUT,9115)  IKE
 9115 FORMAT(' Error encountered in reading extinction coefficients for 
     .time interval ',I3)
      CALL EXIT
      END
