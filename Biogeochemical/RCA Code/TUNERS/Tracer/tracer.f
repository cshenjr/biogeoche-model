C
C***********************************************************************
C 
      SUBROUTINE   TUNER
C 
C***********************************************************************
C 
C                                                         TRACER Model
C                                                       ----------------
C            (Note: this current version is set up for 1 tracer, but may
C                 be readily expanded to 2 or 3 tracers by uncommenting
C                 the appropriate (CC) FORTRAN statements in the code
C                 and changing <NOSYS> in the RCACM PARAMETER statement)
C 
C***********************************************************************
C 
C     SYSTEMS                                                      UNITS
C     -------                                                      -----
C        1 - SALIN   - SALINITY                                      PPT
C        2 - TRACR1  - CONCENTRATION OF TRACER                      MG/L
C        3 - TRACR2  - CONCENTRATION OF TRACER                      MG/L
C        4 - TRACR3  - CONCENTRATION OF TRACER                      MG/L
C
C***********************************************************************
C 
C     CONSTANTS 
C     --------- 
C 
C   NO   NAME      DESCRIPTION                                     UNITS
C   --   --------  ----------------------------------------------  ----- 
C    1   DECAY1    DECAY RATE FOR TRACR 1                             /DAY
C    2   DECAY2    DECAY RATE FOR TRACR 2                             /DAY
C    3   DECAY3    DECAY RATE FOR TRACR 3                             /DAY
C
C***********************************************************************
C 
C     No 2-D PARAMETERS
C     --------------
C 
C***********************************************************************
C 
C     No 3-D PARAMETERS
C     --------------
C 
C***********************************************************************
C 
C     TIME-VARIABLE FUNCTIONS 
C     ----------------------- 
C           NONE
C 
C***********************************************************************
C  
      SAVE

      INCLUDE 'RCACM'
      CHARACTER   GDNAMES(NOSYS)*8,DDNAMES(5,NOSYS)*8

C        STATE-VARIABLES
      REAL
     .      SALIN(NX,NY,NZ)  , TRACR1(NX,NY,NZ)
CC   .    , TRACR2(NX,NY,NZ)
CC   .    , TRACR3(NX,NY,NZ)
      EQUIVALENCE
     .      (CARAY(1,1,1,1),SALIN(1,1,1))
     .    , (CARAY(1,1,1,2),TRACR1(1,1,1))
CC   .    , (CARAY(1,1,1,3),TRACR2(1,1,1))
CC   .    , (CARAY(1,1,1,4),TRACR3(1,1,1))
      REAL
     .      SALIN_DDA(NX,NY,NZ)       , SALIN_DMIN(NX,NY,NZ)  
     .    , SALIN_DMAX(NX,NY,NZ)  
     .    , TRACR_DDA(NX,NY,NZ,NOSYS) , TRACR_DMIN(NX,NY,NZ,NOSYS)
     .    , TRACR_DMAX(NX,NY,NZ,NOSYS) 
     .    , SALIN_GDA(NX,NY,NZ)
     .    , TRACR_GDA(NX,NY,NZ,NOSYS)
     .    , TRACR1_GDA(NX,NY,NZ)
CC   .    , TRACR2_GDA(NX,NY,NZ)
CC   .    , TRACR3_GDA(NX,NY,NZ)
      EQUIVALENCE
     .      (TRACR_GDA(1,1,1,1),TRACR1_GDA(1,1,1))
CC   .    , (TRACR_GDA(1,1,1,2),TRACR2_GDA(1,1,1))
CC   .    , (TRACR_GDA(1,1,1,3),TRACR3_GDA(1,1,1))

C        CONSTANTS
      EQUIVALENCE
     .   (CONST(1),DECAY1)  , (CONST(2),DECAY2) , (CONST(3),DECAY3) 

      INTEGER*2  SYSGDP(40)
C        INITIAL -SYSBY- SETTINGS
      DATA  SYSGDP/40*1/

C        PROVIDE INITIALIZATION, IF FIRST TIME THROUGH -FABLE-
      IF(INITB.EQ.1)   GO TO 50 
C        SET-UP AND WRITE INFORMATION NEEDED BY GDP
        GDNAMES( 1) = 'SALINITY'
        GDNAMES( 2) = 'TRACR1'
CC      GDNAMES( 3) = 'TRACR2'
CC      GDNAMES( 4) = 'TRACR3'
      DO ISYS=1,NOSYS
         SYSGDP(ISYS) = SYSBY(ISYS)
      ENDDO
      REWIND(10)
      WRITE(10)   NX,NY,NZ,NOSYS,NOSYS
      WRITE(10)   GDNAMES
      WRITE(10)   SYSGDP
      WRITE(10)   FSM
C        WRITE DDNAMES TO RCAF12
C         (NOTE: IF SYSTEM 2 AND/OR 3 ACTIVATED ADD APPROPRIATE DDNAMES)
      IF(IDDOPT.EQ.0) THEN
       DDNAMES(1,1)='SALINITY'
       DDNAMES(2,1)='   DUMMY'
       DDNAMES(3,1)='   DUMMY'
       DDNAMES(4,1)='   DUMMY'
       DDNAMES(5,1)='   DUMMY'
       DDNAMES(1,2)='  TRACR1'
       DDNAMES(2,2)='   DUMMY'
       DDNAMES(3,2)='   DUMMY'
       DDNAMES(4,2)='   DUMMY'
       DDNAMES(5,2)='   DUMMY'
CC     DDNAMES(1,3)='  TRACR2'
CC     DDNAMES(2,3)='   DUMMY'
CC     DDNAMES(3,3)='   DUMMY'
CC     DDNAMES(4,3)='   DUMMY'
CC     DDNAMES(5,3)='   DUMMY'
CC     DDNAMES(1,4)='  TRACR3'
CC     DDNAMES(2,4)='   DUMMY'
CC     DDNAMES(3,4)='   DUMMY'
CC     DDNAMES(4,4)='   DUMMY'
CC     DDNAMES(5,4)='   DUMMY'
      ELSE
       DDNAMES(1,1)=' AVE SAL'
       DDNAMES(2,1)=' MAX SAL'
       DDNAMES(3,1)=' MIN SAL'
       DDNAMES(4,1)='   DUMMY'
       DDNAMES(5,1)='   DUMMY'
       DDNAMES(1,2)='AVE TCR1'
       DDNAMES(2,2)='MAX TCR1'
       DDNAMES(3,2)='MIN TCR1'
       DDNAMES(4,2)='   DUMMY'
       DDNAMES(5,2)='   DUMMY'
CC     DDNAMES(1,3)='AVE TCR2'
CC     DDNAMES(2,3)='MAX TCR2'
CC     DDNAMES(3,3)='MIN TCR2'
CC     DDNAMES(4,3)='   DUMMY'
CC     DDNAMES(5,3)='   DUMMY'
CC     DDNAMES(1,4)='AVE TCR3'
CC     DDNAMES(2,4)='MAX TCR3'
CC     DDNAMES(3,4)='MIN TCR3'
CC     DDNAMES(4,4)='   DUMMY'
CC     DDNAMES(5,4)='   DUMMY'
      ENDIF
      WRITE(12)  DDNAMES

      DUMMY=0.0
C        SET INITIAL CONDITIONS FOR SALINITY USING HYDRO MODEL COMPUTATIONS
      DO 10 IZ=1,NZ
       DO 10 IY=1,NY
        DO 10 IX=1,NX
         CARAY(IX,IY,IZ,1) = HYDSAL(IX,IY,IZ)
   10 CONTINUE

C        INITIALIZE ARRAY FOR GLOBAL DUMP AVERAGING, IF REQUIRED
      IF(IGDOPT.EQ.1)  THEN
       DO 15 IZ=1,NZ
        DO 15 IY=1,NY
         DO 15 IX=1,NX
          SALIN_GDA(IX,IY,IZ) = 0.
   15  CONTINUE
       DO 25 ISYS=2,NOSYS
        DO 20 IZ=1,NZ
         DO 20 IY=1,NY
          DO 20 IX=1,NX
            TRACR_GDA(IX,IY,IZ,ISYS-1) = 0.
   20   CONTINUE
   25  CONTINUE
      ENDIF
      IAVGGDCNTR = 0
      IF(IDDOPT.EQ.1) THEN
       DO 30 IZ=1,NZ
        DO 30 IY=1,NY
         DO 30 IX=1,NX
            SALIN_DDA(IX,IY,IZ) = 0.
            SALIN_DMIN(IX,IY,IZ) = 1000.
            SALIN_DMAX(IX,IY,IZ) = -1000.
   30  CONTINUE
       DO 40 ISYS=2,NOSYS
        DO 35 IZ=1,NZ
         DO 35 IY=1,NY
          DO 35 IX=1,NX
            TRACR_DDA(IX,IY,IZ,ISYS-1) = 0.
            TRACR_DMIN(IX,IY,IZ,ISYS-1) = 1000.
            TRACR_DMAX(IX,IY,IZ,ISYS-1) = -1000.
   35   CONTINUE
   40  CONTINUE
      ENDIF
      IAVGDDCNTR = 0


   50 CONTINUE
C        LOOP FOR DETAILED DUMP AVERAGING, IF REQUIRED
      IF(IDDOPT.EQ.1)  THEN
       DO 120 IZ=1,NZ
        DO 120 IY=1,NY
         DO 120 IX=1,NX
          SALIN_DDA(IX,IY,IZ) = SALIN_DDA(IX,IY,IZ) + CARAY(IX,IY,IZ,1)
          SALIN_DMIN(IX,IY,IZ) =
     .          AMIN1(SALIN_DMIN(IX,IY,IZ),CARAY(IX,IY,IZ,1))
          SALIN_DMAX(IX,IY,IZ) =
     .          AMAX1(SALIN_DMAX(IX,IY,IZ),CARAY(IX,IY,IZ,1))
          DO 110 ISYS=2,NOSYS
            TRACR_DDA(IX,IY,IZ,ISYS-1) = TRACR_DDA(IX,IY,IZ,ISYS-1)
     .          + CARAY(IX,IY,IZ,ISYS)
            TRACR_DMIN(IX,IY,IZ,ISYS-1) =
     .          AMIN1(TRACR_DMIN(IX,IY,IZ,ISYS-1),CARAY(IX,IY,IZ,ISYS))
            TRACR_DMAX(IX,IY,IZ,ISYS-1) =
     .          AMAX1(TRACR_DMAX(IX,IY,IZ,ISYS-1),CARAY(IX,IY,IZ,ISYS))
  110   CONTINUE
  120  CONTINUE
       IAVGDDCNTR = IAVGDDCNTR + 1
      ENDIF

C        SYSTEM 1 - SALINITY
      DO 150 IZ=1,NZ
       DO 150 IY=1,NY
        DO 150 IX=1,NX
          CDARAY(IX,IY,IZ,1) = 0.0
  150 CONTINUE
      IF(IDISK.EQ.2 .OR. IDISK.EQ.3)  THEN
        DO 155 IDMP=1,NDMPS
         IX = IFDMPS(IDMP,1)
         IY = IFDMPS(IDMP,2)
         IZ = IFDMPS(IDMP,3)
         IF(IDDOPT.EQ.0)  THEN
           CALL RCAWBUF(1,CARAY(IX,IY,IZ,1),dummy
     .         ,dummy,dummy,dummy)
         ELSE
           CALL RCAWBUF(1,SALIN_DDA(IX,IY,IZ)/IAVGDDCNTR
     .         ,SALIN_DMAX(IX,IY,IZ),SALIN_DMIN(IX,IY,IZ),dummy,dummy)
         ENDIF
  155   CONTINUE
      ENDIF


C        SYSTEM 2(-4) - TRACR1(-3)
      DO 170 ISYS=2,NOSYS
       DO 165 IZ=1,NZ
        DO 165 IY=1,NY
         DO 165 IX=1,NX
          IF(FSM(IX,IY).NE.1.)  GO TO 165
          CDARAY(IX,IY,IZ,ISYS) = -CONST(ISYS-1)*CARAY(IX,IY,IZ,ISYS)
  165  CONTINUE
       IF(IDISK.EQ.2 .OR. IDISK.EQ.3)  THEN
        DO 160 IDMP=1,NDMPS
         IX = IFDMPS(IDMP,1)
         IY = IFDMPS(IDMP,2)
         IZ = IFDMPS(IDMP,3)
         IF(IDDOPT.EQ.0)  THEN
           CALL RCAWBUF(ISYS,CARAY(IX,IY,IZ,ISYS),dummy
     .         ,dummy,dummy,dummy)
         ELSE
           CALL RCAWBUF(ISYS,TRACR_DDA(IX,IY,IZ,ISYS-1)/IAVGDDCNTR
     .         ,TRACR_DMAX(IX,IY,IZ,ISYS-1),TRACR_DMIN(IX,IY,IZ,ISYS-1)
     .         ,dummy,dummy)
         ENDIF
  160   CONTINUE
       ENDIF
  170 CONTINUE

C        CONVERT TO MASS UNITS
c$doacross local(iz,iy,ix,isys)
       DO 200 ISYS=1,NOSYS
        DO 200 IZ=1,NZ
         DO 200 IY=1,NY
          DO 200 IX=1,NX
          CDARAY(IX,IY,IZ,ISYS) = BVOL(IX,IY,IZ)*CDARAY(IX,IY,IZ,ISYS)
  200  CONTINUE

      IF(IDISK.EQ.2 .OR. IDISK.EQ.3)  THEN
C        CLOSE BUFFER AND WRITE TO DISK
       CALL RCAWRIT

C        RE-INITIALIZE ARRAY FOR DETAILED DUMP AVERAGING, IF REQUIRED
       IF(IDDOPT.EQ.1)  THEN
        DO 220 IZ=1,NZ
         DO 220 IY=1,NY
          DO 220 IX=1,NX
           SALIN_DDA(IX,IY,IZ) = 0.
           SALIN_DMIN(IX,IY,IZ) = 1000.
           SALIN_DMAX(IX,IY,IZ) = -1000.
           DO 215 ISYS=2,NOSYS
            TRACR_DDA(IX,IY,IZ,ISYS-1) = 0.
            TRACR_DMIN(IX,IY,IZ,ISYS-1) = 1000.
            TRACR_DMAX(IX,IY,IZ,ISYS-1) = -1000.
  215      CONTINUE
  220   CONTINUE
        IAVGDDCNTR = 0
       ENDIF
      ENDIF

C        PERFORM GLOBAL DUMP AVERAGING, IF REQUIRED
      IF(IGDOPT.EQ.1)  THEN
       DO 245 IZ=1,NZ
        DO 245 IY=1,NY
         DO 245 IX=1,NX
          SALIN_GDA(IX,IY,IZ) = SALIN_GDA(IX,IY,IZ) + CARAY(IX,IY,IZ,1)
          DO 240 ISYS=2,NOSYS
           TRACR_GDA(IX,IY,IZ,ISYS-1) = TRACR_GDA(IX,IY,IZ,ISYS-1)
     .                             + CARAY(IX,IY,IZ,ISYS)
  240   CONTINUE
  245  CONTINUE
       IAVGGDCNTR = IAVGGDCNTR + 1
      ENDIF

C         CHECK IF TIME TO DUMP TO DISK
      IF(IDISK.EQ.0)   RETURN
C        GLOBAL DUMPS
      IF(IDISK.EQ.1 .OR. IDISK.EQ.3)  THEN
       IF(IGDOPT.EQ.0)  THEN
         WRITE(10)   TIME
         WRITE(11)   SALIN
         WRITE(11)   TRACR1
CC       WRITE(11)  TRACR2
CC       WRITE(11)  TRACR3
       ELSE
         DO 360 IZ=1,NZ
          DO 360 IY=1,NY
           DO 360 IX=1,NX
            SALIN_GDA(IX,IY,IZ) = SALIN_GDA(IX,IY,IZ)/
     .                                           FLOAT(IAVGGDCNTR)
            DO 355 ISYS=2,NOSYS
             TRACR_GDA(IX,IY,IZ,ISYS-1) = TRACR_GDA(IX,IY,IZ,ISYS-1)/
     .                                           FLOAT(IAVGGDCNTR)
  355     CONTINUE
  360    CONTINUE
         IF(INITB.EQ.0) WRITE(10)  TIME
         IF(INITB.GE.1) WRITE(10)  TIME - FLOAT(IPRNTGSECS)/86400./2.
         WRITE(11)   SALIN_GDA
         WRITE(11)   TRACR1_GDA
CC       WRITE(11)   TRACR2_GDA
CC       WRITE(11)   TRACR3_GDA
         DO 385 IZ=1,NZ
          DO 385 IY=1,NY
           DO 385 IX=1,NX
            SALIN_GDA(IX,IY,IZ) = 0.
            DO 380 ISYS=2,NOSYS
              TRACR_GDA(IX,IY,IZ,ISYS-1) = 0.
  380       CONTINUE
  385    CONTINUE
         IAVGGDCNTR = 0
       ENDIF
      ENDIF

C        INITIAL CONDITION FILE
      REWIND 15
      WRITE(15)  CARAY

      RETURN
      END
