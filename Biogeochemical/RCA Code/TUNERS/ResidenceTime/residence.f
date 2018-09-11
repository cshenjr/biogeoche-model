C
C***********************************************************************
C 
      SUBROUTINE   TUNER
C 
C***********************************************************************
C 
C          RESIDENCE TIME MODEL - COMPUTE RESIDENCE TIME OF A SEGMENT(s)
C          -------------------------------------------------------------
C               (Note: current version is set up for 1 system, but may
C                  be readily expanded to 2 or 3 systems by uncommenting
C                  the appropriate FORTRAN statements in the code)
C 
C***********************************************************************
C 
C     SYSTEMS                                                      UNITS
C     -------                                                      -----
C        1 - DYE1  - CONCENTRATION OF DYE                           MG/L
C        2 - DYE2  - CONCENTRATION OF DYE                           MG/L
C        3 - DYE3  - CONCENTRATION OF DYE                           MG/L
C
C***********************************************************************
C 
C     CONSTANTS 
C     --------- 
C 
C   NO   NAME      DESCRIPTION                                     UNITS
C   --   --------  ----------------------------------------------  ----- 
C    1   DYEDIS1   DYE DISTRIBUTION OPTION
C                  = 0, ASSIGN INITIAL DYE CONCENTRATIONS UNIFORMLY
C                       THROUGHOUT THE WATER COLUMN (I.E., TOP TO 
C                       BOTTOM) USING PARAM2D(1)
C                       PARAM2D(2) WILL BE USED TO SPECIFY THE DOMAIN
C                       OVER WHICH THE RESIDENCE TIME WILL BE COMPUTED
C                  = 1, INITIAL DYE CONCENTRATIONS WILL BE SPECIFIED
C                       ON A SEGMENT BY SEGMENT BASIS USING PARAM3D(1)
C                       PARAM3D(2) WILL BE USED TO SPECIFY THE DOMAIN
C                       OVER WHICH THE RESIDENCE TIME WILL BE COMPUTED
C    2   TIMEDYE1  TIME TO RELEASE DYE 1
C    3   DECAY1    DECAY RATE FOR DYE 1                             /DAY
C    4   DYEDIS2   DYE DISTRIBUTION OPTION
C                  = 0, ASSIGN INITIAL DYE CONCENTRATIONS UNIFORMLY
C                       THROUGHOUT THE WATER COLUMN (I.E., TOP TO 
C                       BOTTOM) USING PARAM2D(3)
C                       PARAM2D(4) WILL BE USED TO SPECIFY THE DOMAIN
C                       OVER WHICH THE RESIDENCE TIME WILL BE COMPUTED
C                  = 1, INITIAL DYE CONCENTRATIONS WILL BE SPECIFIED
C                       ON A SEGMENT BY SEGMENT BASIS USING PARAM3D(3)
C                       PARAM3D(4) WILL BE USED TO SPECIFY THE DOMAIN
C                       OVER WHICH THE RESIDENCE TIME WILL BE COMPUTED
C    5   TIMEDYE2  TIME TO RELEASE DYE 2
C    6   DECAY2    DECAY RATE FOR DYE 2                             /DAY
C    7   DYEDIS3   DYE DISTRIBUTION OPTION
C                  = 0, ASSIGN INITIAL DYE CONCENTRATIONS UNIFORMLY
C                       THROUGHOUT THE WATER COLUMN (I.E., TOP TO 
C                       BOTTOM) USING PARAM2D(5)
C                       PARAM2D(6) WILL BE USED TO SPECIFY THE DOMAIN
C                       OVER WHICH THE RESIDENCE TIME WILL BE COMPUTED
C                  = 1, INITIAL DYE CONCENTRATIONS WILL BE SPECIFIED
C                       ON A SEGMENT BY SEGMENT BASIS USING PARAM3D(5)
C                       PARAM3D(6) WILL BE USED TO SPECIFY THE DOMAIN
C                       OVER WHICH THE RESIDENCE TIME WILL BE COMPUTED
C    8   TIMEDYE3  TIME TO RELEASE DYE 3
C    9   DECAY3    DECAY RATE FOR DYE 3                             /DAY
C
C***********************************************************************
C 
C     2-D PARAMETERS
C     --------------
C 
C       NO  NAME           DESCRIPTION                             UNITS
C      ---  --------  ------------------------------------------   ----- 
C        1  DYE1IC2   INITIAL CONDITIONS FOR DYE1                   MG/L
C        2  DOMAIN1   DOMAIN FOR DETERMINING RESIDENCE TIME 1
C        3  DYE2IC2   INITIAL CONDITIONS FOR DYE2                   MG/L
C        4  DOMAIN2   DOMAIN FOR DETERMINING RESIDENCE TIME 2
C        5  DYE3IC2   INITIAL CONDITIONS FOR DYE3                   MG/L
C        6  DOMAIN3   DOMAIN FOR DETERMINING RESIDENCE TIME 3
C 
C***********************************************************************
C 
C     3-D PARAMETERS
C     --------------
C 
C       NO  NAME           DESCRIPTION                             UNITS
C      ---  --------  ------------------------------------------   ----- 
C        1  DYE1IC3   INITIAL CONDITIONS FOR DYE1                   MG/L
C        2  DOMAIN1   DOMAIN FOR DETERMINING RESIDENCE TIME 1
C        3  DYE2IC3   INITIAL CONDITIONS FOR DYE2                   MG/L
C        4  DOMAIN2   DOMAIN FOR DETERMINING RESIDENCE TIME 2
C        5  DYE3IC3   INITIAL CONDITIONS FOR DYE3                   MG/L
C        6  DOMAIN3   DOMAIN FOR DETERMINING RESIDENCE TIME 3
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
     .      DYE1(NX,NY,NZ)
CC   .    , DYE2(NX,NY,NZ)
CC   .    , DYE3(NX,NY,NZ)
      EQUIVALENCE
     .      (CARAY(1,1,1,1),DYE1(1,1,1))
CC   .    , (CARAY(1,1,1,2),DYE2(1,1,1))
CC   .    , (CARAY(1,1,1,3),DYE3(1,1,1))
      REAL
     .      DYE_DDA(NX,NY,NZ,NOSYS)  
     .    , DYE_DMIN(NX,NY,NZ,NOSYS)
     .    , DYE_DMAX(NX,NY,NZ,NOSYS) 
     .    , DYE_GDA(NX,NY,NZ,NOSYS)
     .    , DYE1_GDA(NX,NY,NZ)
C    .    , DYE2_GDA(NX,NY,NZ)
C    .    , DYE3_GDA(NX,NY,NZ)
      EQUIVALENCE
     .      (DYE_GDA(1,1,1,1),DYE1_GDA(1,1,1))
C    .    , (DYE_GDA(1,1,1,2),DYE2_GDA(1,1,1))
C    .    , (DYE_GDA(1,1,1,3),DYE3_GDA(1,1,1))

C        CONSTANTS
      EQUIVALENCE
     .   (CONST(1),DYEDIS1)  , (CONST(2),TIMEDYE1) , (CONST(3),DECAY1) ,
     .   (CONST(4),DYEDIS2)  , (CONST(5),TIMEDYE2) , (CONST(6),DECAY2) ,
     .   (CONST(7),DYEDIS3)  , (CONST(8),TIMEDYE3) , (CONST(9),DECAY3) 

      REAL*8   TOTMASIC(NOSYS),TOTMASREG(NOSYS),TOTMASALL(NOSYS)

      INTEGER*2  SYSGDP(40)
C        INITIAL -SYSBY- SETTINGS
      DATA  SYSGDP/40*1/

C        PROVIDE INITIALIZATION, IF FIRST TIME THROUGH -FABLE-
      IF(INITB.EQ.1)   GO TO 50 
C        SET-UP AND WRITE INFORMATION NEEDED BY GDP
        GDNAMES( 1) = 'DYE1'
C       GDNAMES( 2) = 'DYE2'
C       GDNAMES( 3) = 'DYE3'
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
       DDNAMES(1,1)='    DYE1'
       DDNAMES(2,1)='   DUMMY'
       DDNAMES(3,1)='IC MASS1'
       DDNAMES(4,1)='REGNMAS1'
       DDNAMES(5,1)='DOMNMAS1'
       DDNAMES(1,2)='    DYE2'
       DDNAMES(2,2)='   DUMMY'
       DDNAMES(3,2)='IC MASS2'
       DDNAMES(4,2)='REGNMAS2'
       DDNAMES(5,2)='DOMNMAS2'
       DDNAMES(1,3)='    DYE3'
       DDNAMES(2,3)='   DUMMY'
       DDNAMES(3,3)='IC MASS3'
       DDNAMES(4,3)='REGNMAS3'
       DDNAMES(5,3)='DOMNMAS3'
      ELSE
       DDNAMES(1,1)='AVE DYE1'
       DDNAMES(2,1)='MAX DYE1'
       DDNAMES(3,1)='IC MASS1'
       DDNAMES(4,1)='REGNMAS1'
       DDNAMES(5,1)='DOMNMAS1'
       DDNAMES(1,2)='AVE DYE2'
       DDNAMES(2,2)='MAX DYE2'
       DDNAMES(3,2)='IC MASS2'
       DDNAMES(4,2)='REGNMAS2'
       DDNAMES(5,2)='DOMNMAS2'
       DDNAMES(1,3)='AVE DYE3'
       DDNAMES(2,3)='MAX DYE3'
       DDNAMES(3,3)='IC MASS3'
       DDNAMES(4,3)='REGNMAS3'
       DDNAMES(5,3)='DOMNMAS3'
      ENDIF
      WRITE(12)  DDNAMES

C        SET INITIAL CONDITIONS
      DO 30 ISYS=1,NOSYS
       IF(CONST(3*(ISYS-1)+1).EQ.0 .AND. CONST(3*(ISYS-1)+2).EQ.0.) THEN
         DO 10 IZ=1,NZ
          DO 10 IY=1,NY
           DO 10 IX=1,NX
            IF(FSM(IX,IY).EQ.1.)
     .            CARAY(IX,IY,IZ,ISYS) = PARAM2D(IX,IY,2*(ISYS-1)+1)
   10    CONTINUE
C        RESET CONSTANT SO DYE IS NOT PUT INTO SYSTEM AGAIN
         CONST(3*(ISYS-1)+2)=9.99E+9
       ELSEIF(CONST(3*(ISYS-1)+1).EQ.1 .AND. CONST(3*(ISYS-1)+2).EQ.0.)
     .  THEN
         DO 20 IZ=1,NZ
          DO 20 IY=1,NY
           DO 20 IX=1,NX
            IF(FSM(IX,IY).EQ.1.)
     .            CARAY(IX,IY,IZ,ISYS) = PARAM3D(IX,IY,IZ,2*(ISYS-1)+1)
   20    CONTINUE
         CONST(3*(ISYS-1)+2)=9.99E+9
       ENDIF
   30 CONTINUE
      CALL RCA09
      IREC=IREC-1

C        INITIALIZE ARRAY FOR GLOBAL DUMP AVERAGING, IF REQUIRED
      IF(IGDOPT.EQ.1)  THEN
       DO 40 ISYS=1,NOSYS
        DO 35 IZ=1,NZ
         DO 35 IY=1,NY
          DO 35 IX=1,NX
            DYE_GDA(IX,IY,IZ,ISYS) = 0.
            DYE_DDA(IX,IY,IZ,ISYS) = 0.
            DYE_DMIN(IX,IY,IZ,ISYS) = 1000.
            DYE_DMAX(IX,IY,IZ,ISYS) = -1000.
   35   CONTINUE
   40  CONTINUE
       IAVGGDCNTR = 0
       IAVGDDCNTR = 0
       IAVGPPCNTR = 0
      ENDIF
      DUMMY=0.0

   50 CONTINUE

C        CHECK IF TIME TO LOAD ANOTHER DYE
      DO 100 ISYS=1,NOSYS
       IF(CONST(3*(ISYS-1)+1).EQ.0 .AND. TIME.GT.CONST(3*(ISYS-1)+2))
     .  THEN
         WRITE(OUT,1200)  ISYS,TIME
 1200    FORMAT(//10X,'LOADING TRACER #',I2,' AT TIME',F8.2)
         DO 70 IZ=1,NZ
          DO 70 IY=1,NY
           DO 70 IX=1,NX
            CARAY(IX,IY,IZ,ISYS) = PARAM2D(IX,IY,2*(ISYS-1)+1)
   70    CONTINUE
         CONST(3*(ISYS-1)+2)=9.99E+9
         CALL RCA09
         IREC=IREC-1
       ELSEIF(CONST(3*(ISYS-1)+1).EQ.1.AND.TIME.GT.CONST(3*(ISYS-1)+2))
     .  THEN
         WRITE(OUT,1200)  ISYS,TIME
         DO 80 IZ=1,NZ
          DO 80 IY=1,NY
           DO 80 IX=1,NX
            CARAY(IX,IY,IZ,ISYS) = PARAM3D(IX,IY,IZ,2*(ISYS-1)+1)
   80    CONTINUE
         CONST(3*(ISYS-1)+2)=9.99E+9
         CALL RCA09
         IREC=IREC-1
       ENDIF
  100 CONTINUE

C        LOOP FOR DETAILED DUMP AVERAGING, IF REQUIRED
      IF(IDDOPT.EQ.1)  THEN
       DO 125 ISYS=1,NOSYS
        DO 120 IZ=1,NZ
         DO 120 IY=1,NY
          DO 120 IX=1,NX
            DYE_DDA(IX,IY,IZ,ISYS) = DYE_DDA(IX,IY,IZ,ISYS)
     .          + CARAY(IX,IY,IZ,ISYS)
            DYE_DMIN(IX,IY,IZ,ISYS) =
     .          AMIN1(DYE_DMIN(IX,IY,IZ,ISYS),CARAY(IX,IY,IZ,ISYS))
            DYE_DMAX(IX,IY,IZ,ISYS) =
     .          AMAX1(DYE_DMAX(IX,IY,IZ,ISYS),CARAY(IX,IY,IZ,ISYS))
  120   CONTINUE
  125  CONTINUE
        IAVGDDCNTR = IAVGDDCNTR + 1
      ENDIF

C        SYSTEM 1(-3) - DYE1(-3)
      DO 150 ISYS=1,NOSYS
       TOTMASIC(ISYS)=0.
       TOTMASREG(ISYS)=0.
       TOTMASALL(ISYS)=0.
       DO 135 IZ=1,NZ
        DO 135 IY=1,NY
         DO 135 IX=1,NX
          IF(FSM(IX,IY).NE.1.)  GO TO 135
          CDARAY(IX,IY,IZ,ISYS) = -CONST(3*ISYS)*CARAY(IX,IY,IZ,ISYS)
          IF(IDISK.NE.0) THEN
           IF(CONST(3*(ISYS-1)+1).EQ.0.) THEN
            IF(PARAM2D(IX,IY,2*(ISYS-1)+1).GT.0.)  TOTMASIC(ISYS) =
     .           TOTMASIC(ISYS) + BVOL(IX,IY,IZ)*CARAY(IX,IY,IZ,ISYS)
            IF(PARAM2D(IX,IY,2*ISYS).GT.0.)  TOTMASREG(ISYS) =
     .           TOTMASREG(ISYS) + BVOL(IX,IY,IZ)*CARAY(IX,IY,IZ,ISYS)
           ELSE
            IF(PARAM3D(IX,IY,IZ,2*(ISYS-1)+1).GT.0.)  TOTMASIC(ISYS) =
     .           TOTMASIC(ISYS) + BVOL(IX,IY,IZ)*CARAY(IX,IY,IZ,ISYS)
            IF(PARAM3D(IX,IY,IZ,2*ISYS).GT.0.)  TOTMASREG(ISYS) =
     .           TOTMASREG(ISYS) + BVOL(IX,IY,IZ)*CARAY(IX,IY,IZ,ISYS)
           ENDIF
           TOTMASALL(ISYS) =
     .           TOTMASALL(ISYS) + BVOL(IX,IY,IZ)*CARAY(IX,IY,IZ,ISYS)
          ENDIF
  135  CONTINUE
       IF(IDISK.EQ.2 .OR. IDISK.EQ.3)  THEN
        DO 145 IDMP=1,NDMPS
         IX = IFDMPS(IDMP,1)
         IY = IFDMPS(IDMP,2)
         IZ = IFDMPS(IDMP,3)
         IF(IDDOPT.EQ.0)  THEN
           CALL RCAWBUF(ISYS,CARAY(IX,IY,IZ,ISYS),dummy
     .         ,SNGL(TOTMASIC(ISYS)),SNGL(TOTMASREG(ISYS))
     .         ,SNGL(TOTMASALL(ISYS)))
         ELSE
           CALL RCAWBUF(ISYS,DYE_DDA(IX,IY,IZ,ISYS)/IAVGDDCNTR
     .         ,DYE_DMAX(IX,IY,IZ,ISYS)
     .         ,SNGL(TOTMASIC(ISYS)),SNGL(TOTMASREG(ISYS))
     .         ,SNGL(TOTMASALL(ISYS)))
         ENDIF
  145   CONTINUE
       ENDIF
  150 CONTINUE

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

       DO ISYS=1,NOSYS
        WRITE(OUT,2500)   TOTMASIC(ISYS)/1000.,TOTMASREG(ISYS)/1000.
     .                   ,TOTMASALL(ISYS)/1000.
 2500   FORMAT(
     .     /10X,'INITIAL CONDITION REGIONAL TOTAL MASS =',E13.5,' KG'
     .     /10X,'REGIONAL TOTAL MASS =',E13.5,' KG'/
     .     /10X,'DOMAIN TOTAL MASS =',E13.5,' KG'/)
       ENDDO
C        RE-INITIALIZE ARRAY FOR DETAILED DUMP AVERAGING, IF REQUIRED
       IF(IDDOPT.EQ.1)  THEN
        DO 215 ISYS=1,NOSYS
         DO 210 IZ=1,NZ
          DO 210 IY=1,NY
           DO 210 IX=1,NX
             DYE_DDA(IX,IY,IZ,ISYS) = 0.
             DYE_DMIN(IX,IY,IZ,ISYS) = 1000.
             DYE_DMAX(IX,IY,IZ,ISYS) = -1000.
  210    CONTINUE
  215   CONTINUE
        IAVGDDCNTR = 0
       ENDIF
      ENDIF
C        PERFORM GLOBAL DUMP AVERAGING, IF REQUIRED
      IF(IGDOPT.EQ.1)  THEN
       DO 225 ISYS=1,NOSYS
        DO 220 IZ=1,NZ
         DO 220 IY=1,NY
          DO 220 IX=1,NX
            DYE_GDA(IX,IY,IZ,ISYS) = DYE_GDA(IX,IY,IZ,ISYS)
     .                             + CARAY(IX,IY,IZ,ISYS)
  220   CONTINUE
  225  CONTINUE
      ENDIF
      IAVGGDCNTR = IAVGGDCNTR + 1

C         CHECK IF TIME TO DUMP TO DISK
      IF(IDISK.EQ.0)   RETURN
C        GLOBAL DUMPS
      IF(IDISK.EQ.1 .OR. IDISK.EQ.3)  THEN
       IF(IGDOPT.EQ.0)  THEN
         WRITE(10)   TIME
         WRITE(11)   DYE1
CC       WRITE(11)  DYE2
CC       WRITE(11)  DYE3
       ELSE
         DO 365 ISYS=1,NOSYS
          DO 360 IZ=1,NZ
           DO 360 IY=1,NY
            DO 360 IX=1,NX
             DYE_GDA(IX,IY,IZ,ISYS) = DYE_GDA(IX,IY,IZ,ISYS)/
     .                                           FLOAT(IAVGGDCNTR)
  360     CONTINUE
  365    CONTINUE
         IF(INITB.EQ.0)  WRITE(10)   TIME
         IF(INITB.GE.1)  WRITE(10)   TIME - FLOAT(IPRNTGSECS)/86400./2.
         WRITE(11)   DYE1_GDA
CC       WRITE(11)   DYE2_GDA
CC       WRITE(11)   DYE3_GDA
       ENDIF
       IF(IGDOPT.EQ.1)  THEN
        DO 385 ISYS=1,NOSYS
         DO 380 IZ=1,NZ
          DO 380 IY=1,NY
           DO 380 IX=1,NX
            DYE_GDA(IX,IY,IZ,ISYS) = 0.
  380    CONTINUE
  385   CONTINUE
       ENDIF
       IAVGGDCNTR = 0
       IAVGPPCNTR = 0
      ENDIF

C        INITIAL CONDITION FILE
      REWIND 15
      WRITE(15)  CARAY

      RETURN

  900 WRITE(OUT,9990)
 9990 FORMAT(///5X,'INPUT ERROR WHILE READING IGDOPT,IDDOPT'//)
      CALL EXIT
      END
