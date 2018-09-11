C
C*******************************************************************************
C 
      SUBROUTINE   TUNER
C 
C*******************************************************************************
C 
C                           TENS - Perform 10's Test to Check Continuity
C                                  of the Hydrodynamic Model Files
C                           --------------------------------------------
C 
C*******************************************************************************
C 
C     SYSTEMS                                                           UNITS
C     -------                                                           -----
C        1 - CONC  - CONCENTRATION                                       MG/L
C
C*******************************************************************************
C 
C     CONSTANTS 
C     --------- 
C 
C   NO   NAME      DESCRIPTION                                          UNITS
C   --   --------  ---------------------------------------------------  ----- 
C    1   DIFFOPT   DIFFUSER INPUTS OPTION
C                  = 0, NO DIFFUSERS WERE USED IN THE HYDRODYNAMIC
C                       MODEL ... ALL INFLOWS ARE VIA BOUNDARIES
C                  = 1, DIFFUSERS WERE USED IN THE HYDRODYNAMIC MODEL
C                       THEREFORE, LOADS FROM THE DIFFUSERS WILL NEED
C                       TO BE COMPUTED AND ACCOUNTED FOR IN MASS
C                       BALANCE COMPUTATIONS
C    2   DIFCONC   CONCENTRATION TO BE USED IN MASS BALANCE/             MG/L
C                  CONTINUITY CHECK. USUALLY WOULD BE CHOSEN AS
C                  10 MG/L, HENCE, THE NAME "TENS" TEST, BUT THE USER
C                  MAY CHOSE ANY CONCENTRATION DESIRED.  HOWEVER,
C                  IT SHOULD BE CONSISTENT WITH CONCENTRATION
C                  CHOSEN FOR THE BOUNDARY CONDITIONS AND OTHER INPUTS
C    3   ERRLVL    LEVEL OF ERROR (BASED ON PERCENT) AT WHICH "ERROR"       %
C                  MESSAGES ARE OUTPUT FOR USER REVIEW (EXAMPLE, IF A
C                  USER SELECTS 1 PERCENT AND ASSIGNS DIFCONC=10 MG/L
C                  THEN ERROR MESSAGES WILL BE PRINTED IF THE
C                  CONCENTRATION IN ANY SEGMENT IS LESS THAN 9.9 MG/L
C                  OR GREATER THAN 10.1 MG/L)
C                  IF ERRLVL = 0.0 THEN NO "ERROR" PRINTOUT WILL BE
C                  GENERATED .... JUST STANDARD RCA CONCENTRATION
C                  PRINTOUTS
C    4   DIFFTRCK  OPTION TO TRACK (I.E., OUTPUT) DIFFUSER FLOWS
C                  = 0, DO NOT PRINT
C                  = 1, PRINT
C
C*******************************************************************************
C 
C     2-D PARAMETERS - NONE REQUIRED
C     --------------
C 
C*******************************************************************************
C 
C     3-D PARAMETERS - NONE REQUIRED
C     --------------
C 
C*******************************************************************************
C 
C     TIME-VARIABLE FUNCTIONS - NONE REQUIRED
C     ----------------------- 
C 
C*******************************************************************************
C  
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
      CHARACTER   GDNAMES(NOSYS)*8,DDNAMES(5,NOSYS)*8

C        STATE-VARIABLES
      REAL
     .      CONC(NX,NY,NZ)
      EQUIVALENCE
     .      (CARAY(1,1,1,1),CONC(1,1,1))
      REAL
     .      CONC_DDA(NX,NY,NZ) , CONC_DMIN(NX,NY,NZ)
     .    , CONC_DMAX(NX,NY,NZ), CONC_GDA(NX,NY,NZ)

C        CONSTANTS
      EQUIVALENCE
     .   (CONST(1),DIFFOPT)  , (CONST(2),DIFCONC) , (CONST(3),ERRLVL)
     . , (CONST(4),DIFFTRCK)

      REAL*8   TOTMASIC(NOSYS),TOTMASREG(NOSYS),TOTMASALL(NOSYS)

      INTEGER*2  SYSGDP(40)
C        INITIAL -SYSBY- SETTINGS
      DATA  SYSGDP/40*1/

      DATA  ISECSDIFPRTME/0/
!
!-----------------------------------------------------------------------
! Name of this rca file
!-----------------------------------------------------------------------
!
      RCANA='tens.f'
!
!----------------------------------------------------------------------
! Create NetCDF file for restart
!----------------------------------------------------------------------
!
      IF(NETCDFOPT.EQ.1 .AND. rstTflag.EQ.0) THEN
        RSTFILNA = TRIM(OUTFILPX)//'_rst.nc'
        CALL nccrt_rcasys(RSTFILNA)
      ENDIF

C        PROVIDE INITIALIZATION, IF FIRST TIME THROUGH -FABLE-
      IF(INITB.GE.1)   GO TO 50 

C        SET-UP AND WRITE INFORMATION NEEDED BY GDP
        GDNAMES( 1) = 'CONC'
        NGDMP=1
      DO ISYS=1,NOSYS
         SYSGDP(ISYS) = SYSBY(ISYS)
      ENDDO
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
C         (NOTE: IF SYSTEM 2 AND/OR 3 ACTIVATED ADD APPROPRIATE DDNAMES)
      IF(IDDOPT.EQ.0) THEN
       DDNAMES(1,1)='    CONC'
       DDNAMES(2,1)='   DUMMY'
       DDNAMES(3,1)='   DUMMY'
       DDNAMES(4,1)='   DUMMY'
       DDNAMES(5,1)='DOMNMASS'
      ELSE
       DDNAMES(1,1)='AVE CONC'
       DDNAMES(2,1)='MIN CONC'
       DDNAMES(3,1)='MAX CONC'
       DDNAMES(4,1)='   DUMMY'
       DDNAMES(5,1)='DOMNMASS'
      ENDIF
      WRITE(12)  DDNAMES

C        SET INITIAL CONDITIONS
         DO 10 IZ=1,NZ
          DO 10 IY=1,NY
           DO 10 IX=1,NX
            IF(FSM(IX,IY).EQ.1.)
     .            CARAY(IX,IY,IZ,1) = DIFCONC
   10    CONTINUE
      CALL RCA09
      IREC=IREC-1

C        CHECK IF BOUNDARY CONDITION INPUTS SPECIFIED BY USER
      IF(NOBC(1).GT.0) GO TO 25
C          IF NOT, CHECK TO SEE IF THERE ARE BOUNDARIES AND SET
C          THEM TO DIFCONC
      DO 20 IY=1,NY
       DO 20 IX=1,NX
        IF(FSM(IX,IY).EQ.-1. .OR. FSM(IX,IY).EQ.-2) THEN
          IF((NOBC(1)+NZ).GT.MXBC) GO TO 900
          DO 15 IZ=1,NZ
           IBC(1,NOBC(1)+IZ,1) = IX
           IBC(2,NOBC(1)+IZ,1) = IY
           IBC(3,NOBC(1)+IZ,1) = IZ
           BBC(NOBC(1)+IZ,1) = DIFCONC
           CARAY(IX,IY,IZ,1) = DIFCONC
   15     CONTINUE
          NOBC(1)=NOBC(1)+NZ
        ENDIF
   20 CONTINUE

C        INITIALIZE ARRAY FOR GLOBAL DUMP AVERAGING, IF REQUIRED
   25 IF(IGDOPT.EQ.1)  THEN
        DO 35 IZ=1,NZ
         DO 35 IY=1,NY
          DO 35 IX=1,NX
            CONC_GDA(IX,IY,IZ) = 0.
            CONC_DDA(IX,IY,IZ) = 0.
            CONC_DMIN(IX,IY,IZ) = 1000.
            CONC_DMAX(IX,IY,IZ) = -1000.
   35   CONTINUE
       IAVGGDCNTR = 0
       IAVGDDCNTR = 0
      ENDIF
      DUMMY=0.0

   50 CONTINUE

C        LOOP FOR DETAILED DUMP AVERAGING, IF REQUIRED
      IF(IDDOPT.EQ.1)  THEN
       DO 105 ISYS=1,NOSYS
        DO 100 IZ=1,NZ
         DO 100 IY=1,NY
          DO 100 IX=1,NX
            CONC_DDA(IX,IY,IZ) = CONC_DDA(IX,IY,IZ)
     .          + CARAY(IX,IY,IZ,1)
            CONC_DMIN(IX,IY,IZ) =
     .          AMIN1(CONC_DMIN(IX,IY,IZ),CARAY(IX,IY,IZ,1))
            CONC_DMAX(IX,IY,IZ) =
     .          AMAX1(CONC_DMAX(IX,IY,IZ),CARAY(IX,IY,IZ,1))
  100   CONTINUE
  105  CONTINUE
        IAVGDDCNTR = IAVGDDCNTR + 1
      ENDIF

C        SYSTEM 1 - CONC
      DOMNMASS=0.0
      DO 125 IZ=1,NZ
       DO 125 IY=1,NY
        DO 125 IX=1,NX
          CDARAY(IX,IY,IZ,1)=0.0
          DOMNMASS = DOMNMASS + BVOL(IX,IY,IZ)*CARAY(IX,IY,IZ,1)
  125 CONTINUE
      IF(IDIFFOPT.EQ.1) THEN
      IF(ISECSDIFPRTME.NE.NXHYDTSECS) WRITE(OUT,8011) TIME
 8011 FORMAT(10X,'Hydrodynamic Update at Day = ',F11.6)
        DO 135 ID=1,NODIFF
         IX=IDD(ID)
         IY=JDD(ID)
         IF(DIFFTRCK.EQ.1 .AND. QDIFF(ID).GT.0. .AND.
     .      ISECSDIFPRTME.NE.NXHYDTSECS) THEN
              WRITE(OUT,2000)  TIME,IX,IY,QDIFF(ID)
 2000       FORMAT(5X,'TIME = ',F8.3,' DAYS   IX,IY,QDIFF =',2I4,E12.3,
     .             ' M**3/DAY')
         ENDIF
         QSUM=0.
         DO 130 IZ=1,NZ
           CDARAY(IX,IY,IZ,1) = CDARAY(IX,IY,IZ,1) +
     .                   VDDIST(ID,IZ)*QDIFF(ID)*DIFCONC/BVOL(IX,IY,IZ)
           QSUM=QSUM+QX(IX,IY,IZ)-QX(IX+1,IY,IZ)
     .              +QY(IX,IY,IZ)-QY(IX,IY+1,IZ)
  130    CONTINUE
  135   CONTINUE
        IF(ISECSDIFPRTME.NE.NXHYDTSECS) ISECSDIFPRTME=NXHYDTSECS
      ENDIF
      IF(IDISK.EQ.2 .OR. IDISK.EQ.3)  THEN
        DO 145 IDMP=1,NDMPS
         IX = IFDMPS(IDMP,1)
         IY = IFDMPS(IDMP,2)
         IZ = IFDMPS(IDMP,3)
         IF(IDDOPT.EQ.0)  THEN
           CALL RCAWBUF(ISYS,CARAY(IX,IY,IZ,1),dummy
     .         ,dummy,dummy,DOMNMASS)
         ELSE
           CALL RCAWBUF(ISYS,CONC_DDA(IX,IY,IZ)/IAVGDDCNTR
     .         ,CONC_DMIN(IX,IY,IZ),CONC_DMAX(IX,IY,IZ)
     .         ,dummy,DOMNMASS)
         ENDIF
  145   CONTINUE
      ENDIF

C        CONVERT TO MASS UNITS
c$doacross local(iz,iy,ix,isys)
       DO 200 ISYS=1,NOSYS
        ERRVAL= 0.0
        DO 200 IZ=1,NZ
         DO 200 IY=1,NY
          DO 200 IX=1,NX
          CDARAY(IX,IY,IZ,ISYS) = BVOL(IX,IY,IZ)*CDARAY(IX,IY,IZ,ISYS)
          IF(FSM(IX,IY).EQ.1. .AND. ERRLVL.GT.0.0 .AND.
     .     ABS(CARAY(IX,IY,IZ,1)-DIFCONC)/DIFCONC.GT.ERRLVL/100.)
     .        WRITE(OUT,2200) TIME,IX,IY,IZ,CARAY(IX,IY,IZ,1)
 2200         FORMAT(10X,'ERROR AT TIME = ',F8.4,' DAYS IN SEGMENT ',I3,
     .        ',',I3,',',I2,' CONCENTRATION = ',F8.4,' MG/L')
          IF(FSM(IX,IY).EQ.1. .AND. ERRLVL.EQ.0.0 .AND.
     .     ABS(CARAY(IX,IY,IZ,1)-DIFCONC)/DIFCONC.GT.0.001*DIFCONC) THEN
              ERR=ABS(CARAY(IX,IY,IZ,1)-DIFCONC)/DIFCONC
              IF(ERR.GT.ERRVAL) THEN
                ERRVAL=ERR
                IERRX=IX
                IERRY=IY
                IERRZ=IZ
              ENDIF
          ENDIF
  200  CONTINUE
       IF(ERRLVL.EQ.0.0 .AND. ERRVAL.GT.0.0) WRITE(OUT,2300)  TIME,
     .        IERRX,IERRY,IERRZ,CARAY(IERRX,IERRY,IERRZ,1)
 2300  FORMAT(10X,'LARGEST ERROR AT TIME = ',F12.1,' DAYS IN SEGMENT '
     .        ,I3,',',I3,',',I2,' CONCENTRATION = ',F8.4,' MG/L')

      IF(IDISK.EQ.2 .OR. IDISK.EQ.3)  THEN
C        CLOSE BUFFER AND WRITE TO DISK
       CALL RCAWRIT

C        RE-INITIALIZE ARRAY FOR DETAILED DUMP AVERAGING, IF REQUIRED
       IF(IDDOPT.EQ.1)  THEN
        DO 215 ISYS=1,NOSYS
         DO 210 IZ=1,NZ
          DO 210 IY=1,NY
           DO 210 IX=1,NX
             CONC_DDA(IX,IY,IZ) = 0.
             CONC_DMIN(IX,IY,IZ) = 1000.
             CONC_DMAX(IX,IY,IZ) = -1000.
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
            CONC_GDA(IX,IY,IZ) = CONC_GDA(IX,IY,IZ)
     .                             + CARAY(IX,IY,IZ,1)
  220   CONTINUE
  225  CONTINUE
      ENDIF
      IAVGGDCNTR = IAVGGDCNTR + 1

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
         CALL ncwrt_t1dvar(ncID,id10s_TIME,TIME)
         CALL ncwrt_g2dvar(ncID,id10s_ETA ,ETA )
         CALL ncwrt_g3dvar(ncID,id10s_CONC,CONC)
         status=nf90_close(ncID)
         CALL nccheck_status(status,OUTFILNA,RCANA)
         ncIREC=ncIREC+1
        ELSE
         WRITE(10)   TIME
         WRITE(11)   CONC
        ENDIF
       ELSE
         DO 365 ISYS=1,NOSYS
          DO 360 IZ=1,NZ
           DO 360 IY=1,NY
            DO 360 IX=1,NX
             CONC_GDA(IX,IY,IZ) = CONC_GDA(IX,IY,IZ)/
     .                                           FLOAT(IAVGGDCNTR)
  360     CONTINUE
  365    CONTINUE
!
!------------------------------------------------------------
! Write out global dumped state variable into NetCDF files 
!------------------------------------------------------------
!
         IF(NETCDFOPT.EQ.1) THEN
          WRITE(*,'(A,i5.5,A,f15.5,A,i5.5,A,A)'),
     .            'WRITE: TOTAL IREC = ',IREC,' Time = ',time,
     .            ' NetCDF REC = ', ncIREC, ' in ',
     .            TRIM(ADJUSTL(OUTFILNA))
          status=nf90_open(TRIM(ADJUSTL(OUTFILNA)),nf90_write,ncID)
          CALL nccheck_status(status,OUTFILNA,RCANA)
          IF(INITB.EQ.0)
     .    CALL ncwrt_t1dvar(ncID,id10s_TIME    ,TIME)
          IF(INITB.GE.1)
     .    CALL ncwrt_t1dvar(ncID,id10s_TIME    ,
     .                      TIME-(FLOAT(IPRNTGSECS)/86400.)/2.)
          CALL ncwrt_g3dvar(ncID,id10s_CONC_GDA,CONC_GDA)
          status=nf90_close(ncID)
          CALL nccheck_status(status,OUTFILNA,RCANA)
          ncIREC=ncIREC+1
         ELSE
          IF(INITB.EQ.0) WRITE(10) TIME
          IF(INITB.GE.1) WRITE(10) TIME-(FLOAT(IPRNTGSECS)/86400.)/2.
          WRITE(11)   CONC_GDA
         ENDIF
       ENDIF
       IF(IGDOPT.EQ.1)  THEN
        DO 385 ISYS=1,NOSYS
         DO 380 IZ=1,NZ
          DO 380 IY=1,NY
           DO 380 IX=1,NX
            CONC_GDA(IX,IY,IZ) = 0.
  380    CONTINUE
  385   CONTINUE
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
     .           ' WRITE: RESTART tens   IREC = ',IREC,
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
 9990 FORMAT(///5X,'ERROR ... MXBC EXCEEDED IN SUBROUTINE TUNER (TENS VE
     .RSION)'/5X,'RCA TERMINATED')
      CALL EXIT
      END
