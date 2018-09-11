      SUBROUTINE RCA03
C 
C        RCA03 READS THE MODEL GEOMETRY AND INITIAL TRANSPORT FIELDS
C              AS PROVIDED BY THE HYDRODYNAMIC MODEL
C
!=======================================================================
! The file was revised to read ROMS (Region Ocean Modeling System)     !
! grid and initial transport fields in NetCDF format                   !
!                                 YUN LI, UMCES/HPL Oct-13-2008        !
!=======================================================================
!
      USE netcdf
! 
      SAVE
      INCLUDE   'RCACM' 
      INCLUDE   'NetCDFCM'

      INTEGER    IMINX(20),IMINY(20),IMINZ(20),TOTLWATR
      REAL
     .     TPS(NX,NY),MINDT(20),SCRATCH(NX*NY*NZ)
     .   , TIMESTEP(NX,NY,NZ),TIMESTEPS(NX*NY*NZ)
     .   , QBAL(NX,NY),TVDER(NX,NY)
      CHARACTER*4  UNITS
      REAL DUMMY(NX,NY,NZ)

C
      EQUIVALENCE
     .     (TPS(1,1),EZ(1,1,1)) , (SCRATCH(1),TIMESTEP(1,1,1))
     .   , (TIMESTEP(1,1,1),SCRATCH_TRAN(1,1,1,1))

      DATA   MINDT/20*100./ , TOTLWATR/0/
!
!-----------------------------------------------------------------------
! Name of this rca file
!-----------------------------------------------------------------------
!
      RCANA='rca03.f'

      WRITE(OUT,8000) 
 8000 FORMAT(//1X,119('*')//)
  
      READ(IN,1000)  COMMENT
 1000 FORMAT(A)
      READ(IN,1100)  ICOLLOPT,HYDCYCOPT,LNDWTROPT,IDIFFOPT,IECOMVER
 1100 FORMAT(5I10)
      IF(ICOLLOPT.EQ.1)  WRITE(OUT,1102)
 1102 FORMAT(33X,'RCA IS BEING RUN USING A COLLAPSED VERSION OF THE HYDR
     .ODYNAMIC MODEL'/)
      IF(HYDCYCOPT.EQ.1)  WRITE(OUT,1103)
 1103 FORMAT(37X,'HYDRODYNAMIC FILE WILL BE CYCLED IF END-OF-FILE ENCOUN
     .TERED'/)
      IF(LNDWTROPT.EQ.1)  WRITE(OUT,1104)
 1104 FORMAT(35X,'WATER SEGMENTS ONLY READ OPTION SELECTED FOR HYDRODYNA
     .MIC INPUT')
      IF(IDIFFOPT.EQ.1)  WRITE(OUT,1105)
 1105 FORMAT(35X,'USER IS SPECIFYING DIFFUSER INPUT FILES')
      QTOL=0.01
      IF(IDIFFOPT.EQ.1) QTOL=0.001
      IF(IECOMVER.EQ.1)  WRITE(OUT,1106)
 1106 FORMAT(35X,'USER SPECIFIED THAT ECOM RESULTS WILL BE READ USING OL
     .D FILE FORMAT STRUCTURE')
      READ(IN,1000)  COMMENT
      READ(IN,1107)  IHYDDT,UNITS
 1107 FORMAT(I10,6X,A4)
      IF(UNITS.NE.'SECS'.AND.UNITS.NE.'secs'.AND.
     .   UNITS.NE.'MINS'.AND.UNITS.NE.'mins'.AND.
     .   UNITS.NE.'HRS '.AND.UNITS.NE.'hrs '.AND.
     .   UNITS.NE.' HRS'.AND.UNITS.NE.' hrs') THEN
       WRITE(OUT,9001)  UNITS
 9001  FORMAT(/5X,'UNITS CHOSEN FOR HYDRODYNAMIC UPDATE INTERVAL ',
     .     A4,' ARE NOT VALID'/5X,'RCA TERMINATED')
       CALL EXIT
      ENDIF
      IF(UNITS.EQ.'SECS'.OR.UNITS.EQ.'secs')  SCALHYD=1.
      IF(UNITS.EQ.'MINS'.OR.UNITS.EQ.'mins')  SCALHYD=60.
      IF(UNITS.EQ.'HRS '.OR.UNITS.EQ.'hrs ')  SCALHYD=3600.
      IF(UNITS.EQ.' HRS'.OR.UNITS.EQ.' hrs')  SCALHYD=3600.
C        POSSIBLE FUTURE CODE IF HYDRODYNAMIC RECORDS NOT SAVED AT HOURLY
C          OR CONSTANT INTERVALS
C     IF(HYDDT.NE.0.0)  THEN
C       DO I=1,NHYD
C        HYDBRK(I) = FLOAT(I)*HYDDT*SCALHYD
C       ENDDO
C     ELSE
C       READ(IN,1100)  NOHYD
C       IF(NOHYD.GT.NHYD)  GO TO 960
C       READ(IN,1108,ERR=950)  (HYDBRK(I),I=1,NOHYD)
C1108   FORMAT(8F10.0)
C     ENDIF
C     DO I=1,NOHYD
C      HYDBRK(I)=SCALHYD*HYDBRK(I)
C     ENDDO
      IHYDDTSECS = IHYDDT*SCALHYD
      HYDDT = FLOAT(IHYDDTSECS)/86400.0d0
      NXHYDTSECS = 0
      NXHYDT = 0.0

C        CHECK TO MAKE SURE THAT INTEGRATION STEPSIZE(S) ARE EXACT
C        MULTIPLES OF THE HYDRODYNAMIC INTERVAL
      IF(INTGRTYP.EQ.1 .OR. INTGRTYP.EQ.4 .OR. INTGRTYP.EQ.5) THEN
       DO I=1,NSTEP
        IREMAIN=IMOD(IHYDDTSECS,ISTEP(I))
        IF(IREMAIN.NE.0) THEN
         WRITE(OUT,1150)  I,ISTEP(I),IHYDDTSECS
 1150    FORMAT(//5X,'INPUT ERROR ... INTEGRATION STEP NUMBER',I3,
     .    ' (DT =',I6,' SECONDS) IS NOT AN EXACT MULTIPLE OF THE',1X,
     .    'HYDRODYNAMIC UPDATE INTERVAL OF ',I8,' SECONDS'/5X,
     .    'RCA TERMINATED')
         CALL EXIT
        ENDIF
       ENDDO
      ELSE
       IREMAINS=IMOD(IHYDDTSECS,IDTSPLITSECS)
       IREMAINF=IMOD(IHYDDTSECS,IDTFULLSECS)
       IF(IREMAINS.NE.0 .OR. IREMAINF.NE.0) THEN
         IF(IREMAINS.NE.0)  WRITE(OUT,1151) IDTSPLITSECS,IHYDDTSECS
 1151    FORMAT(//5X,'INPUT ERROR ... SPLIT INTEGRATION TIMESTEP',1X,
     .    'IDTSPLIT = ',I6,' SECONDS IS NOT AN EXACT MULTIPLE OF THE'
     .    1X,'HYDRODYNAMIC UPDATE INTERVAL OF ',I8,' SECONDS'/5X,
     .    'RCA TERMINATED')
         IF(IREMAINF.NE.0)  WRITE(OUT,1152) IDTFULLSECS,IHYDDTSECS
 1152    FORMAT(//5X,'INPUT ERROR ... FULL INTEGRATION TIMESTEP',1X,
     .    'IDTFULL = ',I6,' SECONDS IS NOT AN EXACT MULTIPLE OF THE'
     .    1X,'HYDRODYNAMIC UPDATE INTERVAL OF ',I6,' SECONDS'/5X,
     .    'RCA TERMINATED')
         CALL EXIT
       ENDIF
      ENDIF

C     Advective and Dispersive Scale Factors
      READ(IN,1000)   COMMENT
      READ(IN,1200)  SCALRX,SCALRY,SCALRZ
 1200 FORMAT(3F10.0)

      WRITE(OUT,1210)  SCALRX,SCALRY,SCALRZ
 1210 FORMAT(//
     .  40X,'X-DISPERSION SCALE FACTOR =',E11.3/
     .  40X,'Y-DISPERSION SCALE FACTOR =',E11.3/
     .  40X,'Z-DISPERSION SCALE FACTOR =',E11.3)

!
!----------------------------------------------------------------------
! Optionally read horizontal and vertical grid infomation
!!      IBNRHYDOPT=1, RCA reads ECOM binary file 
!!      IBNRHYDOPT=2, RCA reads ROMS NetCDF file 
! similar hereafter.
! NOTE: screen display format was revised, I2->I1.
!----------------------------------------------------------------------
!
!      READ(IN,1000)  COMMENT
!      READ(IN,'(I10)')  IBNRHYDOPT
!      IF(NETCDFOPT.EQ.1) IBNRHYDOPT=2

      IF(NETCDFOPT.EQ.1) THEN
        READ(IN,1000)  COMMENT
        READ(IN,'(A)')  VARNA    ! RCA netcdf variables
        READ(IN,'(A)')  GRDNA    ! ROMS grid file
        READ(IN,'(A)')  RVRNA    ! ROMS river file
        READ(IN,'(A)')  HISNA    ! ROMS history file
        WRITE(OUT,'(//)')
        CALL ncset_ncparams
        CALL ncset_horizgrid     ! load H,DX,DY,FSM
        CALL ncset_vertigrid     ! load DZ,DZZ
        INCLUDE 'ncTid_list.h'
      ELSE
        OPEN(30,FILE='gcm_geom',FORM='UNFORMATTED')
        READ(30)  DZ,DZZ
C        WITH RELEASE OF VERSION 3.0, RCA NOWS EXPECTS TO READ INFO
C        ASSOCIATED WITH DIFFUSER DISCHARGES INCLUDED IN ECOMSED/ECOMSI
        READ(30,ERR=5)  H,DX,DY,FSM
        GO TO 7
    5   READ(30)  H,DX,DY,FSM
    7   IF(ICOLLOPT.EQ.1) READ(30)  XAZ
        CLOSE(30)
      ENDIF

C        OUTPUT LAND MASK ARRAY
      WRITE(OUT,1120) (I,I=1,9)
     .               ,(I,I=1,9),0,(I,I=1,9),0,(I,I=1,9),0
     .               ,(I,I=1,9),0,(I,I=1,9),0,(I,I=1,9),0
     .               ,(I,I=1,9),0,(I,I=1,9),0,(I,I=1,9),0
 1120 FORMAT(//45X,'COMPUTATIONAL GRID (FSM LAND MASK ARRAY)'/
     .          2X,'IY  IX-->'/5X,9(9X,I1)/5X,90I1)
      DO IY=NY,1,-1
      WRITE(OUT,1122)  IY,(IFIX(FSM(IX,IY)),IX=1,NX)
 1122  FORMAT(I4,1X,<NX>I1)
      ENDDO
      WRITE(OUT,1122)

C        GENERATE -ZZ- ARRAY
      ZZ(1) = DZ(1)/2.
      DO  10 IZ=2,NZ
        ZZ(IZ) = ZZ(IZ-1) + DZZ(IZ-1)
   10 CONTINUE

!
!-------------------------------------------------------
! Optionally write information file 
!-------------------------------------------------------
!
      IF(NETCDFOPT.EQ.1) THEN
      ELSE
C        WRITE LAND-MASK ARRAY FOR USE IN GDP
      WRITE(10)   FSM
      ENDIF

C        FIND MIN-MAX WATER SEGMENTS FOR EACH ROW

      DO 20 IY=1,NY
       IXS(IY) = 0
       DO 15 IX=1,NX
        IF(FSM(IX,IY).EQ.1.)  THEN
          TOTLWATR = TOTLWATR+1
          IF(IXS(IY).EQ.0)  IXS(IY) = IX 
          IXE(IY) = IX
        ENDIF
   15  CONTINUE
       IXS(IY) = IXS(IY)-1
       IXE(IY) = IXE(IY)+1
       IF(IXS(IY).LT.1)  IXS(IY)=1
       IF(IXE(IY).GT.NX)  IXE(IY)=NX
   20 CONTINUE
      TOTLWATR = NZ*TOTLWATR
         
C        IF WATER SEGMENTS ONLY READ OPTION SELECTED THEN GET
C        LOCATIONS OF THE 'WET' GRID
      IF(LNDWTROPT.EQ.1)  THEN
        OPEN(30,FILE='wet_grid',FORM='UNFORMATTED')
        READ(30)  IWTRCNT
        DO 25 I=1,IWTRCNT
          READ(30)  IWTRNDX(I),JWTRNDX(I)
   25   CONTINUE
        CLOSE(30)
      ENDIF


      READ(IN,1000) COMMENT
      READ(IN,1100) NOHYDFILNA
      IF(NOHYDFILNA.GT.MXHYDFILES)  THEN
        WRITE(OUT,1260)   NOHYDFILNA,MXHYDFILES
 1260   FORMAT(//
     .   10X,'ERROR ... USER REQUESTED',I5,' HYDRODYNAMIC INPUT FILES'/
     .   10X,'THIS IS GREATER THAN MAXIMUM NUMBER PERMITTED (',I3,')'/
     .   10X,'CORRECT INPUT OR INCREASE RCACM VARIABLE -MXHYDFILES-')
      ENDIF
      DO IHYDFILE=1,NOHYDFILNA
        READ(IN,1300)  HYDFILNA(IHYDFILE),DIFFILNA(IHYDFILE)
 1300   FORMAT(A52,A40) !altered to read 2001 transport (Testa 9/19/2012)
      ENDDO
      IHYDFILE=1
      WRITE(OUT,1310)  HYDFILNA(1),DIFFILNA(1)
 1310 FORMAT(//20X,'OPENING HYDRODYNAMIC TRANPORT FILE = ',A52/
     .         20X,'AND OPENING DIFFUSER DISCHARGE FILE = ',A52/)
!
!----------------------------------------------------------------------
! Optionally open hydrodynamic file
!----------------------------------------------------------------------
!
      IF(NETCDFOPT.EQ.1) THEN
        status=nf90_open(TRIM(ADJUSTL(HYDFILNA(IHYDFILE))),
     .                   nf90_nowrite,ncID_his)
        CALL nccheck_status(status,HYDFILNA(IHYDFILE),RCANA)
      ELSE
        OPEN(30,FILE=HYDFILNA(1),FORM='UNFORMATTED',
     .       STATUS='OLD',ERR=900)
        IF(IDIFFOPT.EQ.1)  THEN
          OPEN(31,FILE=DIFFILNA(1),FORM='UNFORMATTED',STATUS='OLD',
     .         ERR=910)
          READ(31)   NODIFF
          IF(NODIFF.GT.MXDIFF)  GO TO 980
          BACKSPACE(31)
          READ(31)  NODIFF,(IDD(I),I=1,NODIFF),(JDD(I),I=1,NODIFF),
     .     ((VDDIST(I,J),J=1,NZ),I=1,NODIFF)
C       IF(LIST(1).EQ.1) THEN
          WRITE(OUT,1140)
 1140     FORMAT(//30X,'THE FOLLOWING DIFFUSERS WERE UTILIZED IN THE HYD
     .RODYNAMIC MODEL'/10X,'IX   IY   VERTICAL FLOW DISTRIBUTION'/)
          DO I=1,NODIFF
            WRITE(OUT,1142)  IDD(I),JDD(I),(VDDIST(I,J),J=1,NZ)
 1142       FORMAT(9X,I3,I5,3X,10F5.0)
            DO J=1,NZ
              VDDIST(I,J)=VDDIST(I,J)/100.
            ENDDO
          ENDDO
C       ENDIF
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
! Model startup if ROMS NetCDF is used
!-----------------------------------------------------------------------
!
      IF(NETCDFOPT.EQ.1) THEN
        ncIREC=1
        hydTflag=1
        CALL ncset_oceantime(ncID_his,hydTflag,TMID,not_enough_length)
        NXHYDTSECS=int(TMID)
        NXHYDT=TMID/86400.0d0
        HYDTZERO=NXHYDT
      ENDIF

      IF(IDIAGDT.EQ.1 .OR. INPCHCK.EQ.1)   TIME=0.


      ENTRY   RCA03A
!
!-----------------------------------------------------------------------
! Name of this rca file
!-----------------------------------------------------------------------
!
      RCANA='rca03.f'
!
!----------------------------------------------------------------------
! Get model time 'TMID'
!----------------------------------------------------------------------
!
  30  IF(NETCDFOPT.EQ.1) THEN
        CALL ncset_oceantime(ncID_his,hydTflag,TMID,not_enough_length)
        IF(not_enough_length) GO TO 920    ! no record found
        TMID=TMID/86400.0d0                ! change the unit, sec -> day
      ELSE
        READ(30,END=920)  TMID
      ENDIF
!
!----------------------------------------------------------------------
! Create NetCDF file for rca TUNER output
!----------------------------------------------------------------------
!
      IF(NETCDFOPT.EQ.1 .AND. ncIREC.EQ.1) THEN
        nlen_HYDFILNA=LEN(ADJUSTL(TRIM(HYDFILNA(IHYDFILE))))
        OUTFILNA = TRIM(OUTFILPX)
     .           //HYDFILNA(IHYDFILE)(nlen_HYDFILNA-7:nlen_HYDFILNA)
        CALL nccrt_rcasys(OUTFILNA)
      ENDIF

      IF(INITB.EQ.0)  TMID0=TMID
      IF(IDIAGDT.EQ.1 .OR. INPCHCK.EQ.1)   THEN
        INITB=1
        TIME = TMID - TMID0
      ENDIF
      IF(IDIAGDT.EQ.1 .OR. INPCHCK.EQ.1) WRITE(OUT,8011)  TIME, TMID,
     .                  hydTflag
 8011 FORMAT(10X,'Hydrodynamic Update at Day = ',F11.6,' (Hydro Time =)'
     .         ,F11.6, 1X,I3)
      IF(IDIFFOPT.EQ.1)  THEN
       READ(31,END=930)  TMIDD,(QDIFF(I),I=1,NODIFF)
       IF(TMID.NE.TMIDD)  GO TO 990
      ENDIF

C SET FLAG FOR READING SETTLING RATE AND RESUSPENSION FLUX IN "TUNER"
C     ISEDTRAN (=0,DO NOT READ),(=1,READ)
      ISEDTRAN=1
!
!----------------------------------------------------------------------
! OPTIONALLY READ FULL LAND/WATER HYDRODYNAMIC DATA FROM
!          TMID QX/QY/QZ, EX/EY/EZ, ETA, DETA
!----------------------------------------------------------------------
!
      IF(NETCDFOPT.EQ.1) THEN
         CALL ncset_transport(ncID_his)
         hydTflag=hydTflag+1
      ELSE
        IF(LNDWTROPT.EQ.0)  THEN
C        READ FULL LAND/WATER GRID
          READ(30)  QX
          READ(30)  QY
          IF(NZ.GT.1) READ(30)  QZ
          READ(30)  EX
          READ(30)  EY
          IF(NZ.GT.1) READ(30)  EZ
          IF(IECOMVER.EQ.0 .AND. NZ.GT.1)  READ(30) DUMMY
          READ(30)  ETA
          READ(30)  DETA
          IF(NZ.GT.1) READ(30)  HYDSAL
          IF(NZ.GT.1) READ(30)  HYDTEMP
        ELSE
C       Read water only grid and put into full grid
          CALL LNDWTR(QX,SCRATCH,IWTRCNT,NZ)
          CALL LNDWTR(QY,SCRATCH,IWTRCNT,NZ)
          IF (NZ.GT.1)  CALL LNDWTR(QZ,SCRATCH,IWTRCNT,NZ)
C       Check if uncollapsed/collapsed read
          IF(ICOLLOPT.EQ.0)  THEN
            CALL LNDWTR(EX,SCRATCH,IWTRCNT,NZ)
            CALL LNDWTR(EY,SCRATCH,IWTRCNT,NZ)
            IF(NZ.GT.1)  CALL LNDWTR(EZ,SCRATCH,IWTRCNT,NZ)
            IF(IECOMVER.EQ.0 .AND. NZ.GT.1)
     .        CALL LNDWTR(DUMMY,SCRATCH,IWTRCNT,NZ)
            CALL LNDWTR(ETA,SCRATCH,IWTRCNT,1)
            CALL LNDWTR(DETA,SCRATCH,IWTRCNT,1)
          ELSE
            CALL LNDWTR(RX,SCRATCH,IWTRCNT,NZ)
            CALL LNDWTR(RY,SCRATCH,IWTRCNT,NZ)
            IF(NZ.GT.1)  CALL LNDWTR(RZ,SCRATCH,IWTRCNT,NZ)
            CALL LNDWTR(HBAR,SCRATCH,IWTRCNT,1)
            CALL LNDWTR(BVOL,SCRATCH,IWTRCNT,NZ)
            CALL LNDWTR(VDER,SCRATCH,IWTRCNT,NZ)
          ENDIF
          IF(NZ.GT.1)  CALL LNDWTR(HYDSAL,SCRATCH,IWTRCNT,NZ)
          IF(NZ.GT.1)  CALL LNDWTR(HYDTEMP,SCRATCH,IWTRCNT,NZ)
        ENDIF
      ENDIF



C        CHECK FOR NEGATIVE TEMPERATURES IN HYDRO MODEL
      DO IZ=1,NZ
       DO IY=1,NY
        DO IX=1,NX
          IF(HYDTEMP(IX,IY,IZ).LT.0.0) THEN
C           WRITE(OUT,1234)  TIME,IX,IY,IZ,HYDTEMP(IX,IY,IZ)
 1234       FORMAT(/10X,
     .       'WARNING: NEGATIVE TEMPERATURE PROVIDED BY THE HYDRODYNAMIC
     . MODEL AT TIME =',F10.5' DAYS'/10X,
     .       'IN MODEL SEGMENT(IX,IY,IZ) =',3I4,' TEMPERATURE OF',F7.2,
     .       ' DEG C MODIFIED TO 0.0 DEG C')
            HYDTEMP(IX,IY,IZ) = 0.0
          ENDIF
        ENDDO
       ENDDO
      ENDDO

C        POSSIBLE FUTURE CODE IF HYDRODYNAMIC RECORDS NOT SAVED AT HOURLY
C          OR CONSTANT INTERVALS
C     READ(30,END=920)  TMID
C     BACKSPACE(30)
C     NXHYDT=TMID-TMID0
      NXHYDTSECS=NXHYDTSECS+IHYDDTSECS
      NXHYDT=NXHYDTSECS/86400.0d0

C     If uncollapsed grid compute segment volumes  (m^3)

      IF(ICOLLOPT.EQ.0)  THEN
        DO 40 IZ=1,NZ
        DO 40 IY=1,NY
        DO 40 IX=1,NX
C        The following Commented code should not be active if the
C            hydrodynamic model maintains true continuity and time-step
C            round-off is not a problem
C         IF(INITB.GT.0)  THEN
C          DO 36 ISYS=1,NOSYS
C           CARAY(IX,IY,IZ,ISYS) = BVOL(IX,IY,IZ)*CARAY(IX,IY,IZ,ISYS)
C  36      CONTINUE
C         ENDIF
          BVOL(IX,IY,IZ)=DX(IX,IY)*DY(IX,IY)*(H(IX,IY)+ETA(IX,IY))
     .                   *DZ(IZ)+1.E-20
C         IF(INITB.GT.0)  THEN
C          DO 38 ISYS=1,NOSYS
C           CARAY(IX,IY,IZ,ISYS) = CARAY(IX,IY,IZ,ISYS)/BVOL(IX,IY,IZ)
C  38      CONTINUE
C         ENDIF
   40   CONTINUE
      ENDIF

      IF(LIST(1).EQ.1)  THEN
         WRITE(OUT,1400)
 1400    FORMAT(/49X,'SEGMENT VOLUMES (M^3)'//)
         CALL RCAPRNT(BVOL,1,NZ)
      ENDIF

C        IF UNYLAPSED GRID COMPUTE SURFACE AREAS AND
C        CROSS SECTIONAL AREAS

      IF(ICOLLOPT.EQ.0)  THEN
        DO 65 IY=1,NY
        DO 65 IX=1,NX
          XAZ(IX,IY) = DX(IX,IY)*DY(IX,IY)
   65   CONTINUE

        HYDDT2 = HYDDT/2.

        DO 70 IZ=1,NZ
        DO 70 IY=2,NY
        DO 70 IX=2,NX

        IF(FSM(IX,IY).GE.1. OR. FSM(IX,IY).LE.-2)THEN
          IF(FSM(IX-1,IY).GE.1. OR. FSM(IX-1,IY).LE.-2)THEN
            AVEH =(H(IX-1,IY)
     .           +ETA(IX-1,IY)+HYDDT2*DETA(IX-1,IY)*86400.d0
     .           +H(IX,IY)
     .           +ETA(IX,IY)+HYDDT2*DETA(IX,IY)*86400.0d0)/2.0d0
            XAX(IX,IY,IZ) = AVEH*DZ(IZ)*(DY(IX-1,IY)+DY(IX,IY))/2.
          ENDIF

          IF(FSM(IX,IY-1).GE.1. OR. FSM(IX,IY-1).LE.-2)THEN
            AVEH = (H(IX,IY-1)+
     .             ETA(IX,IY-1)+HYDDT2*DETA(IX,IY-1)*86400.0d0
     .           + H(IX,IY)
     .           +ETA(IX,IY)+HYDDT2*DETA(IX,IY)*86400.0d0)/2.0d0
            XAY(IX,IY,IZ) = AVEH*DZ(IZ)*(DX(IX,IY-1)+DX(IX,IY))/2.
          ENDIF
        ENDIF

   70   CONTINUE

C        COMPUTE AVERAGE DEPTHS FOR THE TIME-INTERVAL

        DO 75 IY=1,NY
        DO 75 IX=1,NX
          HBAR(IX,IY) = H(IX,IY)+ETA(IX,IY)
     .                 +HYDDT2*DETA(IX,IY)*86400.0d0
   75   CONTINUE
      ENDIF
     
C        COMPUTE ADVECTIVE AND DISPERSIVE FIELDS
C
      IF (NZ.EQ.1) THEN
        DO 85 IY=1,NY
        DO 85 IX=1,NX
            TPS(IX,IY)=SQRT(DX(IX,IY)*DX(IX,IY)+DY(IX,IY)*DY(IX,IY))
 85     CONTINUE
C
        FMIN=1.E20
        DO 86 IY=1,NY
        DO 86 IX=1,NX
            FMIN=AMIN1(FMIN,SNGL(TPS(IX,IY)))
 86     CONTINUE
      ENDIF

      IF(ICOLLOPT.EQ.0)  THEN
      DO 90 IZ=1,NZ
      DO 90 IY=2,NY
      DO 90 IX=2,NX
 
        IF(FSM(IX,IY).GE.1. OR. FSM(IX,IY).LE.-2)THEN
          IF (NZ.EQ.1) THEN
            IF(FSM(IX-1,IY).GE.1. OR. FSM(IX-1,IY).LE.-2)THEN
              EX(IX,IY,IZ)=1.0*FSM(IX,IY)*FSM(IX-1,IY)*TPS(IX,IY)/FMIN
            ENDIF
            IF(FSM(IX,IY-1).GE.1. OR. FSM(IX,IY-1).LE.-2)THEN
              EY(IX,IY,IZ)=1.0*FSM(IX,IY)*FSM(IX,IY-1)*TPS(IX,IY)/FMIN
            ENDIF
          ENDIF
 
          IF(FSM(IX-1,IY).GE.1. OR. FSM(IX-1,IY).LE.-2)THEN
            RX(IX,IY,IZ) = SCALRX*EX(IX,IY,IZ)*XAX(IX,IY,IZ)
     .                    /(0.5*(DX(IX-1,IY)+DX(IX,IY)))
          ENDIF

          IF(FSM(IX,IY-1).GE.1. OR. FSM(IX,IY-1).LE.-2)THEN
            RY(IX,IY,IZ) = SCALRY*EY(IX,IY,IZ)*XAY(IX,IY,IZ)
     .                    /(0.5*(DY(IX,IY-1)+DY(IX,IY)))
          ENDIF
        ENDIF

        IF(IZ.GT.1.AND.NZ.GT.1.AND.FSM(IX,IY).GE.1) THEN
           RZ(IX,IY,IZ) = SCALRZ*EZ(IX,IY,IZ)*XAZ(IX,IY)
     .                    /(0.5*(DZ(IZ)+DZ(IZ-1))*HBAR(IX,IY))
        ENDIF
   90 CONTINUE

C        COMPUTE VOLUME DERIVATIVE
C        (M^3/SEC * 86400. ==> M^3/DAY)

      DO 100 IZ=1,NZ
      DO 100 IY=1,NY
      DO 100 IX=1,NX
        VDER(IX,IY,IZ) = 86400.0d0*DZ(IZ)*DETA(IX,IY)*XAZ(IX,IY)
  100 CONTINUE
      ENDIF

C        SCALE ADVECTIVE FIELDS

      IF(ICOLLOPT.EQ.1)  THEN
      DO 105 IZ=1,NZ
      DO 105 IY=1,NY
      DO 105 IX=1,NX
        RX(IX,IY,IZ) = SCALRX*RX(IX,IY,IZ)
        RY(IX,IY,IZ) = SCALRY*RY(IX,IY,IZ)
        RZ(IX,IY,IZ) = SCALRZ*RZ(IX,IY,IZ)
  105 CONTINUE
      ENDIF
     
      IF(LIST(1).EQ.0)  GO TO 110

      WRITE(OUT,8000)
      WRITE(OUT,1500)
 1500 FORMAT(44X,'F L O W S   ( M * * 3 / S E C )'//)
      WRITE(OUT,1501)
 1501 FORMAT(20X,'HORIZONTAL OR X-DIRECTION COMPONENT'/)
      CALL RCAPRNT(QX,1,NZ)
      WRITE(OUT,1502)
 1502 FORMAT(20X,'LATERAL OR Y-DIRECTION COMPONENT'/)
      CALL RCAPRNT(QY,1,NZ)
      IF (NZ.GT.1) WRITE(OUT,1503)
 1503 FORMAT(20X,'VERTICAL OR Z-DIRECTION COMPONENT'/)
      IF (NZ.GT.1) CALL RCAPRNT(QZ,1,NZ)
      IF(LIST(1).EQ.1)  THEN
       WRITE(OUT,1550)
 1550  FORMAT(
     .  /20X,'DIFFUSER DISCHARGES'/ 8X,'IX   IY  FLOW RATE (M^3/SEC)')
       DO I=1,NODIFF
        WRITE(OUT,1600) IDD(I),JDD(I),QDIFF(I)
 1600   FORMAT(5X,2I5,2X,10E10.3)
       ENDDO
      ENDIF
      WRITE(OUT,8000)
      WRITE(OUT,2000)
 2000 FORMAT(28X,'E X C H A N G E   C O E F F I C I E N T S ',
     .    ' (M * * 3 / S E C )'//)
      WRITE(OUT,2001)
 2001 FORMAT(20X,'HORIZONTAL OR X-DIRECTION COMPONENT'/)
      CALL RCAPRNT(RX,1,NZ)
      WRITE(OUT,2002)
 2002 FORMAT(20X,'LATERAL OR Y-DIRECTION COMPONENT'/)
      CALL RCAPRNT(RY,1,NZ)
      IF (NZ.GT.1) WRITE(OUT,2003)
 2003 FORMAT(20X,'VERTICAL OR Z-DIRECTION COMPONENT'/)
      IF (NZ.GT.1) CALL RCAPRNT(RZ,1,NZ)
  
C        CONVERT FLOWS AND EXCHANGES TO M^3/DAY

  110 DO 120 IZ=1,NZ
      DO 120 IY=1,NY
      DO 120 IX=1,NX
          QX(IX,IY,IZ) = 86400.0d0*QX(IX,IY,IZ)
          QY(IX,IY,IZ) = 86400.0d0*QY(IX,IY,IZ)
          IF (NZ.GT.1) QZ(IX,IY,IZ) = 86400.0d0*QZ(IX,IY,IZ)
          RX(IX,IY,IZ) = 86400.0d0*RX(IX,IY,IZ)
          RY(IX,IY,IZ) = 86400.0d0*RY(IX,IY,IZ)
          IF (NZ.GT.1) RZ(IX,IY,IZ) = 86400.0d0*RZ(IX,IY,IZ)
          IF (NZ.GT.1.AND.IZ.EQ.1) RZ(IX,IY,IZ) = 0.
          EX(IX,IY,IZ) = 86400.0d0*EX(IX,IY,IZ)
          EY(IX,IY,IZ) = 86400.0d0*EY(IX,IY,IZ)
          IF (NZ.GT.1) EZ(IX,IY,IZ) = 86400.0d0*EZ(IX,IY,IZ)
          IF (FSM(IX,IY).LE.0) EZ(IX,IY,IZ)=0.0
  120 CONTINUE
      IF(IDIFFOPT.EQ.1)  THEN
        DO 130 IDIF=1,NODIFF
          QDIFF(IDIF) = 86400.0d0*QDIFF(IDIF)
  130   CONTINUE
      ENDIF

C        MAKE SURE R ARE ZERO FOR APPROPRIATE LAND AND BC SEGMENTS

      DO 150 IZ=1,NZ
      DO 150 IY=1,NY-1
      DO 150 IX=1,NX-1
         IF(FSM(IX,IY).NE.-1)  GO TO 150
         IF(IX.LT.NX .AND. FSM(IX+1,IY).EQ.1) RX(IX+1,IY,IZ)=0.
         IF(IX.GT.1 .AND. FSM(IX-1,IY).EQ.1) RX(IX,IY,IZ)=0.
         IF(IY.LT.NY .AND. FSM(IX,IY+1).EQ.1) RY(IX,IY+1,IZ)=0.
         IF(IY.GT.1 .AND. FSM(IX,IY-1).EQ.1) RY(IX,IY,IZ)=0.
  150 CONTINUE

C        COMPUTE CRITICAL INTEGRATION STEPSIZES

      DO 160 IZ=1,NZ
      DO 160 IY=1,NY-1
      DO 160 IX=1,NX-1
       IF(FSM(IX,IY).LE.0)  THEN
        TIMESTEP(IX,IY,IZ) = 999.999
       ELSE
        TIMESTEP(IX,IY,IZ) = RX(IX,IY,IZ) + RX(IX+1,IY,IZ)
     .                     + RY(IX,IY,IZ) + RY(IX,IY+1,IZ)
        IF(QX(IX,IY,IZ).LT.0.0)  TIMESTEP(IX,IY,IZ) =
     .                              TIMESTEP(IX,IY,IZ) - QX(IX,IY,IZ)
        IF(QX(IX+1,IY,IZ).GT.0.0)  TIMESTEP(IX,IY,IZ) =
     .                              TIMESTEP(IX,IY,IZ) + QX(IX+1,IY,IZ)
        IF(QY(IX,IY,IZ).LT.0.0)  TIMESTEP(IX,IY,IZ) =
     .                              TIMESTEP(IX,IY,IZ) - QY(IX,IY,IZ)
        IF(QY(IX,IY+1,IZ).GT.0.0)  TIMESTEP(IX,IY,IZ) =
     .                              TIMESTEP(IX,IY,IZ) + QY(IX,IY+1,IZ)
        IF(QZ(IX,IY,IZ).GT.0.0)  TIMESTEP(IX,IY,IZ) =
     .                              TIMESTEP(IX,IY,IZ) + QZ(IX,IY,IZ)
        IF(QZ(IX,IY,IZ+1).LT.0.0)  TIMESTEP(IX,IY,IZ) =
     .                              TIMESTEP(IX,IY,IZ) - QZ(IX,IY,IZ+1)
        TIMESTEP(IX,IY,IZ) = BVOL(IX,IY,IZ)/(TIMESTEP(IX,IY,IZ)+1.E-20)
       ENDIF
  160 CONTINUE
 
      IF(IDIAGDT.EQ.1 .AND. TIME.GT.0.0)  THEN

      IMIN=1
      DO 180 IZ=1,NZ
      DO 180 IY=1,NY
      DO 180 IX=1,NX
        IF(FSM(IX,IY).LE.0)  GO TO 180
        IF(IMIN.LT.21)  THEN
           MINDT(IMIN)=TIMESTEP(IX,IY,IZ)
           IMINX(IMIN) = IX
           IMINY(IMIN) = IY
           IMINZ(IMIN) = IZ
           IMIN = IMIN+1
        ELSE
           STEPMIN=TIMESTEP(IX,IY,IZ)
           IXMIN=IX
           IYMIN=IY
           IZMIN=IZ
           DO 170 IMINR=1,20
             IF(STEPMIN.LT.MINDT(IMINR))   THEN
              TEMPSTOR = MINDT(IMINR)
              MINDT(IMINR)=STEPMIN
              STEPMIN=TEMPSTOR
              ITEMPSTOR=IMINX(IMINR)
              IMINX(IMINR) = IXMIN
              IXMIN=ITEMPSTOR
              ITEMPSTOR=IMINY(IMINR)
              IMINY(IMINR) = IYMIN
              IYMIN=ITEMPSTOR
              ITEMPSTOR=IMINZ(IMINR)
              IMINZ(IMINR) = IZMIN
              IZMIN=ITEMPSTOR
             ENDIF
  170      CONTINUE
        ENDIF
  180 CONTINUE

        WRITE(OUT,4300)  TIME,NINT(24.*TIME),TMID
 4300 FORMAT(///5X,
     .  'TIMESTEP INFORMATION FOR WATER QUALITY MODEL TIME ='
     .  F8.2,' DAYS (= ',I7' HOURS) ... HYDRO MODEL TIME =',F8.2,
     .  ' DAYS')
      WRITE(OUT,4400)  (MINDT(IMIN),IMINX(IMIN),IMINY(IMIN),
     .   IMINZ(IMIN),IMIN=1,20)
 4400 FORMAT(   40X,' MINIMUM TIMESTEP CRITERIA'/
     .          40X,'STEPSIZE      IX  IY  IZ'/
     .          40X,' (DAYS)'/(40X,F7.5,4X,3I4))

      CALL MOVE(NX*NY*NZ,TIMESTEP,TIMESTEPS)
      CALL SORT(NX*NY*NZ,TIMESTEPS)
      IP1 = 0.01*TOTLWATR
      IP5 = 0.05*TOTLWATR
      IP10 = 0.10*TOTLWATR
      IP15 = 0.15*TOTLWATR
      IP20 = 0.20*TOTLWATR
      WRITE(OUT,4600)  1,TIMESTEPS(1),IP1,TIMESTEPS(IP1),
     .     IP5,TIMESTEPS(IP5),IP10,TIMESTEPS(IP10),
     .     IP15,TIMESTEPS(IP15),IP20,TIMESTEPS(IP20)
 4600 FORMAT(/25X,'DISTRIBUTION OF CRITICAL TIMESTEPS'/
     .  15X,' 0 PERCENT OR',I4,' SEGMENT HAS A TIMESTEP OF LESS THAN',
     .  F10.5,' DAYS'/
     .  15X,' 1 PERCENT OR',I4,' SEGMENTS HAVE A TIMESTEP OF LESS THAN',
     .  F8.5,' DAYS'/
     .  15X,' 5 PERCENT OR',I4,' SEGMENTS HAVE A TIMESTEP OF LESS THAN',
     .  F8.5,' DAYS'/
     .  15X,'10 PERCENT OR',I4,' SEGMENTS HAVE A TIMESTEP OF LESS THAN',
     .  F8.5,' DAYS'/
     .  15X,'15 PERCENT OR',I4,' SEGMENTS HAVE A TIMESTEP OF LESS THAN',
     .  F8.5,' DAYS'/
     .  15X,'20 PERCENT OR',I4,' SEGMENTS HAVE A TIMESTEP OF LESS THAN',
     .  F8.5,' DAYS')
      DO 220 I=1,20
  220 MINDT(I) = 100.

      ENDIF

C        FLOW BALANCE CHECK
      IF(INPCHCK.EQ.1)  THEN
       DO 230 IY=1,NY
        DO 230 IX=1,NX
         IF(FSM(IX,IY).EQ.1.)  THEN
          QBAL(IX,IY)=0.0
          TVDER(IX,IY)=0.0
          DO 225 IZ=1,NZ
           QBAL(IX,IY) = QBAL(IX,IY) + QX(IX,IY,IZ)-QX(IX+1,IY,IZ)
     .                               + QY(IX,IY,IZ)-QY(IX,IY+1,IZ)
           TVDER(IX,IY) = TVDER(IX,IY) + VDER(IX,IY,IZ)
  225     CONTINUE
          IF(NODIFF.GT.0)  THEN
           DO IDIF=1,NODIFF
            IXD=IDD(IDIF)
            IYD=JDD(IDIF)
            IF(IX.EQ.IXD .AND. IY.EQ.IYD)  THEN
             DO IZ=1,NZ
              QBAL(IX,IY) = QBAL(IX,IY)+VDDIST(IDIF,IZ)*QDIFF(IDIF)
             ENDDO
            ENDIF
           ENDDO
          ENDIF
          IF(ABS(QBAL(IX,IY)-TVDER(IX,IY)).LE.QTOL*ABS(TVDER(IX,IY)))
     .                                                    GO TO 230
          WRITE(OUT,5000)  IX,IY,TIME,NINT(24.*TIME),TMID,
     .                     QBAL(IX,IY),TVDER(IX,IY)
 5000     FORMAT(15X,'CONTINUITY ERROR FOUND IN CELL',I4,I4,
     .       ' AT WATER QUALITY MODEL TIME = ',F8.3,' DAYS (=',I7,
     .       ' HOURS)'/20X,'HYDRODYNAMIC TIME = ',F8.3,' DAYS',5X,
     .       'QBAL =',E11.4,'   VDER =',E11.4)
         ENDIF
  230  CONTINUE
      ENDIF

      IF(IDIAGDT.EQ.0 .AND. INPCHCK.EQ.0)   RETURN

      DO 250 IZ=1,NZ
       DO 250 IY=1,NY
        DO 250 IX=1,NX
         TIMESTEP(IX,IY,IZ) = 100.
  250 CONTINUE
   
      GO TO 30

C        EXIT POINT IF PERFORMING TIME-STEP DIAGNOSTIC RUN
  300 CALL EXIT

C        HYDRODYNAMIC FILE DOES NOT EXIST ...
C          ERROR MESSAGE AND EXIT
  900 WRITE(OUT,9000)   HYDFILNA(1)
 9000 FORMAT(///5X,'ERROR...THE FOLLOWING HYDRODYNAMIC FILE DOES NOT EXI
     .ST OR IS EMPTY'/5X,A40/5X,'RCA TERMINATED')
      CALL EXIT
C        DIFFUSER FILE DOES NOT EXIST ...
C          ERROR MESSAGE AND EXIT
  910 WRITE(OUT,9100)   DIFFILNA(1)
 9100 FORMAT(///5X,'ERROR...THE FOLLOWING DIFFUSER FILE DOES NOT EXIST O
     .R IS EMPTY'/5X,A40/5X,'RCA TERMINATED')
      CALL EXIT

C        END OF HYDRODYNAMIC FILE ENCOUNTERED ...
C          CHECK IF RECYCLE OPTION ENABLED
  920 IF(HYDCYCOPT.EQ.1)  THEN
C        IF SO REWIND FILE AND CONTINUE SIMULATION
        REWIND(30)
        IF(IDIFFOPT.EQ.1)  THEN
          REWIND(31)
          READ(31)
        ENDIF
        GO TO 30
C          OR IF "LINK-FILE" OPTION ENABLED
      ELSE IF(NOHYDFILNA.GT.1)  THEN
C        IF SO CLOSE PRESENT "GCM" FILE AND OPEN NEW "GCM" FILE
!
!----------------------------------------------------------------------
! Optionally close hydrodynamic file
!----------------------------------------------------------------------
!
        IF(NETCDFOPT.EQ.1) THEN
          status=nf90_close(ncID_his)
          CALL nccheck_status(status,HYDFILNA(IHYDFILE),RCANA)
        ELSE
          CLOSE(30)
        ENDIF

        IF(IDIFFOPT.EQ.1)  CLOSE(31)
        IHYDFILE=IHYDFILE+1
        IF(IHYDFILE.GT.NOHYDFILNA)  GO TO 970
        WRITE(OUT,1310)  HYDFILNA(IHYDFILE),DIFFILNA(IHYDFILE)
!
!----------------------------------------------------------------------
! Optionally open hydrodynamic file
!----------------------------------------------------------------------
!
        IF(NETCDFOPT.EQ.1) THEN
          status=nf90_open(TRIM(ADJUSTL(HYDFILNA(IHYDFILE))),
     .                   nf90_nowrite,ncID_his)
          CALL nccheck_status(status,HYDFILNA(IHYDFILE),RCANA)
          hydTflag=1       ! ready for beginning of new file
          ncIREC=1         ! first record of new output file
        ELSE
          OPEN(30,FILE=HYDFILNA(IHYDFILE),FORM='UNFORMATTED',
     .         STATUS='OLD',ERR=900)
          IF(IDIFFOPT.EQ.1)  THEN
            OPEN(31,FILE=DIFFILNA(IHYDFILE),
     .           FORM='UNFORMATTED',STATUS='OLD',ERR=910)
            READ(31)
          ENDIF
        ENDIF  
        GO TO 30

      ELSE IF(INPCHCK.EQ.1) THEN
        WRITE(OUT,9150)
 9150   FORMAT(//,20X,'HYDRODYNAMIC FILE CHECKING COMPLETED'//)
        RETURN
      ENDIF
C        ERROR MESSAGE AND EXIT
      WRITE(OUT,9200)   HYDFILNA(IHYDFILE),DIFFILNA(IHYDFILE)
 9200 FORMAT(///5X,'END-OF-FILE ENCOUNTERED DURING READ OF'/
     .     5X,'HYDRODYNAMIC FILE =',A40,' -or-'/
     .     5X,'DIFFUSER FILE =',A40/5X,' ...RCA TERMINATED')
      CALL EXIT

C        ERROR MESSAGE AND EXIT
  930 WRITE(OUT,9300)   HYDFILNA(IHYDFILE),DIFFILNA(IHYDFILE)
 9300 FORMAT(///5X,'APPARENT SYNCHRONIZATION ERROR BETWEEN HYDRODYNAMIC
     .AND DIFFUSER FILES'/5X,'END-OF-FILE ENCOUNTERED DURING READ OF DIF
     .FUSER FILE'/5X,A40/5X,' ...RCA TERMINATED')
      CALL EXIT

  950 CALL FMTER
      CALL EXIT

  960 WRITE(OUT,9600)
 9600 FORMAT(///10X,
     .  'ERROR...USER REQUESTED MORE HYDRODYNAMIC TIME BREAKS'/
     .  ' THEN DIMENSIONED FOR IN -RCACM- ... RCA TERMINATED')
      CALL EXIT
  970 IF(INPCHCK.EQ.1) RETURN
      WRITE(OUT,9700)
 9700 FORMAT(///
     .  5X,'ERROR...END-OF-FILE ENCOUNTERED IN LAST HYDRODYNAMIC FILE'/
     .  5X,'NO MORE USER-SPECIFIED HYDRODYNAMIC FILES EXIST'/
     .  5X,' ... RCA TERMINATED')
      CALL EXIT
  980 WRITE(OUT,9800)   NODIFF,MXDIFF
 9800 FORMAT(///
     .  5X,'ERROR...USER REQUESTED MORE DIFFUSERS (',I5,') THEN RCA DIME
     .NSIONED FOR (',I5,')'/5X,'CHANGE MXDIFF IN RCA AND RECOMPILE'/
     .  5X,'RCA TERMINATED'/)
      CALL EXIT
  990 WRITE(OUT,9900)   TMID,TMIDD
 9900 FORMAT(///
     .  5X,'ERROR...SYNC ERROR BETWEEN HYDRODYNAMIC AND DIFFUSER FILES'
     . /5X,'HYDRO TMID (',F10.4,') .NE. TO DIFFUSER TMID (',F10.4,')'
     . /5X,'RCA TERMINATED'/)
      CALL EXIT

      END 

      SUBROUTINE LNDWTR(ARRAY,SCRATCHA,ICNT,JCNT)
      INCLUDE 'RCACM'
      REAL   ARRAY(NX,NY,JCNT),SCRATCHA(ICNT,JCNT)
      READ(30)   SCRATCHA
      DO 50 KK = 1,JCNT
       DO 50 JJ = 1,NY
        DO 50 II = 1,NX
          ARRAY(II,JJ,KK) = 0.
  50  CONTINUE
      DO 100 IZ=1,JCNT
        DO 100 I=1,ICNT
          ARRAY(IWTRNDX(I),JWTRNDX(I),IZ) = SCRATCHA(I,IZ)
  100 CONTINUE
      RETURN
      END

      SUBROUTINE MOVE(N,FROM,TO)
      REAL  FROM(N),TO(N)
      DO 10 I=1,N
        TO(I) = FROM(I)
   10 CONTINUE
      RETURN
      END

      SUBROUTINE SORT(N,RA)
      DIMENSION RA(N)
      L=N/2+1
      IX=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IX)
          RA(IX)=RA(1)
          IX=IX-1
          IF(IX.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IX)THEN
          IF(J.LT.IX)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IX+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END
