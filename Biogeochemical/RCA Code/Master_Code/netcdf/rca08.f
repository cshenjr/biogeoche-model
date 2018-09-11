      SUBROUTINE  RCA08 
C 
C      RCA08 IS USED TO SET UP SYNCRONIZE DATA FILES FOR "B" VECTORS
C            IF SIMULATION STARTS OTHER THAN TIME EQUAL TO ZERO 
C 
!=======================================================================
! The file was revised to syncronize data in NetCDF format             ! 
!                                 YUN LI, UMCES/HPL Feb-18-2011        !
!=======================================================================
!
      USE netcdf
!
      SAVE
      INCLUDE 'RCACM' 
      INCLUDE 'NetCDFCM'
 
      REAL   TVAL(600)
!
!-----------------------------------------------------------------------
! Name of this rca file
!-----------------------------------------------------------------------
!
      RCANA='rca08.f'
      WRITE(OUT,'(////1X,119(''*'')//)')
      WRITE(OUT,'(42X,''S Y N C R O N I Z E   D A T A''//)')

C     Transport

   10 IF(TZERO.GT.NXHYDT)   THEN
C        Dummy reads for TMID/QX/QY/QZ/EX/EY/EZ/ETA/DETA/HYDTEMP/HYDSAL
!
!----------------------------------------------------------------------
! Optionally read model time 'TMID'
!----------------------------------------------------------------------
!
       IF(NETCDFOPT.EQ.1) THEN
         CALL ncset_oceantime(ncID_his,hydTflag,TMID,not_enough_length)
         IF(not_enough_length) GO TO 900  ! no record found
         TMID=TMID/86400.                 ! change the unit, sec -> day   
         hydTflag=hydTflag+1
       ELSE
         READ(30,END=900) TMID
         IF(LNDWTROPT.EQ.0)  THEN
C          READ FULL LAND/WATER GRID
           READ(30) 
           READ(30) 
           IF (NZ.GT.1) READ(30) 
           READ(30) 
           READ(30) 
           IF (NZ.GT.1) READ(30) 
           IF(IECOMVER.EQ.0 .AND. NZ.GT.1)  READ(30)
           READ(30) 
           READ(30) 
           IF (NZ.GT.1) READ(30) 
           IF (NZ.GT.1) READ(30) 
         ELSE
C          Read water only grid and put into full grid
           READ(30)
           READ(30)
           IF (NZ.GT.1)    READ(30)
C          Check if uncollapsed/collapsed read
           IF(ICOLLOPT.EQ.0)  THEN
             READ(30)
             READ(30)
             IF (NZ.GT.1)    READ(30)
             IF(IECOMVER.EQ.0 .AND. NZ.GT.1)  READ(30)
             READ(30)
             READ(30)
           ELSE
             READ(30)
             READ(30)
             IF (NZ.GT.1)    READ(30)
             READ(30)
             READ(30)
             READ(30)
           ENDIF
           READ(30)
           READ(30)
         ENDIF
C        Dummy reads for TMID/QDIFF
         IF(IDIFFOPT.EQ.1)  READ(31,END=900)
       ENDIF
!
!----------------------------------------------------------------------
! Show time info when rca syncronizes hydro transport
!----------------------------------------------------------------------
!
         WRITE(OUT,'(A,A,f12.5,A,f12.5,A,I4)') TRIM(ADJUSTL(RCANA)),
     .             ':   Pass Hydro data until ', TZERO,
     .             ',   Hydro clock now is ', NXHYDT,
     .             ',   Hydro record is ', hydTflag

         NXHYDTSECS=NXHYDTSECS+IHYDDTSECS
         NXHYDT=NXHYDTSECS/86400.
         GO TO 10
      ENDIF

      CALL RCA03A

C     Boundary Conditions

      IF(BCFILNA.EQ.'NULL') GO TO 175
  120 IF(IBCOPT.EQ.1)  GO TO 175
  130 IF(TZERO.GE.NXBCT)  THEN

 2300  FORMAT(10X,7F10.0)
       CALL RCA11
!
!----------------------------------------------------------------------
! Show time info when rca syncronizes boundary condition
!----------------------------------------------------------------------
!
       RCANA='rca08.f'
       WRITE(OUT,'(A,A,f12.5,A,f12.5,A,I4)') TRIM(ADJUSTL(RCANA)),
     .             ':   Pass Bndry data until ', TZERO,
     .             ',   Bndry clock now is ', NXBCT,
     .             ',   Bndry record is ', bryTflag 
       GO TO 120

      ENDIF

C     Put Boundary Conditions into -CARAY-

      DO 170 ISYS=1,NOSYS
        IF(NOBC(ISYS).EQ.0)  GO TO 170
        DO 160 I=1,NOBC(ISYS)
          CARAY(IBC(1,I,ISYS),IBC(2,I,ISYS),IBC(3,I,ISYS),ISYS)
     .    =BBC(I,ISYS)
  160   CONTINUE
  170 CONTINUE

C        POINT SOURCE LOADS

  175 IF(PSFILNA.EQ.'NULL') GO TO 205
  180 IF(IPSOPT.EQ.1)  GO TO 205
  185 IF(TZERO.GE.NXPST)  THEN

 2750   FORMAT(10X,7F10.0)
        CALL RCA10(SPS,BPS,MXWK,NOPS,NXPST,33,IPSPWLOPT,SCALPS)
        GO TO 185

      ENDIF

C     Nonpoint source loads

  205 IF(NPSFILNA.EQ.'NULL') GO TO 235
  210 IF(INPSOPT.EQ.1)  GO TO 235
  215 IF(TZERO.GE.NXNPST)  THEN

        CALL RCA10(SNPS,BNPS,MXWK,NONPS,NXNPST,34,INPSPWLOPT,SCALNPS)
        GO TO 215

      ENDIF

C     Fall-line loads

  235 IF(FLFILNA.EQ.'NULL') GO TO 275
  240 IF(IFLOPT.EQ.1)  GO TO 275
  245 IF(TZERO.GE.NXFLT)  THEN

        CALL RCA10(SFL,BFL,MXWK,NOFL,NXFLT,35,IFLPWLOPT,SCALFL)
        GO TO 245

      ENDIF

C     Atmospheric loads

  275 IF(ATMFILNA.EQ.'NULL') GO TO 305
  280 IF(IATMOPT.EQ.1)  GO TO 305
  285 IF(TZERO.GE.NXATMT)  THEN

        CALL RCA10(SATM,BATM,NX*NY,NOATM,NXATMT,36,IATMPWLOPT,
     .             SCALATM)
        GO TO 285

      ENDIF

C     Misc time functions

  305 IF(PCFILNA.EQ.'NULL') RETURN
  310 IF(NOFUNC.EQ.0)  RETURN
      DO I=1,NOFUNC
       IF(TZERO.GE.NXFUNT(I))  THEN
         DO 320 IBRK=1,MXFUNCT
           IF(TZERO.LE.TIMEMTF(I,IBRK))   GO TO 330
  320    CONTINUE
         WRITE(OUT,2000)  TZERO,I
 2000    FORMAT(///5X,'ERROR ... REQUESTED START TIME (TZERO) =',F10.3,
     .    ' GREATER THAN MAXIMUM TIME BREAK SPECIFIED FOR',
     .    ' MISCELLANEOUS TIME FUNCTION NUMBER',I3/5X,'RCA TERMINATED')
  330    ITIMF(I) = IBRK
       ENDIF
      ENDDO

      RETURN

C     End of hydrodynamic file encountered ...
C     Check if recycle option enabled
  900 IF(HYDCYCOPT.EQ.1)  THEN
C        IF SO REWIND FILE AND CONTINUE SIMULATION
        REWIND(30)
        IF(IDIFFOPT.EQ.1) THEN
          REWIND(31)
          READ(31)
        ENDIF
        GO TO 10
C          OR IF "LINK-FILE" OPTION ENABLED
      ELSE IF(HYDCYCOPT.EQ.2)  THEN
C        IF SO CLOSE PRESENT "GCM-" FILES AND OPEN NEW "GCM-" FILES
!
!----------------------------------------------------------------------
! Optionally close hydrodynamic file
!----------------------------------------------------------------------
!
        IF(NETCDFOPT.EQ.1) THEN
          status=nf90_close(ncID_his)
          CALL nccheck_status(status,HYDFILNA(IHYDFILE),SUBNA)
        ELSE
          CLOSE(30)
        ENDIF

        IF(IDIFFOPT.EQ.1)  CLOSE(31)
        IHYDFILE=IHYDFILE+1
        IF(IHYDFILE.GT.MXHYDFILES)  GO TO 970
        WRITE(OUT,1310)  HYDFILNA(IHYDFILE),DIFFILNA(IHYDFILE)
 1310   FORMAT(//20X,'OPENING HYDRODYNAMIC TRANPORT FIELD FILE =',A40/
     .       16X,'AND OPENING DIFFUSER DISCHARGE FILE =',A40)
!
!----------------------------------------------------------------------
! Optionally open hydrodynamic file
!----------------------------------------------------------------------
!
        IF(NETCDFOPT.EQ.1) THEN 
          status=nf90_open(TRIM(ADJUSTL(HYDFILNA(IHYDFILE))),
     .                   nf90_nowrite,ncID_his)
          CALL nccheck_status(status,HYDFILNA(IHYDFILE),SUBNA)
          hydTflag=1       ! ready for beginning of new file
        ELSE
          OPEN(30,FILE=HYDFILNA(IHYDFILE),FORM='UNFORMATTED',
     .         STATUS='OLD',ERR=980)
          IF(IDIFFOPT.EQ.1) THEN
            OPEN(31,FILE=DIFFILNA(IHYDFILE),
     .           FORM='UNFORMATTED',STATUS='OLD',ERR=990)
            READ(31)
          ENDIF
        ENDIF
        GO TO 10
      ELSE
C        ERROR MESSAGE AND EXIT
        WRITE(OUT,9100)   HYDFILNA(IHYDFILE)
 9100   FORMAT(///5X,'END-OF-FILE ENCOUNTERED DURING READ OF'/
     .     5X,A40)
      ENDIF
      CALL EXIT

  950 CALL FMTER
      CALL EXIT 
  970 WRITE(OUT,9700)
 9700 FORMAT(///
     .  5X,'ERROR...END-OF-FILE ENCOUNTERED IN LAST HYDRODYNAMIC FILE'/
     .  5X,'NO MORE USER-SPECIFIED HYDRODYNAMIC FILES EXIST'/
     .  5X,' ... RCA TERMINATED')
      CALL EXIT
  980 WRITE(OUT,9800)  HYDFILNA(IHYDFILE)
 9800 FORMAT(///
     .  5X,'ERROR ENCOUNTERED OPENING HYDRODYNAMIC FILE ',A40/
     .  5X,' ... RCA TERMINATED')
      CALL EXIT
  990 WRITE(OUT,9900)  DIFFILNA(IHYDFILE)
 9900 FORMAT(///
     .  5X,'ERROR ENCOUNTERED OPENING DIFFUSER FILE ',A40/
     .  5X,' ... RCA TERMINATED')
      CALL EXIT
      END 
