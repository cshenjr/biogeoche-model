      SUBROUTINE RCA07
C 
C        RCA07 READS THE INITIAL CONDITIONS
C 
!=======================================================================
! The file was revised to read initial conditions in NetCDF format     !
!                                 YUN LI, UMCES/HPL Oct-13-2008        !
!=======================================================================
!
      USE netcdf
!
      SAVE
      INCLUDE  'RCACM'
      INCLUDE  'NetCDFCM'

      REAL   CARAYSL(NX,NY,NSL),SIGDEPTH(NZ)
      REAL   CSTD(NSL),CSIGMA(NZ)
!
!-----------------------------------------------------------------------
! Name of this rca file
!-----------------------------------------------------------------------
!
      RCANA='rca07.f'

      WRITE(OUT,8000) 
 8000 FORMAT(////1X,119('*')//)
      WRITE(OUT,1050) 
 1050 FORMAT(42X,'I N I T I A L   C O N D I T I O N S'//)

      READ(IN,1000)  COMMENT
 1000 FORMAT(A)
      READ(IN,1150)  ICFILNA,IBNRYRDOPT
 1150 FORMAT(A40,I10)
      IF(ICFILNA.EQ.'NULL' .OR. ICFILNA.EQ.'null' .OR. 
     .                          ICFILNA.EQ.'Null')  THEN
        ICFILNA='NULL'
        WRITE(OUT,1156)
        GO TO 183
      ENDIF
 1156 FORMAT(/36X,
     .   'USER DID NOT SPECIFY AN INITIAL FILE TO BE READ')
      IF(IBNRYRDOPT.EQ.0)  THEN
        OPEN(38,FILE=ICFILNA,FORM='FORMATTED')
      ELSEIF(IBNRYRDOPT.EQ.1)  THEN
        OPEN(38,FILE=ICFILNA,FORM='UNFORMATTED')
!
! Open NetCDF format initial condition file
!
      ELSEIF(IBNRYRDOPT.EQ.2)  THEN
        status=nf90_open(TRIM(ADJUSTL(ICFILNA)),nf90_nowrite,ncID_ini)
        CALL nccheck_status(status,ICFILNA,RCANA)
        CALL ncini_time(ncID_ini)  ! reset TZERO
        IF(TZERO.LT.HYDTZERO) THEN
          WRITE(OUT,'(A,f12.5,A,f12.5,A/7X,10A)') 
     .                 'ERROR: the 1st hydro time = ', HYDTZERO, 
     .                 ' DAYS, but TZERO = ', TZERO, ' DAYS',
     .                 TRIM(ADJUSTL(RCANA)),
     .                 ' RCA model failed to start, because the',
     .                 ' available hydrodynamic time in NetCDF files',
     .                 ' is larger than TZERO'
          CALL EXIT
        ENDIF
        IF(TZERO.LT.BCTZERO) THEN
          WRITE(OUT,'(A,f12.5,A,f12.5,A/7X,10A)') 
     .                 'ERROR: the 1st boundary time = ', BCTZERO, 
     .                 ' DAYS, but TZERO = ', TZERO, ' DAYS',
     .                 TRIM(ADJUSTL(RCANA)),
     .                 ' RCA model failed to start, because the',
     .                 ' available boundary time in NetCDF files',
     .                 ' is larger than TZERO'
          CALL EXIT
        ENDIF
      ELSE
        WRITE(OUT,'(3X,A,I3)'), 'ERROR: Invalid IBNRYRDOPT ', IBNRYRDOPT
        WRITE(OUT,'(3X,A)'),    '       0 -- formatted   or'
        WRITE(OUT,'(3X,A)'),    '       1 -- unformatted or'
        WRITE(OUT,'(3X,A)'),    '       2 -- NetCDF'
        CALL EXIT
      ENDIF


      IF(CYCLE.EQ.0)   THEN

C        READ SIGMA-LEVEL/STANDARD LEVEL OPTION
       IF(IBNRYRDOPT.EQ.0)  THEN
         READ(38,1000)  COMMENT
         READ(38,1175,ERR=938)  ICOPT
       ELSEIF(IBNRYRDOPT.EQ.1)  THEN
         READ(38)  ICOPT
!
! NetCDF uses sigma-level input 
!
       ELSEIF(IBNRYRDOPT.EQ.2)  THEN
         ICOPT=0
       ENDIF
 1175  FORMAT(I10)
       IF(ICOPT.EQ.0)  THEN
C        SIGMA-LEVEL INPUT
!
!----------------------------------------------
! NetCDF only reads SIGMA-LEVEL input
!----------------------------------------------
!
        IF(IBNRYRDOPT.EQ.2)  THEN
         CALL ncini_CARAY(ncID_ini)
        ELSE
         DO 100 ISYS=1,NOSYS 
          WRITE(OUT,1400)   ISYS
          WRITE(OUT,2150)   SYNAME(ISYS)
          IF(IBNRYRDOPT.EQ.0)   READ(38,1000)  COMMENT
           DO 100 IZ=1,NZ
            DO 100 IY=1,NY
            IF(IBNRYRDOPT.EQ.0)  THEN
              READ(38,1200,ERR=940)   (CARAY(IX,IY,IZ,ISYS),IX=1,NX) 
 1200         FORMAT(10X,7F10.0)
            ELSEIF(IBNRYRDOPT.EQ.1)  THEN
              READ(38)   (CARAY(IX,IY,IZ,ISYS),IX=1,NX) 
            ENDIF
  100    CONTINUE
        ENDIF

       ELSE

C        STANDARD LEVEL INPUT
         DO 170 ISYS=1,NOSYS
           WRITE(OUT,1400)   ISYS
 1400      FORMAT(/44X,'INITIAL CONDITIONS FOR SYSTEM',I3)
           WRITE(OUT,2150)   SYNAME(ISYS)
 2150      FORMAT(41X,11('-'),3X,A8,3X,11('-')/) 
C          READ STANDARD LEVEL DEPTHS
           IF(IBNRYRDOPT.EQ.0)  THEN
             READ(38,1000)  COMMENT
             READ(38,1175,ERR=938)  NLVLS
             IF(NLVLS.GT.NSL)  GO TO 970
             READ(38,1000)  COMMENT
             READ(38,1200,ERR=938)  (SLDEPTH(I,ISYS),I=1,NLVLS)
           ELSE
             READ(38)  NLVLS
             IF(NLVLS.GT.NSL)  GO TO 970
             READ(38)  (SLDEPTH(I,ISYS),I=1,NLVLS)
           ENDIF
           IF(LIST(5).EQ.1)   WRITE(OUT,1300)  ISYS,
     .          (I,I=1,10),(SLDEPTH(I,ISYS),I=1,NLVLS)
 1300      FORMAT(42X,'STANDARD LEVEL DEPTHS FOR SYSTEM',I3/
     .         17X,10(I2,8X)/(11X,10F10.2))
           DO 105 I=1,NLVLS
             SLDEPTH(I,ISYS) = -SLDEPTH(I,ISYS)
  105      CONTINUE

           IF(LIST(5).EQ.1)  WRITE(OUT,1350)   (I,I=1,10)
 1350      FORMAT(/43X,'CONCENTRATIONS AT STANDARD LEVELS'/
     .        2X,'ROW COL',8X,10(I2,8X))

           DO 125 IZ=1,NLVLS
            DO 115 IY=1,NY
             IF(IBNRYRDOPT.EQ.0)  THEN
               READ(38,1200,ERR=945)  (CARAYSL(IX,IY,IZ),IX=1,NX)
             ELSE
               READ(38)  (CARAYSL(IX,IY,IZ),IX=1,NX)
             ENDIF
  115       CONTINUE
             IF(LIST(5).EQ.1)   CALL RCAPRNT(CARAYSL(1,1,1),IZ,IZ)
  125      CONTINUE

           DO 150 IY=1,NY
            DO 150 IX=1,NX
             DO 130 IZ=1,NZ
              SIGDEPTH(IZ) = -ZZ(IZ)*HBAR(IX,IY)
              CSIGMA(IZ) = 0.0
  130        CONTINUE
             DO 140 IZ=1,NLVLS
              CSTD(IZ) = CARAYSL(IX,IY,IZ)
  140        CONTINUE

             CALL SINTER(SLDEPTH(1,ISYS),CSTD,SIGDEPTH
     .             ,CSIGMA,NLVLS,NZ)

             DO 145 IZ=1,NZ
              CARAY(IX,IY,IZ,ISYS) = CSIGMA(IZ)
  145        CONTINUE
  150      CONTINUE

  170    CONTINUE
       ENDIF

      ELSE
!
!---------------------------------------------------------------------
! Optionally read initial condition when restart
!---------------------------------------------------------------------
!
         IF(NETCDFOPT.EQ.1) THEN
! NetCDF uses sigma-level input
           CALL ncini_CARAY(ncID_ini)
         ELSE
           READ(15)  CARAY
           REWIND(15)
         ENDIF
      ENDIF

      IF(ICOPT.EQ.1)  WRITE(OUT,1380) 
 1380 FORMAT(//33X,'INITIAL CONDITIONS AFTER TRANSFORMATION TO SIGMA-LEV
     .ELS'/)

      IF(LIST(5).EQ.1)  THEN
       DO 180 ISYS=1,NOSYS
        WRITE(OUT,1400)   ISYS
        WRITE(OUT,2150)   SYNAME(ISYS)
        CALL RCAPRNT(CARAY(1,1,1,ISYS),1,NZ)
 180   CONTINUE
      ENDIF

 183  CONTINUE
      DO 185 ISYS=1,NOSYS
      DO 185 IZ=1,NZ
      DO 185 IY=1,NY
      DO 185 IX=1,NX
        IF(FSM(IX,IY).EQ.0.) CARAY(IX,IY,IZ,ISYS)=0.0
        IF(FSM(IX,IY).EQ.-1. .OR. FSM(IX,IY).EQ.-2.) THEN
          IF(NOBC(ISYS).EQ.0) THEN
            CARAY(IX,IY,IZ,ISYS)=0.0
          ELSE
            DO I=1,NOBC(ISYS)
             IF(IX.EQ.IBC(1,I,ISYS) .AND. IY.EQ.IBC(2,I,ISYS) .AND.
     .          IZ.EQ.IBC(3,I,ISYS) )GO TO 185
            ENDDO
            CARAY(IX,IY,IZ,ISYS)=0.0
          ENDIF
        ENDIF
  185 CONTINUE

C        MOVE BOUNDARY CONDITIONS TO -CARAY-
      DO 200 ISYS=1,NOSYS
      IF(NOBC(ISYS).EQ.0)  GO TO 200
      DO 190 I=1,NOBC(ISYS)
       CARAY(IBC(1,I,ISYS),IBC(2,I,ISYS),IBC(3,I,ISYS),ISYS)=BBC(I,ISYS)
  190 CONTINUE
  200 CONTINUE

 
      READ(IN,1000)  COMMENT
      READ(IN,1200,ERR=950)   (CMAX(I),I=1,NOSYS) 
      READ(IN,1000)  COMMENT
      READ(IN,1200,ERR=950)   (CMIN(I),I=1,NOSYS) 
      WRITE(OUT,3000)   (I,CMAX(I),CMIN(I),I=1,NOSYS) 
 3000 FORMAT(///32X,'STABILITY AND ACCURACY CRITERIA FOR NUMERICAL SOLUT
     .ION'//33X,'SYSTEM   MAXIMUM CONCENTRATION',10X,'EPSILON'/ 
     .   (34X,I3,10X,E11.4,14X,E11.4))

      RETURN

  940 IF(MOD(NX,7).EQ.0)  THEN
        LINES = 2+(ISYS-1)*NZ*(NY*(NX/7)+1)+(IY)*(NX/7)+1
        LINEE = 2+(ISYS-1)*NZ*(NY*(NX/7)+1)+(IY+1)*(NX/7)+1
      ELSE
        LINES = 2+(ISYS-1)*NZ*(NY*(NX/7+1)+1)+(IY)*(NX/7+1)+1
        LINEE = 2+(ISYS-1)*NZ*(NY*(NX/7+1)+1)+(IY+1)*(NX/7+1)+1
      ENDIF
      WRITE(OUT,9400)   ICFILNA,ISYS,LINES,LINEE
 9400 FORMAT(///
     . 25X,'ERROR ... ENCOUNTERED READING INITIAL CONDITIONS FILE ',A40/
     . 30X,'ERROR ENCOUNTERED READING INITIAL CONDITIONS FOR SYSTEM',I3/
     . 30X,'BETWEEN INPUT RECORDS',I6,' AND',I6)
      IN=38
      GO TO 950
  945 IF(MOD(NX,7).EQ.0)  THEN
        LINES = 2+(ISYS-1)*NLVLS*(NY*(NX/7)+1)+(IY)*(NX/7)+1
        LINEE = 2+(ISYS-1)*NLVLS*(NY*(NX/7)+1)+(IY+1)*(NX/7)+1
      ELSE
        LINES = 2+(ISYS-1)*NLVLS*(NY*(NX/7+1)+1)+(IY)*(NX/7+1)+1
        LINEE = 2+(ISYS-1)*NLVLS*(NY*(NX/7+1)+1)+(IY+1)*(NX/7+1)+1
      ENDIF
      WRITE(OUT,9400)   ICFILNA,ISYS,LINES,LINEE
  938 IN=38
  950 CALL FMTER
      CALL EXIT 
  970 WRITE(OUT,9700)  
 9700 FORMAT(//11X,'ERROR...USER REQUESTING MORE STANDARD LEVELS THAN DI
     .MENSIONED FOR IN -RCACM-'/11X,'RCA TERMINATED'//)
      CALL EXIT 
      END 
