      SUBROUTINE RCA06
C 
C        RCA06 READS PARAMS, CONST, MISC TIME FUNCTIONS AND MISC FILE NAMES
C 
!=======================================================================
! The file was revised to read ROMS 2D/3D params in NetCDF format      !
!                                 YUN LI, UMCES/HPL May-12-2011        !
!=======================================================================
!
      USE netcdf
! 
      SAVE
      INCLUDE   'RCACM'
      INCLUDE   'NetCDFCM'
      CHARACTER  PNAME*20,CNAME(MXCONS)*10
      REAL   PSCAL(30)
!
!-----------------------------------------------------------------------
! Name of this rca file
!-----------------------------------------------------------------------
!
      RCANA='rca06.f'

C         P A R A M E T E R S 

      WRITE(OUT,8000) 
 8000 FORMAT(/////1X,119('*')//)
      TWARPTVF='SECS'

      READ(IN,1000)  COMMENT
 1000 FORMAT(A)
      READ(IN,1004)  PCFILNA,IBNRYRDOPT
 1004 FORMAT(A40,I10)
      IF(PCFILNA.EQ.'NULL')  THEN
        WRITE(OUT,1006)
        GO TO 450
      ENDIF
 1006 FORMAT(///30X,
     .      'USER SPECIFIES THAT THERE IS NO PC FILE TO BE READ')
      WRITE(OUT,1005)  PCFILNA
 1005 FORMAT(10X,'OPENING PCFILNA = ',A40//)
      IF(IBNRYRDOPT.EQ.0)  THEN
        OPEN(37,FILE=PCFILNA,FORM='FORMATTED')
      ELSE
        OPEN(37,FILE=PCFILNA,FORM='UNFORMATTED')
      ENDIF

C        DEPTH INDEPENDENT (2-D)

      IF(IBNRYRDOPT.EQ.0)  THEN
        READ(37,1000) COMMENT
        READ(37,1001,ERR=950) NOPAM, I2DNCOPT
 1001   FORMAT(5I10) 
      ELSE
        READ(37) NOPAM 
      ENDIF
      IF(NOPAM.GT.MXPARM2D)  GO TO 960
      IF(NOPAM.EQ.0)   GO TO 120
!
!--------------------------------------------------------------
! Optionally open model 2D parameter file
!-------------------------------------------------------------
!
      IF(I2DNCOPT.EQ.0) THEN
      ELSEIF(I2DNCOPT.EQ.1) THEN
       READ(37,1000) COMMENT
       READ(37,'(A80)',ERR=950) PARAMFILNA
       status=nf90_open(TRIM(ADJUSTL(PARAMFILNA)),nf90_nowrite,ncID_prm)
       CALL nccheck_status(status,PARAMFILNA,RCANA)
      ELSE
       WRITE(OUT,'( A)') 'ERROR: WRONG I2DNCOPT,'
       WRITE(OUT,'(2A)') '       I2DNCOPT=0, format same as ',
     .                          'RCA IBNRYRDOPT defined'
       WRITE(OUT,'( A)') '       I2DNCOPT=1, read NetCDF 2D-params'
       CALL EXIT
      ENDIF

      IF(IBNRYRDOPT.EQ.0)  THEN
        READ(37,1000) COMMENT
        READ(37,1100,ERR=950)    (PSCAL(I),I=1,NOPAM) 
 1100   FORMAT(8F10.0)
      ELSE
        READ(37)    (PSCAL(I),I=1,NOPAM) 
      ENDIF
      WRITE(OUT,1002)    (PSCAL(I),I=1,NOPAM) 
 1002 FORMAT(//51X,'LAYER INDEPENDENT (2-D) PARAMETERS'
     .   //8X,'SCALE'/6X,'FACTORS',5(6X,E14.5)/13X,5(6X,E14.5))
      DO 100 I=1,NOPAM
        IF(IBNRYRDOPT.EQ.0)  THEN
          READ(37,1200,ERR=950)   PNAME
 1200     FORMAT(A20)
        ELSE
          READ(37)   PNAME
        ENDIF
        WRITE(OUT,1110)   PNAME
 1110   FORMAT(//20X,A20)
!
!--------------------------------------------------------------
! Optionally read model 2D parameters
!      I2DNCOPT=0, format same as RCA IBNRYRDOPT defined
!      I2DNCOPT=1, RCA reads NetCDF file 
!-------------------------------------------------------------
!
        IF(I2DNCOPT.EQ.0) THEN
         DO 74 IY=1,NY
           IF(IBNRYRDOPT.EQ.0)  THEN
             READ(37,1100,ERR=950)   (PARAM2D(IX,IY,I),IX=1,NX)
           ELSE
             READ(37)   (PARAM2D(IX,IY,I),IX=1,NX)
           ENDIF
   74    CONTINUE
        ELSEIF(I2DNCOPT.EQ.1) THEN
         CALL ncprm_set_n2dvar(ncID_prm,id2Dpm(I),PARAM2D(:,:,I))
        ENDIF

        DO 75 IY=1,NY
          IF(LIST(4).EQ.1)  THEN
           WRITE(OUT,1052) 
 1052      FORMAT(10X,'COL',43X,'ROW'/25X,'2',15X,'3',14X,'...'/)
           WRITE(OUT,1202)   IY,(PARAM2D(IX,IY,I),IX=1,NX)
 1202      FORMAT(10X,I4,1X,8(2X,E13.5)/15X,8(2X,E13.5)/
     .            15X,8(2X,E13.5))
          ENDIF
C        INCORPORATE SCALE FACTOR
          DO 50 IX=1,NX
            PARAM2D(IX,IY,I) = PARAM2D(IX,IY,I)*PSCAL(I) 
   50     CONTINUE
   75   CONTINUE
  100 CONTINUE
!
!--------------------------------------------------------------
!  Close Parameter File after 2D parameters are done
!-------------------------------------------------------------
!
      IF(I2DNCOPT.EQ.1) THEN
        status=nf90_close(ncID_prm)
        CALL nccheck_status(status,PARAMFILNA,RCANA)
      ENDIF

C        DEPTH DEPENDENT (3-D)

  120 IF(IBNRYRDOPT.EQ.0)  THEN
        READ(37,1000) COMMENT
        READ(37,1001,ERR=950) NOPAM, I3DNCOPT 
      ELSE
        READ(37) NOPAM 
      ENDIF
      IF(NOPAM.GT.MXPARM3D)  GO TO 962
      IF(NOPAM.EQ.0)   GO TO 220
!
!--------------------------------------------------------------
! Optionally open model 3D parameter file
!-------------------------------------------------------------
!
      IF(I3DNCOPT.EQ.0) THEN
      ELSEIF(I3DNCOPT.EQ.1) THEN
       READ(37,1000) COMMENT
       READ(37,'(A80)',ERR=950) PARAMFILNA
       status=nf90_open(TRIM(ADJUSTL(PARAMFILNA)),nf90_nowrite,ncID_prm)
       CALL nccheck_status(status,PARAMFILNA,RCANA)
      ELSE
       WRITE(OUT,'( A)') 'ERROR: WRONG I3DNCOPT,'
       WRITE(OUT,'(2A)') '       I3DNCOPT=0, format same as ',
     .                          'RCA IBNRYRDOPT defined'
       WRITE(OUT,'( A)') '       I3DNCOPT=1, read NetCDF 3D-params'
       CALL EXIT
      ENDIF
      IF(IBNRYRDOPT.EQ.0)  THEN
        READ(37,1000) COMMENT
        READ(37,1100,ERR=950)    (PSCAL(I),I=1,NOPAM) 
      ELSE
        READ(37)    (PSCAL(I),I=1,NOPAM) 
      ENDIF
      WRITE(OUT,1003)    (PSCAL(I),I=1,NOPAM) 
 1003 FORMAT(//51X,'LAYER DEPENDENT (3-D) PARAMETERS'
     .   //8X,'SCALE'/6X,'FACTORS',5(6X,E14.5)/13X,5(6X,E14.5))
      DO 200 I=1,NOPAM
        IF(IBNRYRDOPT.EQ.0)  THEN
          READ(37,1200,ERR=950)   PNAME
        ELSE
          READ(37)   PNAME
        ENDIF
        WRITE(OUT,1110)   PNAME
!
!--------------------------------------------------------------
! Optionally read model 3D parameters
!      I3DNCOPT=0, format same as RCA IBNRYRDOPT defined
!      I3DNCOPT=1, RCA reads NetCDF file 
!-------------------------------------------------------------
!
        IF(I3DNCOPT.EQ.0) THEN
          DO 179 IZ=1,NZ
          DO 179 IY=1,NY
            IF(IBNRYRDOPT.EQ.0)  THEN
              READ(37,1100,ERR=950)   (PARAM3D(IX,IY,IZ,I),IX=1,NX)
            ELSE
              READ(37)   (PARAM3D(IX,IY,IZ,I),IX=1,NX)
            ENDIF
  179     CONTINUE
        ELSEIF(I3DNCOPT.EQ.1) THEN
          CALL ncprm_set_n3dvar(ncID_prm,id3Dpm(I),PARAM3D(:,:,:,I))
        ENDIF
        DO 180 IZ=1,NZ
          DO 160 IY=1,NY
            IF(LIST(4).EQ.1)  THEN
              WRITE(OUT,1053) 
 1053         FORMAT(4X,'LAYER COL',43X,'ROW'/25X,'2',15X,'3',
     .               14X,'...'/)
              WRITE(OUT,1203)   IZ,IY,(PARAM3D(IX,IY,IZ,I),IX=1,NX)
 1203         FORMAT(6X,I2,2X,I4,1X,8(2X,E13.5)/15X,8(2X,E13.5)/
     .          15X,8(2X,E13.5)/15X,8(2X,E13.5)/15X,8(2X,E13.5)/
     .          15X,8(2X,E13.5)/15X,8(2X,E13.5)/15X,8(2X,E13.5))
            ENDIF
C        INCORPORATE SCALE FACTOR
            DO 140 IX=1,NX
              PARAM3D(IX,IY,IZ,I) = PARAM3D(IX,IY,IZ,I)*PSCAL(I) 
  140       CONTINUE
  160     CONTINUE
  180   CONTINUE
  200 CONTINUE
!
!--------------------------------------------------------------
!  Close Parameter File after 3D parameters are done
!-------------------------------------------------------------
!
      IF(I3DNCOPT.EQ.1) THEN
        status=nf90_close(ncID_prm)
        CALL nccheck_status(status,PARAMFILNA,RCANA)
      ENDIF

C         C O N S T A N T S 

  220 IF(IBNRYRDOPT.EQ.0)  THEN
        READ(37,1000) COMMENT
        READ(37,1001,ERR=950)   NOCONS 
      ELSE
        READ(37)   NOCONS 
      ENDIF
      IF(NOCONS.GT.MXCONS)  GO TO 964
      IF(NOCONS.EQ.0)   GO TO 250
      IF(IBNRYRDOPT.EQ.0)  THEN
        DO I=1,NOCONS,8
         JL=I
         JM=MIN(I+7,NOCONS)
         READ(37,1300,ERR=950)   (CNAME(J),J=JL,JM)
         READ(37,1350,ERR=950)   (CONST(J),J=JL,JM)
 1300    FORMAT(8A10)
 1350    FORMAT(8F10.0)
        ENDDO
      ELSE
        READ(37)   (CNAME(I),CONST(I),I=1,NOCONS)
      ENDIF
      IF(LIST(4).EQ.1)  THEN
        WRITE(OUT,1400)   (CNAME(I),CONST(I),I=1,NOCONS)
 1400   FORMAT(//55X,'CONSTANTS'//8X,4(2X,A10,'=',E13.5)/ 
     .        (10X,A10,'=',E13.5,2X,A10,'=',E13.5,2X,A10,'=',E13.5,2X,
     .         A10,'=',E13.5)) 
      ENDIF

C        M I S C   T I M E   F U N C T I O N S

  250 IF(IBNRYRDOPT.EQ.0)  THEN
        READ(37,1000) COMMENT
        READ(37,1001,ERR=950)  NOFUNC,ITVFPWLOPT,ITVNCOPT
      ELSE
        READ(37)  NOFUNC,ITVFPWLOPT
      ENDIF
      IF(NOFUNC.GT.MXFUNC)  GO TO 966
      IF(NOFUNC.EQ.0)   GO TO 350
      WRITE(OUT,1600)
 1600 FORMAT(//45X,'MODEL DEPENDENT TIME FUNCTIONS')
      IF(ITVFPWLOPT.EQ.0)   THEN
        WRITE(OUT,1610)
 1610   FORMAT(/37X,'THE USER SELECTED STEP FUNCTION INTERPOLATION')
      ELSE
        WRITE(OUT,1620)
 1620   FORMAT(/35X,'THE USER SELECTED PIECEWISE LINEAR INTERPOLATION')
      ENDIF
      WRITE(OUT,1650) 
 1650 FORMAT(/6X,5(7X,'VAL(T)',6X,'T',2X)/) 
!
!--------------------------------------------------------------
! Optionally open model time-dependent parameter file
!-------------------------------------------------------------
!
      IF(ITVNCOPT.EQ.0) THEN
      ELSEIF(ITVNCOPT.EQ.1) THEN
       READ(37,1000) COMMENT
       READ(37,'(A80)',ERR=950) PARAMFILNA
       status=nf90_open(TRIM(ADJUSTL(PARAMFILNA)),nf90_nowrite,ncID_prm)
       CALL nccheck_status(status,PARAMFILNA,RCANA)
      ELSE
       WRITE(OUT,'( A)') 'ERROR: WRONG ITVNCOPT,'
       WRITE(OUT,'(2A)') '       ITVNCOPT=0, format same as ',
     .                          'RCA IBNRYRDOPT defined'
       WRITE(OUT,'(2A)') '       ITVNCOPT=1, read NetCDF ',
     .                          'time-dependent params'
       CALL EXIT
      ENDIF

      NXCALL13T=999.
      DO 300 I=1,NOFUNC
!
!--------------------------------------------------------------
! Optionally read model time-dependent parameters
!      ITVNCOPT=0, format same as RCA IBNRYRDOPT defined
!      ITVNCOPT=1, RCA reads NetCDF file 
!-------------------------------------------------------------
!
        IF(ITVNCOPT.EQ.0) THEN
          IF(IBNRYRDOPT.EQ.0)  THEN
            READ(37,1000) COMMENT
            READ(37,1700,ERR=950)   PNAME,NOBRK,TWARPTVF
 1700       FORMAT(A10,I10,6X,A4)
            IF(NOBRK.GT.MXFUNCT)  GO TO 968
            READ(37,1800,ERR=950)   (VALMTF(I,J),TIMEMTF(I,J),J=1,NOBRK)
 1800       FORMAT(8F10.0) 
          ELSE
            READ(37)   PNAME,NOBRK,TWARPTVF
            IF(NOBRK.GT.MXFUNCT)  GO TO 968
            READ(37)   (VALMTF(I,J),TIMEMTF(I,J),J=1,NOBRK)
          ENDIF
        ELSEIF(ITVNCOPT.EQ.1) THEN
          READ(37,1000) COMMENT
          READ(37,1700,ERR=950)   PNAME
          CALL ncprm_set_t1dvar(ncID_prm,idTVpm_time(I),idTVpm(I),
     .                          NOBRK,TWARPTVF,VALMTF(I,:),TIMEMTF(I,:))
          IF(NOBRK.GT.MXFUNCT)  GO TO 968
        ENDIF
        ISCALTVF=IUNITCHECK(TWARPTVF,'TWARPTVF')
c       DO J=1,NOBRK
c        TVAL(J)=ISCALTVF*TVAL(J)
c       ENDDO
          IF(LIST(4).EQ.1)  THEN
            WRITE(OUT,2000)   PNAME,(VALMTF(I,J),TIMEMTF(I,J),J=1,NOBRK)
 2000       FORMAT(1X,A10,5(E15.5,F7.2)/(6X,E15.5,F7.2,E15.5,F7.2
     .          ,E15.5,F7.2,E15.5,F7.2,E15.5,F7.2)) 
          ENDIF
          IF(ITVFPWLOPT.EQ.0)  THEN
            MFUNC(I)=0.
            BFUNC(I)=VALMTF(I,1) 
          ELSE
            MFUNC(I)=(VALMTF(I,2)-VALMTF(I,1))
     .              /(TIMEMTF(I,2)-TIMEMTF(I,1))
            BFUNC(I)=VALMTF(I,2)
          ENDIF
          NXFUNT(I) = TIMEMTF(I,2)
          ITIMF(I) = 2
          IF(NXFUNT(I).LT.NXCALL13T)  NXCALL13T = NXFUNT(I)
  300 CONTINUE
!
!--------------------------------------------------------------
!  Close Parameter File after time-depedent parameters are done
!-------------------------------------------------------------
!
      IF(ITVNCOPT.EQ.1) THEN
        status=nf90_close(ncID_prm)
        CALL nccheck_status(status,PARAMFILNA,RCANA)
      ENDIF

  350 READ(37,1000) COMMENT
      READ(37,1001,ERR=950)  NOKINFILNA
      IF(NOKINFILNA.GT.MXKINFILES) GO TO 970
      IF(NOKINFILNA.EQ.0) GO TO 450
      WRITE(OUT,2100)  NOKINFILNA
 2100 FORMAT(//28X,
     . 'THE FOLLOWING',I3,' FILES WILL BE UTILIZED BY THE KINETIC SUBR')
      DO 375 I=1,NOKINFILNA
       READ(37,1004,ERR=950)  KINFILNA(I)
       WRITE(OUT,2200)  KINFILNA(I)
 2200  FORMAT(10X,A40)
  375 CONTINUE
  
  450 RETURN

  950 IN=37
      CALL FMTER
      CALL EXIT 
  960 WRITE(OUT,9600)  
 9600 FORMAT(//11X,'ERROR...USER REQUESTING MORE 2-D PARAMETERS THAN DIM
     .ENSIONED FOR IN -RCACM-'/11X,'RCA TERMINATED'//)
      CALL EXIT 
  962 WRITE(OUT,9602)  
 9602 FORMAT(//11X,'ERROR...USER REQUESTING MORE 3-D PARAMETERS THAN DIM
     .ENSIONED FOR IN -RCACM-'/11X,'RCA TERMINATED'//)
      CALL EXIT 
  964 WRITE(OUT,9604)  
 9604 FORMAT(//11X,'ERROR...USER REQUESTING MORE CONSTANTS THAN DIMENSIO
     .NED FOR IN -RCACM-'/11X,'RCA TERMINATED'//)
      CALL EXIT 
  966 WRITE(OUT,9606)  
 9606 FORMAT(//11X,'ERROR...USER REQUESTING MORE TIME FUNCTIONS THAN DIM
     .ENSIONED FOR IN -RCACM-'/11X,'RCA TERMINATED'//)
      CALL EXIT 
  968 WRITE(OUT,9608)  
 9608 FORMAT(//11X,'ERROR...USER REQUESTING MORE TIME BREAKS THAN DIMENS
     .IONED FOR IN -RCACM-'/11X,'RCA TERMINATED'//)
      CALL EXIT 
  970 WRITE(OUT,9610)  
 9610 FORMAT(//11X,'ERROR...USER REQUESTING MORE KINETIC INPUT FILE NAME
     .S THAN DIMENSIONED FOR IN -RCACM-'/11X,'RCA TERMINATED'//)
      CALL EXIT 
      END 
