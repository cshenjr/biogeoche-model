      SUBROUTINE RCA11
C 
C        RCA11 UPDATES THE BOUNDARY CONCENTRATIONS 
C
!=======================================================================
! The file was revised to read boundary file in NetCDF format          !
!                                 YUN LI, UMCES/HPL Feb-19-2011        !
!====================================================================== 
!
      USE netcdf
! 
      SAVE
      INCLUDE  'RCACM'
      INCLUDE  'NetCDFCM'
      REAL  SLBC(NSL),SIGDEPTH(NZ),SIGBC(NZ)
!
!-----------------------------------------------------------------------
! Name of this rca file
!-----------------------------------------------------------------------
!
      RCANA='rca11.f'
!
!-----------------------------------------------------------------------
! Optially read boundary time 
!-----------------------------------------------------------------------
!
      IF(BCFILNA.EQ.'NULL')  GO TO 930
      IF(IBCPWLOPT.EQ.1)  THEN
        OLDBCT = NXBCT
        IF(IBNRYRDOPTS(1).EQ.0) THEN
          READ(32,1050,ERR=950,END=960)  NXBCTSECS
 1050     FORMAT(10X,I10)
        ELSEIF(IBNRYRDOPTS(1).EQ.1) THEN
          READ(32,ERR=950,END=960)  NXBCTSECS
        ELSEIF(IBNRYRDOPTS(1).EQ.2) THEN
          bryTflag=bryTflag+1
          CALL ncbry_time(ncID_bry)
        ENDIF
        NXBCTSECS = ISCALBC*NXBCTSECS
        NXBCT=NXBCTSECS/86400.
      ENDIF

      IF(IBCPWLOPT.EQ.0)  THEN
        WRITE(OUT,1000)   NXBCT
 1000   FORMAT(//40X,'BOUNDARY CONDITIONS AT TIME =',F8.2,' DAYS')
      ELSE
        IF(INITB.NE.0) WRITE(OUT,1000)   OLDBCT
      ENDIF
      IF(IBCPWLOPT.EQ.0)  GO TO 50

      DO 20 ISYS=1,NOSYS
        IF(NOBC(ISYS).EQ.0)  GO TO 20
        DO 10 I=1,NOBC(ISYS)
          SBC(I,ISYS) = BBC(I,ISYS)
   10   CONTINUE
   20 CONTINUE

   50 DO 200 ISYS=1,NOSYS 
      IF(NOBC(ISYS).EQ.0)   GO TO 200

      IF(INITB.NE.0) THEN
        WRITE(OUT,1200) ISYS
 1200   FORMAT(41X,'BOUNDARY CONCENTRATIONS FOR SYSTEM',I3) 
        WRITE(OUT,1250)   SYNAME(ISYS)
 1250   FORMAT(41X,11('-'),3X,A8,3X,12('-')) 
      ENDIF

      IF(IBCOPT.EQ.4)  GO TO 100

C        SIGMA-LEVELS
      IF(LIST(2).EQ.1)   WRITE(OUT,2000)
 2000 FORMAT(4(6X,'ROW COL LYR',6X,'CONC',3X))
!
!-----------------------------------------------------------------------
! Optionally read boundary indices and values
!-----------------------------------------------------------------------
!
      IF (IBNRYRDOPTS(1).EQ.0) THEN
        READ(32,2100,ERR=950,END=960)  (BBC(I,ISYS),I=1,NOBC(ISYS))
      ELSEIF (IBNRYRDOPTS(1).EQ.1) THEN
        READ(32,ERR=950,END=960)  (BBC(I,ISYS),I=1,NOBC(ISYS))
      ELSEIF (IBNRYRDOPTS(1).EQ.2) THEN
        CALL ncbry_IBCBBC(ncID_bry,bryTflag)
      ENDIF
 2100 FORMAT(10X,7F10.0) 
      IF(IBCPWLOPT.EQ.0) THEN
        IF(LIST(2).EQ.1)   WRITE(OUT,2200)  ((IBC(J,I,ISYS),J=1,3),
     .       BBC(I,ISYS),I=1,NOBC(ISYS))
 2200   FORMAT(4(5X,3I4,E13.3))
      ELSE
        IF(INITB.NE.0 .AND. LIST(2).EQ.1) WRITE(OUT,2200)
     .               ((IBC(J,I,ISYS),J=1,3),SBC(I,ISYS),I=1,NOBC(ISYS))
      ENDIF

      DO 75 I=1,NOBC(ISYS)
        BBC(I,ISYS) = SCALBC(ISYS)*BBC(I,ISYS) 
   75 CONTINUE

      GO TO 190

C        STANDARD LEVELS
  100 NBCRD = NOBC(ISYS)/NZ
      JBC = 0
      DO 135 I=1,NBCRD
        IX = IBC(1,NZ*(I-1)+1,ISYS)
        IY = IBC(2,NZ*(I-1)+1,ISYS)
        IF (IBNRYRDOPTS(1).EQ.0) THEN
          READ(32,1100,ERR=950,END=960)  (SLBC(J),J=1,NOBCSL(I,ISYS))
        ELSE
          READ(32,ERR=950,END=960)  (SLBC(J),J=1,NOBCSL(I,ISYS))
        ENDIF
 1100   FORMAT(10X,7F10.0) 
        DO 125 IZ=1,NZ
          SIGDEPTH(IZ) = -ZZ(IZ)*HBAR(IX,IY)
  125   CONTINUE
        CALL SINTER(SLDEPTH(1,ISYS),SLBC,SIGDEPTH,SIGBC
     .     ,NOBCSL(I,ISYS),NZ)
        DO 130 IZ=1,NZ
          JBC = JBC + 1
          BBC(JBC,ISYS) = SCALBC(ISYS)*SIGBC(IZ)
  130   CONTINUE
  135 CONTINUE

      IF(IBCPWLOPT.EQ.0) THEN
        IF(LIST(2).EQ.1)  THEN
          WRITE(OUT,3300) 
 3300     FORMAT(//31X,
     .     'BOUNDARY CONCENTRATIONS AFTER SIGMA-LEVEL TRANSFORMATION'/)
          WRITE(OUT,1700) 
 1700     FORMAT(12X,4(8X,'BC  ROW COL LYR ')) 
          WRITE(OUT,1900)  (BBC(I,ISYS),(IBC(J,I,ISYS),J=1,3),
     .                                                  I=1,NOBC(ISYS))
 1900     FORMAT(11X,E12.4,3I4,E12.4,3I4,E12.4,3I4,E12.4,3I4)
        ENDIF
      ELSE
        IF(INITB.NE.0 .AND. LIST(2).EQ.1)  THEN
          WRITE(OUT,3300) 
          WRITE(OUT,1700) 
          WRITE(OUT,1900)  (SBC(I,ISYS),(IBC(J,I,ISYS),J=1,3),
     .                                                  I=1,NOBC(ISYS))
        ENDIF
      ENDIF

  190 IF(LIST(2).EQ.1)  WRITE(OUT,3400)
 3400 FORMAT(//)
  200 CONTINUE

      IF(INPCHCK.EQ.1)  RETURN

      IF(IBCPWLOPT.EQ.0)  GO TO 250

C        CALCULATE SLOPES FOR PIECEWISE LINEAR BC OPTION
      DO 230 ISYS=1,NOSYS
        IF(NOBC(ISYS).EQ.0)  GO TO 230
        DO 215 I=1,NOBC(ISYS)
          SBC(I,ISYS) = (SBC(I,ISYS)-BBC(I,ISYS))/(OLDBCT-NXBCT)
  215   CONTINUE
  230 CONTINUE
      OLDBCT = NXBCT
C        MOVE BOUNDARY CONDITIONS TO -CARAY-
  250 DO 300 ISYS=1,NOSYS
       IF(NOBC(ISYS).EQ.0)  GO TO 300
       DO 275 I=1,NOBC(ISYS)
        IF(IBCPWLOPT.EQ.0)  THEN
          CARAY(IBC(1,I,ISYS),IBC(2,I,ISYS),IBC(3,I,ISYS),ISYS)=
     .                        BBC(I,ISYS)
        ELSE
          CARAY(IBC(1,I,ISYS),IBC(2,I,ISYS),IBC(3,I,ISYS),ISYS)=
     .                       (TIME-NXBCT)*SBC(I,ISYS) + BBC(I,ISYS)
        ENDIF
  275  CONTINUE
  300 CONTINUE
!
!-----------------------------------------------------------------------
! Optionally read boundary time 
!-----------------------------------------------------------------------
!
      IF(IBCPWLOPT.EQ.0)  THEN
        IF (IBNRYRDOPTS(1).EQ.0) THEN
          READ(32,1050,ERR=950,END=960)  NXBCTSECS
        ELSEIF (IBNRYRDOPTS(1).EQ.1) THEN
          READ(32,ERR=950,END=960)  NXBCTSECS
        ELSEIF (IBNRYRDOPTS(1).EQ.2) THEN
          bryTflag=bryTflag+1
          CALL ncbry_time(ncID_bry)
        ENDIF
        NXBCTSECS = ISCALBC*NXBCTSECS
        NXBCT=NXBCTSECS/86400.
      ENDIF

  930 RETURN

  950 IN=32
      BACKSPACE(32)
      CALL FMTER
      CALL EXIT 
  960 WRITE(OUT,9600)
 9600 FORMAT(//5X,
     .  'END OF FILE ENCOUNTERED WHILE READING BOUNDARY CONDITION FILE')
      INPCHCK=2
      RETURN
      END 
