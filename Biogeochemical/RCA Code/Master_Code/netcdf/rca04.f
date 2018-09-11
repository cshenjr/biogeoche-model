      SUBROUTINE RCA04
C 
C        RCA04 READS THE BOUNDARY CONCENTRATIONS 
C
!=======================================================================
! The file was revised to read boundary file in NetCDF format          !
!                                 YUN LI, UMCES/HPL Feb-19-2011        !
!======================================================================= 
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
      RCANA='rca04.f'

      WRITE(OUT,8000) 
 8000 FORMAT(////1X,119('*')///)
      TWARPBC='SECS'

      READ(IN,1000)   COMMENT
 1000 FORMAT(A)
      READ(IN,1005)  BCFILNA,IBNRYRDOPTS(1)
 1005 FORMAT(A40,I10)
      IF(BCFILNA.EQ.'NULL' .OR. BCFILNA.EQ.'null' .OR. 
     .                          BCFILNA.EQ.'Null')  THEN
        BCFILNA='NULL'
        WRITE(OUT,1006)
 1006   FORMAT(10X,
     ,   'USER DID NOT SPECIFY A BOUNDARY CONDITION FILE TO BE READ')
        DO ISYS=1,NOSYS 
          NOBC(ISYS)=0
        ENDDO
        GO TO 930
      ENDIF
      WRITE(OUT,1010)  BCFILNA
 1010 FORMAT(10X,'OPENING BCFILNA = ',A40)
!
!-----------------------------------------------------------------------
! Optionally open boundary file and get BC options
!-----------------------------------------------------------------------
!
      IF(IBNRYRDOPTS(1).EQ.0)  THEN
        OPEN(32,FILE=BCFILNA,FORM='FORMATTED')
        READ(32,1000)  COMMENT
        READ(32,1020,ERR=950) IBCOPT,IBCPWLOP
      ELSEIF(IBNRYRDOPTS(1).EQ.1)  THEN
        OPEN(32,FILE=BCFILNA,FORM='UNFORMATTED')
        READ(32) IBCOPT,IBCPWLOPT
      ELSEIF(IBNRYRDOPTS(1).EQ.2)  THEN
        status=nf90_open(TRIM(ADJUSTL(BCFILNA)),nf90_nowrite,ncID_bry)
        CALL nccheck_status(status,BCFILNA,RCANA)
        CALL ncbry_BCOPT(ncID_bry,idIBCOPT,IBCOPT)
        CALL ncbry_BCOPT(ncID_bry,idIBCPWLOPT,IBCPWLOPT)
      ELSE
        WRITE(OUT,'(3X,A,I3)'),'ERROR: Invalid IBNRYRDOPTS ',IBNRYRDOPTS
        WRITE(OUT,'(3X,A)'),   '       0 -- FORMATTED  or'
        WRITE(OUT,'(3X,A)'),   '       1 -- BINDARY or'
        WRITE(OUT,'(3X,A)'),   '       2 -- NetCDF'
        CALL EXIT
      ENDIF

 1020 FORMAT(3I10)
      WRITE(OUT,1030)
 1030 FORMAT(//42X,'BOUNDARY CONDITION OPTIONS SELECTED')
      IF(IBCOPT.EQ.1 .OR. IBCOPT.EQ.2)  WRITE(OUT,1040)
 1040 FORMAT(36X,'BOUNDARY CONCENTRATIONS INPUT USING SIGMA-LEVELS')
      IF(IBCOPT.EQ.3 .OR. IBCOPT.EQ.4)  WRITE(OUT,1050)
 1050 FORMAT(34X,'BOUNDARY CONCENTRATIONS INPUT USING STANDARD LEVELS')
      IF(IBCPWLOPT.EQ.0)  WRITE(OUT,1080)
 1080 FORMAT(45X,'STEP FUNCTION OPTION SELECTED')
      IF(IBCPWLOPT.EQ.1)  WRITE(OUT,1090)
 1090 FORMAT(43X,'PIECEWISE LINEAR OPTION SELECTED')

!
!-----------------------------------------------------------------------
! Optionally read boundary time
!-----------------------------------------------------------------------
!
      IF(IBCOPT.EQ.2 .OR. IBCOPT.EQ.4)  THEN
        IF(IBNRYRDOPTS(1).EQ.0)  THEN
          READ(32,1000)  COMMENT
          READ(32,1110,ERR=950)  NXBCTSECS,TWARPBC
 1110     FORMAT(10X,I10,6X,A4)
        ELSEIF(IBNRYRDOPTS(1).EQ.1)  THEN
          READ(32)  NXBCTSECS,TWARPBC
        ELSEIF(IBNRYRDOPTS(1).EQ.2)  THEN
          bryTflag=1  ! scan bry_FILE from the 1st record
          CALL ncbry_time(ncID_bry)
        ENDIF

 1100   FORMAT(10X,7F10.0) 
        ISCALBC=IUNITCHECK(TWARPBC,'TWARPBC ')
        NXBCTSECS=ISCALBC*NXBCTSECS
        NXBCT=NXBCTSECS/86400.
        WRITE(OUT,1120)   NXBCT
 1120   FORMAT(//39X,'BOUNDARY CONDITIONS AT TIME =',F8.2,' DAYS')
        BCTZERO=NXBCT
      ENDIF

      NOBCALL = 0
      DO 500 ISYS=1,NOSYS 
!-----------------------------------------------------------------------
! Optionally read number of boundary cells
!-----------------------------------------------------------------------
!
      IF(IBNRYRDOPTS(1).EQ.0)  THEN
        READ(32,1000)  COMMENT
        READ(32,1020,ERR=950) NOBC(ISYS)
      ELSEIF(IBNRYRDOPTS(1).EQ.1)  THEN
        READ(32) NOBC(ISYS)
      ELSEIF(IBNRYRDOPTS(1).EQ.2)  THEN
        CALL ncbry_NOBC
      ENDIF

      IF(NOBC(ISYS).EQ.0)   THEN
        WRITE(OUT,1200)    ISYS 
 1200   FORMAT(///42X,'NO BOUNDARY CONDITIONS FOR SYSTEM',I3) 
        WRITE(OUT,1250)   SYNAME(ISYS)
 1250   FORMAT(41X,11('-'),3X,A8,3X,11('-')) 
        GO TO 500 
      ENDIF

      IF(NOBC(ISYS).GT.MXBC)  GO TO 960

      WRITE(OUT,1300) ISYS
 1300 FORMAT(/41X,'BOUNDARY CONCENTRATIONS FOR SYSTEM',I3) 
      WRITE(OUT,1250)   SYNAME(ISYS)
      WRITE(OUT,1400) NOBC(ISYS)
 1400 FORMAT(48X,'NUMBER OF BC''S READ',I5)

!
!-----------------------------------------------------------------------
! Read scale number of boundary value
!-----------------------------------------------------------------------
!
      IF(IBNRYRDOPTS(1).EQ.0)  THEN
        READ(32,1000)  COMMENT
        READ(32,1500,ERR=950)  SCALBC(ISYS)
 1500   FORMAT(E10.3) 
      ELSEIF(IBNRYRDOPTS(1).EQ.1)  THEN
        READ(32)  SCALBC(ISYS)
      ELSEIF(IBNRYRDOPTS(1).EQ.2)  THEN
        SCALBC(ISYS)=Fscale(idTvar(ISYS))
      ENDIF

      WRITE(OUT,1600)  SCALBC(ISYS)
 1600 FORMAT(47X,'SCALE FACTOR',E13.3/) 

      GO TO (100,200,300,300),IBCOPT

C        CONSTANT - SIGMA LEVEL
  100 IF(LIST(2).EQ.1)  WRITE(OUT,1700) 
 1700 FORMAT(12X,4(8X,'BC  ROW COL LYR ')) 
      IF(IBNRYRDOPTS(1).EQ.0)  THEN
        READ(32,1000)  COMMENT
        READ(32,1800,ERR=950)  (BBC(I,ISYS),(IBC(J,I,ISYS),J=1,3),
     .     I=1,NOBC(ISYS)) 
 1800   FORMAT(4(F10.0,1X,3I3)) 
      ELSE
        READ(32)  (BBC(I,ISYS),(IBC(J,I,ISYS),J=1,3),I=1,NOBC(ISYS)) 
      ENDIF
      IF(LIST(2).EQ.1)  WRITE(OUT,1900)  (BBC(I,ISYS),(IBC(J,I,ISYS),
     .   J=1,3),I=1,NOBC(ISYS))
 1900 FORMAT(11X,E12.4,3I4,E12.4,3I4,E12.4,3I4,E12.4,3I4)

      DO 115 I=1,NOBC(ISYS)
      IX=IBC(1,I,ISYS)
      IY=IBC(2,I,ISYS)
C        CHECK IF VALID BOUNDARY SEGMENT NUMBER
      IF(FSM(IX,IY).EQ.-1 .OR. FSM(IX,IY).EQ.-2)  GO TO 105
       WRITE(OUT,1950)  ISYS,I,IX,IY,IBC(3,I,ISYS)
 1950  FORMAT(//10X,'ERROR IN BOUNDARY INPUT FILE FOR SYSTEM ',I4//
     .   10X,'BOUNDARY NUMBER ',I3,'  IBC(1)=',I3,'  IBC(2)=',I3,
     .   '  IBC(3)=',I3,' IS NOT A VALID BOUNDARY SEGMENT'/
     .   10X,'RCA TERMINATED')
       CALL EXIT
  105 IF(NOBCALL.EQ.0) THEN
         NOBCALL = 1
         IBCALL(1,1) = IBC(1,I,ISYS)
         IBCALL(2,1) = IBC(2,I,ISYS)
         IBCALL(3,1) = IBC(3,I,ISYS)
      ELSE
         DO 110 J=1,NOBCALL
         IF(IBC(1,I,ISYS).EQ.IBCALL(1,J) .AND.
     .        IBC(2,I,ISYS).EQ.IBCALL(2,J) .AND.
     .          IBC(3,I,ISYS).EQ.IBCALL(3,J))  GO TO 113
  110    CONTINUE
         NOBCALL = NOBCALL+1
         IF(NOBCALL.GT.MXBC)  GO TO 960
         IBCALL(1,NOBCALL) = IBC(1,I,ISYS)
         IBCALL(2,NOBCALL) = IBC(2,I,ISYS)
         IBCALL(3,NOBCALL) = IBC(3,I,ISYS)
      ENDIF
C          APPLY SCALE FACTOR
  113 BBC(I,ISYS) = SCALBC(ISYS)*BBC(I,ISYS)
  115 CONTINUE

      GO TO 500

C        TIME-VARYING - SIGMA LEVEL
  200 IF(LIST(2).EQ.1)   WRITE(OUT,2000)  
 2000 FORMAT(//4(6X,'ROW COL LYR',6X,'CONC',3X)/)
!
!-----------------------------------------------------------------------
! Read boundary indices and values (NetCDF uses mere sigma-level)
!-----------------------------------------------------------------------
!
      IF(IBNRYRDOPTS(1).EQ.0)  THEN
        READ(32,1000)   COMMENT
        READ(32,2100,ERR=950)  ((IBC(J,I,ISYS),J=1,3),I=1,NOBC(ISYS))
 2100   FORMAT(10X,I4,2I3,I4,2I3,I4,2I3,I4,2I3,I4,2I3,I4,2I3,I4,2I3)
        READ(32,1100,ERR=950)  (BBC(I,ISYS),I=1,NOBC(ISYS))
      ELSEIF(IBNRYRDOPTS(1).EQ.1)  THEN
        READ(32)  ((IBC(J,I,ISYS),J=1,3),I=1,NOBC(ISYS))
        READ(32)  (BBC(I,ISYS),I=1,NOBC(ISYS))
      ELSEIF(IBNRYRDOPTS(1).EQ.2)  THEN
        CALL ncbry_IBCBBC(ncID_bry)
      ENDIF

      IF(LIST(2).EQ.1)   WRITE(OUT,2200)  ((IBC(J,I,ISYS),J=1,3),
     .     BBC(I,ISYS),I=1,NOBC(ISYS))
 2200 FORMAT(4(5X,3I4,E13.3))

      DO 225 I=1,NOBC(ISYS)
      IX=IBC(1,I,ISYS)
      IY=IBC(2,I,ISYS)
      IF(FSM(IX,IY).EQ.-1 .OR. FSM(IX,IY).EQ.-2)  GO TO 215
       WRITE(OUT,1950)  ISYS,I,IX,IY,IBC(3,I,ISYS)
       CALL EXIT
  215 IF(NOBCALL.EQ.0) THEN
         NOBCALL = 1
         IBCALL(1,1) = IBC(1,I,ISYS)
         IBCALL(2,1) = IBC(2,I,ISYS)
         IBCALL(3,1) = IBC(3,I,ISYS)
      ELSE
         DO 220 J=1,NOBCALL
         IF(IBC(1,I,ISYS).EQ.IBCALL(1,J) .AND.
     .      IBC(2,I,ISYS).EQ.IBCALL(2,J) .AND.
     .          IBC(3,I,ISYS).EQ.IBCALL(3,J))  GO TO 222
  220    CONTINUE
         NOBCALL = NOBCALL+1
         IBCALL(1,NOBCALL) = IBC(1,I,ISYS)
         IBCALL(2,NOBCALL) = IBC(2,I,ISYS)
         IBCALL(3,NOBCALL) = IBC(3,I,ISYS)
      ENDIF
  222 BBC(I,ISYS) = SCALBC(ISYS)*BBC(I,ISYS) 
  225 CONTINUE

      GO TO 500

  300 IF(IBNRYRDOPTS(1).EQ.0)  THEN
C        READ STANDARD LEVEL DEPTHS
        READ(32,1000)  COMMENT
        READ(32,1020,ERR=950)  NLVLS
      ELSE
        READ(32)  NLVLS
      ENDIF
      IF(NLVLS.GT.NSL)  GO TO 970
      IF(IBNRYRDOPTS(1).EQ.0)  THEN
        READ(32,1000)  COMMENT
        READ(32,1100,ERR=950)  (SLDEPTH(I,ISYS),I=1,NLVLS)
      ELSE
        READ(32)  (SLDEPTH(I,ISYS),I=1,NLVLS)
      ENDIF
      IF(LIST(2).EQ.1)   WRITE(OUT,3000)  ISYS,
     .          (I,I=1,10),(SLDEPTH(I,ISYS),I=1,NLVLS)
 3000 FORMAT(42X,'STANDARD LEVEL DEPTHS FOR SYSTEM',I3/
     .    17X,10(I2,8X)/(11X,10F10.2))
      DO 305 I=1,NLVLS
        SLDEPTH(I,ISYS) = -SLDEPTH(I,ISYS)
  305 CONTINUE

      IF(IBCOPT.EQ.4)  GO TO 400

C        CONSTANT - STANDARD LEVELS
      IF(LIST(2).EQ.1)  WRITE(OUT,3100)   (I,I=1,10)
 3100 FORMAT(/43X,'CONCENTRATIONS AT STANDARD LEVELS'/
     .   2X,'ROW COL',8X,10(I2,8X))
      JBC = 0
      DO 325 I=1,NOBC(ISYS)
        IF(IBNRYRDOPTS(1).EQ.0)  THEN
          READ(32,1000)  COMMENT
          READ(32,1020,ERR=950)  IX,IY,NOBCSL(I,ISYS)
          READ(32,1000)  COMMENT
          READ(32,1100,ERR=950)  (SLBC(J),J=1,NOBCSL(I,ISYS))
        ELSE
          READ(32)  IX,IY,NOBCSL(I,ISYS)
          READ(32)  (SLBC(J),J=1,NOBCSL(I,ISYS))
        ENDIF
        IF(LIST(2).EQ.1)   WRITE(OUT,3200)  IX,IY,
     .          (SLBC(J),J=1,NOBCSL(I,ISYS))
 3200   FORMAT(1X,2I4,2X,10E10.3)

        DO 310 IZ=1,NZ
          SIGDEPTH(IZ) = -ZZ(IZ)*HBAR(IX,IY)
  310   CONTINUE

        CALL SINTER(SLDEPTH(1,ISYS),SLBC,SIGDEPTH,SIGBC
     .     ,NOBCSL(I,ISYS),NZ)
        DO 320 IZ=1,NZ
          JBC = JBC + 1
          IF(JBC.GT.MXBC)  GO TO 960
          IBC(1,JBC,ISYS) = IX
          IBC(2,JBC,ISYS) = IY
          IBC(3,JBC,ISYS) = IZ
          BBC(JBC,ISYS) = SIGBC(IZ)
  320   CONTINUE
  325 CONTINUE

      NOBC(ISYS) = JBC

      IF(LIST(2).EQ.1)  WRITE(OUT,3300) 
 3300 FORMAT(//30X,'BOUNDARY CONCENTRATIONS AFTER TRANSFORMATION TO SIGM
     .A-LEVELS'/)
      IF(LIST(2).EQ.1)  WRITE(OUT,1700) 
      IF(LIST(2).EQ.1)  WRITE(OUT,1900)  (BBC(I,ISYS),(IBC(J,I,ISYS),
     .   J=1,3),I=1,NOBC(ISYS))

      DO 335 I=1,NOBC(ISYS)
      IF(NOBCALL.EQ.0) THEN
         NOBCALL = 1
         IBCALL(1,1) = IBC(1,I,ISYS)
         IBCALL(2,1) = IBC(2,I,ISYS)
         IBCALL(3,1) = IBC(3,I,ISYS)
      ELSE
         DO 330 J=1,NOBCALL
         IF(IBC(1,I,ISYS).EQ.IBCALL(1,J) .AND.
     .        IBC(2,I,ISYS).EQ.IBCALL(2,J) .AND.
     .          IBC(3,I,ISYS).EQ.IBCALL(3,J))  GO TO 333
  330    CONTINUE
         NOBCALL = NOBCALL+1
         IF(NOBCALL.GT.MXBC)  GO TO 960
         IBCALL(1,NOBCALL) = IBC(1,I,ISYS)
         IBCALL(2,NOBCALL) = IBC(2,I,ISYS)
         IBCALL(3,NOBCALL) = IBC(3,I,ISYS)
      ENDIF
C          APPLY SCALE FACTOR
  333 BBC(I,ISYS) = SCALBC(ISYS)*BBC(I,ISYS)
  335 CONTINUE

      GO TO 500

C        TIME-VARIABLE - STANDARD LEVELS
  400 IF(IBNRYRDOPTS(1).EQ.0)  THEN
        READ(32,1000)  COMMENT
        READ(32,2100,ERR=950)  (IBC(1,I,ISYS),IBC(2,I,ISYS),
     .     NOBCSL(I,ISYS),I=1,NOBC(ISYS))
      ELSE
        READ(32)  (IBC(1,I,ISYS),IBC(2,I,ISYS),NOBCSL(I,ISYS),
     .     I=1,NOBC(ISYS))
      ENDIF

      IF(LIST(2).EQ.1)  WRITE(OUT,3100)   (I,I=1,10)
      JBC = 0
      DO 415 I=1,NOBC(ISYS)
        IF(IBNRYRDOPTS(1).EQ.0)  THEN
          READ(32,1100,ERR=950)  (SLBC(J),J=1,NOBCSL(I,ISYS))
        ELSE
          READ(32)  (SLBC(J),J=1,NOBCSL(I,ISYS))
        ENDIF
        IF(LIST(2).EQ.1)   WRITE(OUT,3200)  IBC(1,I,ISYS),IBC(2,I,ISYS),
     .          (SLBC(J),J=1,NOBCSL(I,ISYS))

        DO 405 IZ=1,NZ
          SIGDEPTH(IZ) = -ZZ(IZ)*HBAR(IBC(1,I,ISYS),IBC(2,I,ISYS))
  405   CONTINUE

        CALL SINTER(SLDEPTH(1,ISYS),SLBC,SIGDEPTH,SIGBC
     .     ,NOBCSL(I,ISYS),NZ)
        DO 410 IZ=1,NZ
          JBC = JBC + 1
          IF(JBC.GT.MXBC)  GO TO 960
          BBC(JBC,ISYS) = SIGBC(IZ)
  410   CONTINUE
  415 CONTINUE

      DO 430 I=1,NOBC(ISYS)
        IFROM = (NOBC(ISYS)+1)-I
        ITO = NZ*(NOBC(ISYS)-I)
        DO 425 IZ=1,NZ
          ITO = ITO+1
          IBC(1,ITO,ISYS) = IBC(1,IFROM,ISYS)
          IBC(2,ITO,ISYS) = IBC(2,IFROM,ISYS)
          IBC(3,ITO,ISYS) = IZ
  425   CONTINUE
  430 CONTINUE

      NOBC(ISYS) = NZ*NOBC(ISYS)
      IF(LIST(2).EQ.1)  WRITE(OUT,3300) 
      IF(LIST(2).EQ.1)  WRITE(OUT,1700) 
      IF(LIST(2).EQ.1)  WRITE(OUT,1900)  (BBC(I,ISYS),(IBC(J,I,ISYS),
     .   J=1,3),I=1,NOBC(ISYS))

      DO 445 I=1,NOBC(ISYS)
      IF(NOBCALL.EQ.0) THEN
         NOBCALL = 1
         IBCALL(1,1) = IBC(1,I,ISYS)
         IBCALL(2,1) = IBC(2,I,ISYS)
         IBCALL(3,1) = IBC(3,I,ISYS)
      ELSE
         DO 440 J=1,NOBCALL
         IF(IBC(1,I,ISYS).EQ.IBCALL(1,J) .AND.
     .      IBC(2,I,ISYS).EQ.IBCALL(2,J) .AND.
     .          IBC(3,I,ISYS).EQ.IBCALL(3,J))  GO TO 443
  440    CONTINUE
         NOBCALL = NOBCALL+1
         IBCALL(1,NOBCALL) = IBC(1,I,ISYS)
         IBCALL(2,NOBCALL) = IBC(2,I,ISYS)
         IBCALL(3,NOBCALL) = IBC(3,I,ISYS)
      ENDIF
  443 BBC(I,ISYS) = SCALBC(ISYS)*BBC(I,ISYS) 
  445 CONTINUE

  500 CONTINUE

!
!-----------------------------------------------------------------------
! Optionally read next boundary time 
!-----------------------------------------------------------------------
!
      IF(IBCOPT.EQ.1 .OR. IBCOPT.EQ.3)  RETURN
      IF(INPCHCK.EQ.0)  THEN
       IF(IBCOPT.EQ.2 .OR. IBCOPT.EQ.4)  THEN
        IF(IBCPWLOPT.EQ.0)  THEN
          IF(IBNRYRDOPTS(1).EQ.0)  THEN
            READ(32,1110,ERR=950)  NXBCTSECS
          ELSEIF(IBNRYRDOPTS(1).EQ.1)  THEN
            READ(32)  NXBCTSECS
          ELSEIF(IBNRYRDOPTS(1).EQ.2)  THEN
            bryTflag=bryTflag+1
            CALL ncbry_time(ncID_bry)
          ENDIF
          NXBCTSECS=ISCALBC*NXBCTSECS
          NXBCT=NXBCTSECS/86400.
        ELSE
          CALL RCA11
        ENDIF
       ENDIF

      ELSE

  550  IF(IBCPWLOPT.EQ.0)  THEN
         IF(IBNRYRDOPTS(1).EQ.0)  THEN
           READ(32,1110,ERR=950)  NXBCTSECS
         ELSE
           READ(32)  NXBCTSECS
         ENDIF
         NXBCTSECS=ISCALBC*NXBCTSECS
         NXBCT=NXBCTSECS/86400.
       ENDIF
       CALL RCA11
       IF(INPCHCK.EQ.1)  GO TO 550
       INPCHCK=1
      ENDIF

  930 RETURN

  950 IN=32
      CALL FMTER
      CALL EXIT 
  960 WRITE(OUT,9600)  
 9600 FORMAT(//11X,'ERROR...USER REQUESTING MORE BOUNDARY CONDITIONS THA
     .N DIMENSIONED FOR IN -RCACM-'/11X,'RCA TERMINATED'//)
      CALL EXIT 
  970 WRITE(OUT,9700)  
 9700 FORMAT(//11X,'ERROR...USER REQUESTING MORE STANDARD LEVELS THAN DI
     .MENSIONED FOR IN -RCACM-'/11X,'RCA TERMINATED'//)
      CALL EXIT 
      END 
