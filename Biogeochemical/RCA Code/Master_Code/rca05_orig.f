      SUBROUTINE RCA05
C 
C        RCA05 READS THE FORCING FUNCTIONS
C 
      SAVE
      INCLUDE  'RCACM'
      CHARACTER  LOCATION(MXWK)*50

C        NOTE: SINCE LOADS WILL BE DIVIDED BY THE SEGMENT VOLUMES IN
C        FORMULATING THE MASS BALANCE DERIVATIVE EQUATIONS THE FOLLOWING
C        CONVERSION IS NECESSARY 
C           KG/DAY * 10^6 MG/KG * M^3/10^3 L ==> MG*M^3/DAY-L

      WRITE(OUT,8000) 
 8000 FORMAT(////1X,119('*')//)
      TWARPPS='SECS'
      TWARPNPS='SECS'
      TWARPFL='SECS'
      TWARPATM='SECS'

C             POINT SOURCE LOADS

      READ(IN,1000)  COMMENT
 1000 FORMAT(A)
      WRITE(OUT,1001)
 1001 FORMAT(//43X,'---- POINT SOURCE LOADINGS ----')
      READ(IN,1005)  PSFILNA,IBNRYRDOPTS(2)
 1005 FORMAT(A40,I10)
      IF(PSFILNA.EQ.'NULL' .OR. PSFILNA.EQ.'null' .OR. 
     .                          PSFILNA.EQ.'Null')  THEN
        PSFILNA='NULL'
        WRITE(OUT,1006)
 1006   FORMAT(38X,'USER DID NOT SPECIFY A PS FILE TO BE READ')
        DO ISYS=1,NOSYS 
          NOPS(ISYS)=0
        ENDDO
        GO TO 100
      ENDIF

      WRITE(OUT,1010)  PSFILNA
 1010 FORMAT(//40X,'OPENING PSFILNA = ',A40)
      IF(IBNRYRDOPTS(2).EQ.0)  THEN
        OPEN(33,FILE=PSFILNA,FORM='FORMATTED')
        READ(33,1000) COMMENT
        READ(33,1050,ERR=933) IPSOPT,IPSPWLOPT
 1050   FORMAT(3I10) 
        READ(33,1000) COMMENT
        WRITE(OUT,1020)
 1020   FORMAT(//41X,'POINT SOURCE LOADING TABLE'/3X,'LOAD'/
     .   2X,'NUMBER   IX   IY',1X,20('-'),' LOCATION ',20('-'),1X,
     .   '---- VERTICAL DISTRIBUTION OF LOADS (FRACTION/LAYER) ---')
        DO 20 IWK=1,MXWK
         READ(33,1030,ERR=933)  (IPSTABL(I,IWK),I=1,3)
     .                         ,LOCATION(IWK)
 1030    FORMAT(3I5,A50) 
         IF(IPSTABL(1,IWK).LT.0)  GO TO 25
         READ(33,1035,ERR=933)  (ZFRACPS(IZ,IWK),IZ=1,NZ)
 1035    FORMAT(10X,10F6.3)
         WRITE(OUT,1040)  (IPSTABL(I,IWK),I=1,3)
     .                   ,LOCATION(IWK),(ZFRACPS(IZ,IWK),IZ=1,NZ)
 1040    FORMAT(I7,I6,I5,A50,1X,10F5.2/(69X,10F5.2))
   20   CONTINUE
   25   IF(IPSOPT.EQ.2) THEN
          READ(33,1000) COMMENT
          READ(33,1100,ERR=933)  NXPSTSECS,TWARPPS
 1100     FORMAT(10X,I10,6X,A4)
        ENDIF
      ELSE
        OPEN(33,FILE=PSFILNA,FORM='UNFORMATTED')
        READ(33) IPSOPT,IPSPWLOPT
        WRITE(OUT,1020)
        DO 30 IWK=1,MXWK
         READ(33,ERR=943)  (IPSTABL(I,IWK),I=1,3),LOCATION(IWK)
         IF(IPSTABL(1,IWK).LT.0)  GO TO 32
         READ(33,ERR=943)  (ZFRACPS(IZ,IWK),IZ=1,NZ)
         WRITE(OUT,1040)  (IPSTABL(I,IWK),I=1,3)
     .                   ,LOCATION(IWK),(ZFRACPS(IZ,IWK),IZ=1,NZ)
   30   CONTINUE
   32   IF(IPSOPT.EQ.2)  READ(33)  NXPSTSECS,TWARPPS
      ENDIF
      IF(IPSOPT.EQ.2) THEN
       ISCALPS=IUNITCHECK(TWARPPS,'TWARPPS ')
       NXPSTSECS=ISCALPS*NXPSTSECS
       NXPST=NXPSTSECS/86400.
      ENDIF
      WRITE(OUT,1120)  NXPST
 1120 FORMAT(//40X,'POINT SOURCE LOADS AT TIME =',F8.2,' DAYS')

      DO 80 ISYS=1,NOSYS 
      IF(IBNRYRDOPTS(2).EQ.0)  THEN
        READ(33,1000) COMMENT
        READ(33,1050,ERR=933) NOPS(ISYS)
      ELSE
        READ(33) NOPS(ISYS)
      ENDIF
      IF(NOPS(ISYS).GT.MXWK)  GO TO 960
      IF(NOPS(ISYS).EQ.0)   GO TO 70
      WRITE(OUT,2010) ISYS
 2010 FORMAT(//42X,'POINT SOURCE LOADINGS FOR SYSTEM',I3) 
      WRITE(OUT,1150)   SYNAME(ISYS)
 1150 FORMAT(41X,11('-'),3X,A8,3X,11('-'))
      WRITE(OUT,2020) IPSOPT,NOPS(ISYS)
 2020 FORMAT(39X,'PS OPTION',I3,' USED',5X,'NO.OF WK''S READ',I5)
      IF(IPSPWLOPT.EQ.0)  WRITE(OUT,2030)
 2030 FORMAT(45X,'STEP FUNCTION OPTION SELECTED')
      IF(IPSPWLOPT.EQ.1)  WRITE(OUT,2040)
 2040 FORMAT(43X,'PIECEWISE LINEAR OPTION SELECTED')
      IF(IBNRYRDOPTS(2).EQ.0)  THEN
        READ(33,1000) COMMENT
        READ(33,1500,ERR=933)  SCALPS(ISYS)
 1500   FORMAT(E10.3) 
      ELSE
        READ(33)  SCALPS(ISYS)
      ENDIF
      WRITE(OUT,1600)  SCALPS(ISYS)
 1600 FORMAT(47X,'SCALE FACTOR',E13.3/) 

      IF(IPSOPT.EQ.1)  THEN

        IF(LIST(3).EQ.1)  WRITE(OUT,2050) 
 2050   FORMAT(30X,'LOAD     IX  IY      LOCATION '/) 
        IF(IBNRYRDOPTS(2).EQ.0)  THEN
          READ(33,1000) COMMENT
          READ(33,1110,ERR=933)  (BPS(I,ISYS),IPS(4,I,ISYS)
     .        ,I=1,NOPS(ISYS)) 
 1110   FORMAT(5(F10.0,I5)) 
        ELSE
          READ(33)  (BPS(I,ISYS),IPS(4,I,ISYS),I=1,NOPS(ISYS))
        ENDIF
        DO 35 I=1,NOPS(ISYS)
         DO IWK=1,MXWK
          IF(IPS(4,I,ISYS).EQ.IPSTABL(1,IWK))  THEN
            IPS(1,I,ISYS)=IPSTABL(2,IWK)
            IPS(2,I,ISYS)=IPSTABL(3,IWK)
            IPS(3,I,ISYS)=IWK
            GO TO 35
          ENDIF
         ENDDO
C        ERROR --- SPECIFIED LOAD NOT FOUND IN LOADING TABLE --- EXIT
         GO TO 965
   35   CONTINUE
        IF(LIST(3).EQ.1)  WRITE(OUT,1700)   (BPS(I,ISYS),
     .     (IPS(J,I,ISYS),J=1,2),LOCATION(IPS(3,I,ISYS)),I=1,NOPS(ISYS))
 1700   FORMAT(25X,E12.4,2I4,2X,A50)
        TOTLD = 0.
        DO 40 I=1,NOPS(ISYS)
          BPS(I,ISYS) = 1000.*SCALPS(ISYS)*BPS(I,ISYS)
          TOTLD = TOTLD+BPS(I,ISYS)
   40   CONTINUE
        WRITE(OUT,2060)  TOTLD/1000.
 2060   FORMAT(/29X,'TOTAL POINT SOURCE LOADING =',E13.4,' KG/DAY')

      ELSE

        IF(LIST(3).EQ.1)   WRITE(OUT,2050)
        IF(IBNRYRDOPTS(2).EQ.0)  THEN
          READ(33,2700,ERR=933)  (IPS(4,I,ISYS),I=1,NOPS(ISYS))
 2700     FORMAT(10X,7I10)
          READ(33,2750,ERR=933)  (BPS(I,ISYS),I=1,NOPS(ISYS))
 2750     FORMAT(10X,7F10.0)
        ELSE
          READ(33)  (IPS(4,I,ISYS),I=1,NOPS(ISYS))
          READ(33)  (BPS(I,ISYS),I=1,NOPS(ISYS))
        ENDIF
        DO 45 I=1,NOPS(ISYS)
         DO IWK=1,MXWK
          IF(IPS(4,I,ISYS).EQ.IPSTABL(1,IWK))  THEN
            IPS(1,I,ISYS)=IPSTABL(2,IWK)
            IPS(2,I,ISYS)=IPSTABL(3,IWK)
            IPS(3,I,ISYS)=IWK
            GO TO 45
          ENDIF
         ENDDO
C        ERROR --- SPECIFIED LOAD NOT FOUND IN LOADING TABLE --- EXIT
         GO TO 965
   45   CONTINUE
        IF(LIST(3).EQ.1)  WRITE(OUT,1700)   (BPS(I,ISYS),
     .     (IPS(J,I,ISYS),J=1,2),LOCATION(IPS(3,I,ISYS)),I=1,NOPS(ISYS))
        TOTLD = 0.0
        DO 50 I=1,NOPS(ISYS)
          BPS(I,ISYS) = 1000.*SCALPS(ISYS)*BPS(I,ISYS) 
          TOTLD = TOTLD+BPS(I,ISYS)
   50   CONTINUE
        WRITE(OUT,2080)  TIME,TOTLD/1000.
 2080   FORMAT(/29X,'TOTAL POINT SOURCE LOADING AT TIME,',F7.2,
     .    ' =',E13.4,' KG/DAY')

      ENDIF

C        CHECK FOR ERRONEOUS SEGMENT NUMBER
      DO 60 I=1,NOPS(ISYS)
       IX=IPS(1,I,ISYS)
       IY=IPS(2,I,ISYS)
       IF(IX.LE.0 .OR. IX.GT.NX .OR.
     .    IY.LE.0 .OR. IY.GT.NY)  GO TO 980
       IF(FSM(IX,IY).NE.1.)  GO TO 990
   60 CONTINUE

      GO TO 80 

   70 WRITE(OUT,2090)    ISYS 
 2090 FORMAT(///40X,'NO POINT SOURCE LOADINGS FOR SYSTEM',I3) 
      WRITE(OUT,1150)   SYNAME(ISYS)

   80 CONTINUE

      IF(IPSOPT.EQ.1)  GO TO 100
      IF(INPCHCK.EQ.0)  THEN
       IF(IPSOPT.EQ.2)  THEN
        IF(IPSPWLOPT.EQ.0) THEN
         IF(IBNRYRDOPTS(2).EQ.0)  THEN
           READ(33,1100,ERR=933)  NXPSTSECS
         ELSE
           READ(33)  NXPSTSECS
         ENDIF
         NXPSTSECS=ISCALPS*NXPSTSECS
         NXPST=NXPSTSECS/86400.
        ELSE
         CALL RCA10(SPS,BPS,MXWK,NOPS,NXPST,33,IPSPWLOPT,SCALPS)
        ENDIF
       ENDIF

      ELSE

   90  IF(IPSPWLOPT.EQ.0)  THEN
         IF(IBNRYRDOPTS(2).EQ.0)  THEN
           READ(33,1100,ERR=933)  NXPSTSECS
         ELSE
           READ(33)   NXPSTSECS
         ENDIF
         NXPSTSECS=ISCALPS*NXPSTSECS
         NXPST=NXPSTSECS/86400.
       ENDIF
       CALL RCA10(SPS,BPS,MXWK,NOPS,NXPST,33,IPSPWLOPT,SCALPS)
       BACKSPACE(33)
       IF(INPCHCK.EQ.1)  GO TO  90
       INPCHCK=1

      ENDIF

  100 CONTINUE

C             NONPOINT SOURCE LOADS

      READ(IN,1000)  COMMENT
      WRITE(OUT,1002)
 1002 FORMAT(//42X,'---- NONPOINT SOURCE LOADINGS ----')
      READ(IN,1005)  NPSFILNA,IBNRYRDOPTS(3)
      IF(NPSFILNA.EQ.'NULL' .OR. NPSFILNA.EQ.'null' .OR. 
     .                          NPSFILNA.EQ.'Null')  THEN
        NPSFILNA='NULL'
        WRITE(OUT,3006)
 3006   FORMAT(38X,'USER DID NOT SPECIFY A NPS FILE TO BE READ')
        DO ISYS=1,NOSYS 
          NONPS(ISYS)=0
        ENDDO
        GO TO 200
      ENDIF
      WRITE(OUT,3010)  NPSFILNA
 3010 FORMAT(//40X,'OPENING NPSFILNA = ',A40)
      IF(IBNRYRDOPTS(3).EQ.0)  THEN
        OPEN(34,FILE=NPSFILNA,FORM='FORMATTED')
        READ(34,1000) COMMENT
        READ(34,1050,ERR=934) INPSOPT,INPSPWLOPT
        READ(34,1000) COMMENT
        WRITE(OUT,3020)
 3020   FORMAT(//39X,'NONPOINT SOURCE LOADING TABLE'/3X,'LOAD'/
     .   2X,'NUMBER   IX   IY',1X,20('-'),' LOCATION ',20('-'),1X,
     .   '---- VERTICAL DISTRIBUTION OF LOADS (FRACTION/LAYER) ---')
        DO 120 IWK=1,MXWK
         READ(34,1030,ERR=934)  (INPSTABL(I,IWK),I=1,3)
     .                         ,LOCATION(IWK)
         IF(INPSTABL(1,IWK).LT.0)  GO TO 125
         READ(34,1035,ERR=934)  (ZFRACNPS(IZ,IWK),IZ=1,NZ)
         WRITE(OUT,1040)  (INPSTABL(I,IWK),I=1,3)
     .                   ,LOCATION(IWK),(ZFRACNPS(IZ,IWK),IZ=1,NZ)
  120   CONTINUE
  125   IF(INPSOPT.EQ.2)  THEN
          READ(34,1000) COMMENT
          READ(34,1100,ERR=934)  NXNPSTSECS,TWARPNPS
        ENDIF
      ELSE
        OPEN(34,FILE=NPSFILNA,FORM='UNFORMATTED')
        READ(34) INPSOPT,INPSPWLOPT
        WRITE(OUT,3020)
        DO 130 IWK=1,MXWK
         READ(34,ERR=944)  (INPSTABL(I,IWK),I=1,3),LOCATION(IWK)
         IF(INPSTABL(1,IWK).LT.0)  GO TO 132
         READ(34,ERR=944)  (ZFRACNPS(IZ,IWK),IZ=1,NZ)
         WRITE(OUT,1040)  (INPSTABL(I,IWK),I=1,3)
     .                   ,LOCATION(IWK),(ZFRACNPS(IZ,IWK),IZ=1,NZ)
  130   CONTINUE
  132   IF(INPSOPT.EQ.2)  READ(34)  NXNPSTSECS,TWARPNPS
      ENDIF
      IF(INPSOPT.EQ.2) THEN
       ISCALNPS=IUNITCHECK(TWARPNPS,'TWARPNPS')
       NXNPSTSECS=ISCALNPS*NXNPSTSECS
       NXNPST=NXNPSTSECS/86400.
      ENDIF
      WRITE(OUT,1121)  NXNPST
 1121 FORMAT(//38X,'NONPOINT SOURCE LOADS AT TIME =',F8.2,' DAYS')

      DO 180 ISYS=1,NOSYS 
      IF(IBNRYRDOPTS(3).EQ.0)  THEN
        READ(34,1000) COMMENT
        READ(34,1050,ERR=934) NONPS(ISYS)
      ELSE
        READ(34) NONPS(ISYS)
      ENDIF
      IF(NONPS(ISYS).GT.MXWK)  GO TO 962
      IF(NONPS(ISYS).EQ.0)   GO TO 170
      WRITE(OUT,2011) ISYS
 2011 FORMAT(//41X,'NONPOINT SOURCE LOADINGS FOR SYSTEM',I3) 
      WRITE(OUT,1150)   SYNAME(ISYS)
      WRITE(OUT,2021) INPSOPT,NONPS(ISYS)
 2021 FORMAT(39X,'NPS OPTION',I3,' USED',5X,'NO.OF WKS''S READ',I5)
      IF(INPSPWLOPT.EQ.0)  WRITE(OUT,2030)
      IF(INPSPWLOPT.EQ.1)  WRITE(OUT,2040)
      IF(IBNRYRDOPTS(3).EQ.0)  THEN
        READ(34,1000) COMMENT
        READ(34,1500,ERR=934)  SCALNPS(ISYS)
      ELSE
        READ(34)  SCALNPS(ISYS)
      ENDIF
      WRITE(OUT,1600)  SCALNPS(ISYS)

      IF(INPSOPT.EQ.1)  THEN
        IF(LIST(3).EQ.1)  WRITE(OUT,2050) 
        IF(IBNRYRDOPTS(3).EQ.0)  THEN
          READ(34,1000) COMMENT
          READ(34,1110,ERR=934)  (BNPS(I,ISYS),INPS(4,I,ISYS),
     .        I=1,NONPS(ISYS)) 
        ELSE
          READ(34)  (BNPS(I,ISYS),INPS(J,I,ISYS),
     .        I=1,NONPS(ISYS))
        ENDIF
        DO 135 I=1,NONPS(ISYS)
         DO IWK=1,MXWK
          IF(INPS(4,I,ISYS).EQ.INPSTABL(1,IWK))  THEN
            INPS(1,I,ISYS)=INPSTABL(2,IWK)
            INPS(2,I,ISYS)=INPSTABL(3,IWK)
            INPS(3,I,ISYS)=IWK
            GO TO 135
          ENDIF
         ENDDO
C        ERROR --- SPECIFIED LOAD NOT FOUND IN LOADING TABLE --- EXIT
         GO TO 970
  135   CONTINUE
        IF(LIST(3).EQ.1)  WRITE(OUT,1700)   (BNPS(I,ISYS),
     .     (INPS(J,I,ISYS),J=1,2),LOCATION(INPS(3,I,ISYS)),
     .                                          I=1,NONPS(ISYS))
C        MULTIPLY BY SCALE FACTOR AND APPLY CONVERSION FACTOR
        TOTLD = 0.0
        DO 140 I=1,NONPS(ISYS)
          BNPS(I,ISYS) = 1000.*SCALNPS(ISYS)*BNPS(I,ISYS)
          TOTLD = TOTLD+BNPS(I,ISYS)
  140   CONTINUE
        WRITE(OUT,3060)  TOTLD/1000.
 3060   FORMAT(/29X,'TOTAL NONPOINT SOURCE LOADING =',E13.4,' KG/DAY')

      ELSE

        IF(LIST(3).EQ.1)   WRITE(OUT,2050)
        IF(IBNRYRDOPTS(3).EQ.0)  THEN
          READ(34,2700,ERR=934)  (INPS(4,I,ISYS),I=1,NONPS(ISYS))
          READ(34,2750,ERR=934)  (BNPS(I,ISYS),I=1,NONPS(ISYS))
        ELSE
          READ(34)  (INPS(4,I,ISYS),I=1,NONPS(ISYS))
          READ(34)  (BNPS(I,ISYS),I=1,NONPS(ISYS))
        ENDIF
        DO 145 I=1,NONPS(ISYS)
         DO IWK=1,MXWK
          IF(INPS(4,I,ISYS).EQ.INPSTABL(1,IWK))  THEN
            INPS(1,I,ISYS)=INPSTABL(2,IWK)
            INPS(2,I,ISYS)=INPSTABL(3,IWK)
            INPS(3,I,ISYS)=IWK
            GO TO 145
          ENDIF
         ENDDO
C        ERROR --- SPECIFIED LOAD NOT FOUND IN LOADING TABLE --- EXIT
         GO TO 970
  145   CONTINUE
        IF(LIST(3).EQ.1)  WRITE(OUT,1700)   (BNPS(I,ISYS),
     .     (INPS(J,I,ISYS),J=1,2),LOCATION(INPS(3,I,ISYS)),
     .                                          I=1,NONPS(ISYS))
        TOTLD = 0.0
        DO 150 I=1,NONPS(ISYS)
          BNPS(I,ISYS) = 1000.*SCALNPS(ISYS)*BNPS(I,ISYS) 
          TOTLD = TOTLD+BNPS(I,ISYS)
  150   CONTINUE
        WRITE(OUT,3080)  TIME,TOTLD/1000.
 3080   FORMAT(/29X,'TOTAL NONPOINT SOURCE LOADING AT TIME,',F7.2,
     .    ' =',E13.4,' KG/DAY')

      ENDIF

C        CHECK FOR ERRONEOUS SEGMENT NUMBER
      DO 160 I=1,NONPS(ISYS)
       IX=INPS(1,I,ISYS)
       IY=INPS(2,I,ISYS)
       IF(IX.LE.0 .OR. IX.GT.NX .OR.
     .    IY.LE.0 .OR. IY.GT.NY)  GO TO 980
       IF(FSM(IX,IY).NE.1.)  GO TO 990
 160  CONTINUE

      GO TO 180 

  170 WRITE(OUT,3090)    ISYS 
 3090 FORMAT(///40X,'NO NONPOINT SOURCE LOADINGS FOR SYSTEM',I3) 
      WRITE(OUT,1150)   SYNAME(ISYS)

  180 CONTINUE

      IF(INPSOPT.EQ.1)  GO TO 200
      IF(INPCHCK.EQ.0)  THEN
       IF(INPSOPT.EQ.2)  THEN
        IF(INPSPWLOPT.EQ.0)  THEN
         IF(IBNRYRDOPTS(3).EQ.0)  THEN
           READ(34,1100,ERR=934)  NXNPSTSECS
         ELSE
           READ(34)  NXNPSTSECS
         ENDIF
         NXNPSTSECS=ISCALNPS*NXNPSTSECS
         NXNPST=NXNPSTSECS/86400.
        ELSE
         CALL RCA10(SNPS,BNPS,MXWK,NONPS,NXNPST,34,INPSPWLOPT,SCALNPS)
        ENDIF
       ENDIF

      ELSE

  190  IF(INPSPWLOPT.EQ.0)  THEN
         IF(IBNRYRDOPTS(3).EQ.0)  THEN
           READ(34,1100,ERR=934)  NXNPSTSECS
         ELSE
           READ(34)   NXNPSTSECS
         ENDIF
         NXNPSTSECS=ISCALNPS*NXNPSTSECS
         NXNPST=NXNPSTSECS/86400.
       ENDIF
       CALL RCA10(SNPS,BNPS,MXWK,NONPS,NXNPST,34,INPSPWLOPT,SCALNPS)
       BACKSPACE(34)
       IF(INPCHCK.EQ.1)  GO TO 190
       INPCHCK=1

      ENDIF

  200 CONTINUE

C             FALL-LINE LOADS

      READ(IN,1000)  COMMENT
      WRITE(OUT,1003)
 1003 FORMAT(//42X,'---- FALL-LINE SOURCE LOADINGS ----')
      READ(IN,1005)  FLFILNA,IBNRYRDOPTS(4)
      IF(FLFILNA.EQ.'NULL' .OR. FLFILNA.EQ.'null' .OR. 
     .                          FLFILNA.EQ.'Null')  THEN
        FLFILNA='NULL'
        WRITE(OUT,4006)
 4006   FORMAT(38X,'USER DID NOT SPECIFY A FL FILE TO BE READ')
        DO ISYS=1,NOSYS 
          NOFL(ISYS)=0
        ENDDO
        GO TO 300
      ENDIF
      WRITE(OUT,4010)  FLFILNA
 4010 FORMAT(//40X,'OPENING FLFILNA = ',A40)
      IF(IBNRYRDOPTS(4).EQ.0)  THEN
        OPEN(35,FILE=FLFILNA,FORM='FORMATTED')
        READ(35,1000) COMMENT
        READ(35,1050,ERR=935) IFLOPT,IFLPWLOPT
        READ(35,1000) COMMENT
        WRITE(OUT,4020)
 4020   FORMAT(//42X,'FALL LINE LOADING TABLE'/3X,'LOAD'/
     .   2X,'NUMBER   IX   IY',1X,20('-'),' LOCATION ',20('-'),1X,
     .   '---- VERTICAL DISTRIBUTION OF LOADS (FRACTION/LAYER) ---')
        DO 320 IWK=1,MXWK
        READ(35,1030,ERR=935)  (IFLTABL(I,IWK),I=1,3)
     .                        ,LOCATION(IWK)
        IF(IFLTABL(1,IWK).LT.0)  GO TO 325
         READ(35,1035,ERR=935)  (ZFRACFL(IZ,IWK),IZ=1,NZ)
         WRITE(OUT,1040)  (IFLTABL(I,IWK),I=1,3)
     .                   ,LOCATION(IWK),(ZFRACFL(IZ,IWK),IZ=1,NZ)
  320   CONTINUE
  325   IF(IFLOPT.EQ.2)  THEN
          READ(35,1000) COMMENT
          READ(35,1100,ERR=935)  NXFLTSECS,TWARPFL
        ENDIF
      ELSE
        OPEN(35,FILE=FLFILNA,FORM='UNFORMATTED')
        READ(35) IFLOPT,IFLPWLOPT
        WRITE(OUT,4020)
        DO 330 IWK=1,MXWK
        READ(35,ERR=945)  (IFLTABL(I,IWK),I=1,3),LOCATION(IWK)
        IF(IFLTABL(1,IWK).LT.0)  GO TO 332
         READ(35,ERR=945)  (ZFRACFL(IZ,IWK),IZ=1,NZ)
         WRITE(OUT,1040)  (IFLTABL(I,IWK),I=1,3)
     .                   ,LOCATION(IWK),(ZFRACFL(IZ,IWK),IZ=1,NZ)
  330   CONTINUE
  332   IF(IFLOPT.EQ.2)  READ(35)  NXFLTSECS,TWARPFL
      ENDIF
      IF(IFLOPT.EQ.2) THEN
       ISCALFL=IUNITCHECK(TWARPFL,'TWARPFL ')
       NXFLTSECS=ISCALFL*NXFLTSECS
       NXFLT=NXFLTSECS/86400.
      ENDIF
      WRITE(OUT,1122)  NXFLT
 1122 FORMAT(//41X,'FALL-LINE LOADS AT TIME =',F8.2,' DAYS')

      DO 280 ISYS=1,NOSYS 

      IF(IBNRYRDOPTS(4).EQ.0)  THEN
        READ(35,1000) COMMENT
        READ(35,1050,ERR=935) NOFL(ISYS)
      ELSE
        READ(35) NOFL(ISYS)
      ENDIF
      IF(NOFL(ISYS).GT.MXWK)  GO TO 964
      IF(NOFL(ISYS).EQ.0)   GO TO 270
      WRITE(OUT,4030) ISYS
 4030 FORMAT(//43X,'FALL-LINE LOADINGS FOR SYSTEM',I3) 
      WRITE(OUT,1150)   SYNAME(ISYS)
      WRITE(OUT,4040) IFLOPT,NOFL(ISYS)
 4040 FORMAT(39X,'FL OPTION',I3,' USED',5X,'NO.OF FL''S READ',I5)
      IF(IFLPWLOPT.EQ.0)  WRITE(OUT,2030)
      IF(IFLPWLOPT.EQ.1)  WRITE(OUT,2040)
      IF(IBNRYRDOPTS(4).EQ.0)  THEN
        READ(35,1000) COMMENT
        READ(35,1500,ERR=935)  SCALFL(ISYS)
      ELSE
        READ(35)  SCALFL(ISYS)
      ENDIF
      WRITE(OUT,1600)  SCALFL(ISYS)

      IF(IFLOPT.EQ.1)  THEN

        IF(LIST(3).EQ.1)  WRITE(OUT,4050) 
 4050   FORMAT(12X,4(8X,'FL  ROW COL LYR ')/) 
        IF(IBNRYRDOPTS(4).EQ.0)  THEN
          READ(35,1000) COMMENT
          READ(35,1110,ERR=935)  (BFL(I,ISYS),IFL(4,I,ISYS),
     .        I=1,NOFL(ISYS)) 
        ELSE
          READ(35)  (BFL(I,ISYS),IFL(4,I,ISYS),I=1,NOFL(ISYS)) 
        ENDIF
        DO 235 I=1,NOFL(ISYS)
         DO IWK=1,MXWK
          IF(IFL(4,I,ISYS).EQ.IFLTABL(1,IWK))  THEN
             IFL(1,I,ISYS)=IFLTABL(2,IWK)
             IFL(2,I,ISYS)=IFLTABL(3,IWK)
             IFL(3,I,ISYS)=IWK
             GO TO 235
          ENDIF
         ENDDO
C        ERROR --- SPECIFIED LOAD NOT FOUND IN LOADING TABLE --- EXIT
         GO TO 975
  235   CONTINUE
        IF(LIST(3).EQ.1)  WRITE(OUT,1700)   (BFL(I,ISYS),
     .     (IFL(J,I,ISYS),J=1,2),LOCATION(IFL(3,I,ISYS)),I=1,NOFL(ISYS))
C        MULTIPLY BY SCALE FACTOR AND APPLY CONVERSION FACTOR
        TOTLD = 0.0
        DO 240 I=1,NOFL(ISYS)
          BFL(I,ISYS) = 1000.*SCALFL(ISYS)*BFL(I,ISYS)
          TOTLD = TOTLD+BFL(I,ISYS)
  240  CONTINUE
        WRITE(OUT,4060)  TOTLD/1000.
 4060   FORMAT(/29X,'TOTAL FALL-LINE SOURCE LOADING =',E13.4,' KG/DAY')

      ELSE

        IF(LIST(3).EQ.1)   WRITE(OUT,2050)
        IF(IBNRYRDOPTS(4).EQ.0)  THEN
          READ(35,2700,ERR=935)  (IFL(4,I,ISYS),I=1,NOFL(ISYS))
          READ(35,2750,ERR=935)  (BFL(I,ISYS),I=1,NOFL(ISYS))
        ELSE
          READ(35)  (IFL(4,I,ISYS),I=1,NOFL(ISYS))
          READ(35)  (BFL(I,ISYS),I=1,NOFL(ISYS))
        ENDIF
        DO 245 I=1,NOFL(ISYS)
         DO IWK=1,MXWK
          IF(IFL(4,I,ISYS).EQ.IFLTABL(1,IWK))  THEN
             IFL(1,I,ISYS)=IFLTABL(2,IWK)
             IFL(2,I,ISYS)=IFLTABL(3,IWK)
             IFL(3,I,ISYS)=IWK
             GO TO 245
          ENDIF
         ENDDO
C        ERROR --- SPECIFIED LOAD NOT FOUND IN LOADING TABLE --- EXIT
         GO TO 975
  245   CONTINUE
        IF(LIST(3).EQ.1)  WRITE(OUT,1700)   (BFL(I,ISYS),
     .     (IFL(J,I,ISYS),J=1,2),LOCATION(IFL(3,I,ISYS)),I=1,NOFL(ISYS))
        TOTLD = 0.0
        DO 250 I=1,NOFL(ISYS)
          BFL(I,ISYS) = 1000.*SCALFL(ISYS)*BFL(I,ISYS) 
          TOTLD = TOTLD+BFL(I,ISYS)
  250   CONTINUE
        WRITE(OUT,4080)  TIME,TOTLD/1000.
 4080   FORMAT(/29X,'TOTAL FALL-LINE SOURCE LOADING AT TIME,',F7.2,
     .    ' =',E13.4,' KG/DAY')

      ENDIF

C        CHECK FOR ERRONEOUS SEGMENT NUMBER
      DO 260 I=1,NOFL(ISYS)
       IX=IFL(1,I,ISYS)
       IY=IFL(2,I,ISYS)
       IF(IX.LE.0 .OR. IX.GT.NX .OR.
     .    IY.LE.0 .OR. IY.GT.NY)  GO TO 980
       IF(FSM(IX,IY).NE.1.)  GO TO 990
 260  CONTINUE

      GO TO 280 

  270 WRITE(OUT,4090)    ISYS 
 4090 FORMAT(///40X,'NO FALL-LINE LOADINGS FOR SYSTEM',I3) 
      WRITE(OUT,1150)   SYNAME(ISYS)

  280 CONTINUE

      IF(IFLOPT.EQ.1)  GO TO 300
      IF(INPCHCK.EQ.0)  THEN
       IF(IFLOPT.EQ.2)  THEN
        IF(IFLPWLOPT.EQ.0)  THEN
         IF(IBNRYRDOPTS(4).EQ.0)  THEN
           READ(35,1100,ERR=935)  NXFLTSECS
         ELSE
           READ(35)  NXFLTSECS
         ENDIF
         NXFLTSECS=ISCALFL*NXFLTSECS
         NXFLT=NXFLTSECS/86400.
        ELSE
         CALL RCA10(SFL,BFL,MXWK,NOFL,NXFLT,35,IFLPWLOPT,SCALFL)
        ENDIF
       ENDIF

      ELSE

  290  IF(IFLPWLOPT.EQ.0)  THEN
         IF(IBNRYRDOPTS(4).EQ.0)  THEN
           READ(35,1100,ERR=935)  NXFLTSECS
         ELSE
           READ(35)   NXFLTSECS
         ENDIF
         NXFLTSECS=ISCALFL*NXFLTSECS
         NXFLT=NXFLTSECS/86400.
       ENDIF
       CALL RCA10(SFL,BFL,MXWK,NOFL,NXFLT,35,IFLPWLOPT,SCALFL)
       BACKSPACE(35)
       IF(INPCHCK.EQ.1)  GO TO 290
       INPCHCK=1

      ENDIF

  300 CONTINUE

C             ATMOSPHERIC LOADS

      READ(IN,1000)  COMMENT
      WRITE(OUT,1004)
 1004 FORMAT(//41X,'---- ATMOSPHERIC SOURCE LOADINGS ----')
      READ(IN,1005)  ATMFILNA,IBNRYRDOPTS(5)
      IF(ATMFILNA.EQ.'NULL' .OR. ATMFILNA.EQ.'null' .OR. 
     .                          ATMFILNA.EQ.'Null')  THEN
        ATMFILNA='NULL'
        WRITE(OUT,5010)
 5010   FORMAT(38X,'USER DID NOT SPECIFY AN ATM FILE TO BE READ')
        DO ISYS=1,NOSYS 
          NOATM(ISYS)=0
        ENDDO
        GO TO 400
      ENDIF
      WRITE(OUT,5020)  ATMFILNA
 5020 FORMAT(//40X,'OPENING ATMFILNA = ',A40)
      IF(IBNRYRDOPTS(5).EQ.0)  THEN
        OPEN(36,FILE=ATMFILNA,FORM='FORMATTED')
        READ(36,1000) COMMENT
        READ(36,1050,ERR=936) IATMOPT,IATMPWLOPT
        IF(IATMOPT.EQ.2)  THEN
          READ(36,1000) COMMENT
          READ(36,1100,ERR=936)  NXATMTSECS,TWARPATM
        ENDIF
      ELSE
        OPEN(36,FILE=ATMFILNA,FORM='UNFORMATTED')
        READ(36) IATMOPT,IATMPWLOPT
        IF(IATMOPT.EQ.2)  READ(36)  NXATMTSECS,TWARPATM
      ENDIF
      IF(IATMOPT.EQ.2) THEN
       ISCALATM=IUNITCHECK(TWARPATM,'TWARPATM')
       NXATMTSECS=ISCALATM*NXATMTSECS
       NXATMT=NXATMTSECS/86400.
      ENDIF
      WRITE(OUT,1123)  NXATMT
 1123 FORMAT(//40X,'ATMOSPHERIC LOADS AT TIME =',F8.2,' DAYS')

      DO 380 ISYS=1,NOSYS 

      IF(IBNRYRDOPTS(5).EQ.0)  THEN
        READ(36,1000) COMMENT
        READ(36,1050,ERR=936) NOATM(ISYS)
      ELSE
        READ(36) NOATM(ISYS)
      ENDIF
      IF(NOATM(ISYS).EQ.0)   GO TO 370
      WRITE(OUT,5030) ISYS
 5030 FORMAT(//43X,'ATMOSPHERIC LOADINGS FOR SYSTEM',I3) 
      WRITE(OUT,1150)   SYNAME(ISYS)
      WRITE(OUT,5040) IATMOPT,NOATM(ISYS)
 5040 FORMAT(39X,'IATMOPT OPTION',I3,' USED'
     .        5X,'NOATM OPTION',I5,' USED')
      IF(IATMPWLOPT.EQ.0)  WRITE(OUT,2030)
      IF(IATMPWLOPT.EQ.1)  WRITE(OUT,2040)
      IF(IBNRYRDOPTS(5).EQ.0)  THEN
        READ(36,1000) COMMENT
        READ(36,1500,ERR=936)  SCALATM(ISYS)
      ELSE
        READ(36)  SCALATM(ISYS)
      ENDIF
      WRITE(OUT,1600)  SCALATM(ISYS)

      IF(IATMOPT.EQ.1)  THEN
C        CONSTANT ATMOSPHERIC LOADINGS
        IF(IBNRYRDOPTS(5).EQ.0)  THEN
          IF(NOATM(ISYS).EQ.1)  THEN
            READ(36,1000) COMMENT
            READ(36,2750,ERR=936)  BATMC
            DO IY=1,NY
             DO IX=1,NX
              IF(FSM(IX,IY).EQ.1.)  BATM(IX,IY,ISYS)=BATMC
             ENDDO
            ENDDO
          ELSE
            READ(36,1000) COMMENT
            DO IY=1,NY
             READ(36,2750,ERR=936) (BATM(IX,IY,ISYS),IX=1,NX)
             DO IX=1,NX
              IF(FSM(IX,IY).NE.1.)  BATM(IX,IY,ISYS)=0.0
             ENDDO
            ENDDO
          ENDIF
        ELSE
          IF(NOATM(ISYS).EQ.1)  THEN
            READ(36)  BATMC
            DO IY=1,NY
             DO IX=1,NX
              IF(FSM(IX,IY).EQ.1.)  BATM(IX,IY,ISYS)=BATMC
             ENDDO
            ENDDO
          ELSE
            DO IY=1,NY
             READ(36) (BATM(IX,IY,ISYS),IX=1,NX)
             DO IX=1,NX
              IF(FSM(IX,IY).NE.1.)  BATM(IX,IY,ISYS)=0.0
             ENDDO
            ENDDO
          ENDIF
        ENDIF
        IF(LIST(3).EQ.1)  THEN
         IF(NOATM(ISYS).EQ.1) THEN
           WRITE(OUT,5050)  BATMC
 5050      FORMAT(42X,'ATMOSPHERIC LOADING =',E13.3)
         ELSE
           WRITE(OUT,5055)
 5055      FORMAT(42X,'SPATIALLY DEPENDENT ATMOSPHERIC LOADINGS'/
     .             3X,'IY')
           DO IY=1,NY
            WRITE(OUT,5060)  IY,(BATM(IX,IY,ISYS),IX=1,NX)
 5060       FORMAT(I5,10E11.3/(5X,10E11.3))
           ENDDO
         ENDIF
        ENDIF
C        MULTIPLY BY SCALE FACTOR AND APPLY CONVERSION FACTOR
C          KG/M^2-DAY ==> MG*M^3/DAY-M^2-L
        DO 335 IY=1,NY
         DO 335 IX=1,NX
           BATM(IX,IY,ISYS) = 1000.*SCALATM(ISYS)*BATM(IX,IY,ISYS)
  335   CONTINUE

        TOTLD = 0.0
        DO 340 IY=1,NY
         DO 340 IX=1,NX
           TOTLD = TOTLD + XAZ(IX,IY)*BATM(IX,IY,ISYS)
  340   CONTINUE
        WRITE(OUT,5065)  TOTLD/1000.
 5065   FORMAT(/36X,'TOTAL ATMOSPHERIC LOADING =',E13.3,' KG/DAY')

      ELSE

C        TIME-VARIABLE ATMOSPHERIC LOADINGS
        IF(IBNRYRDOPTS(5).EQ.0)  THEN
          IF(NOATM(ISYS).EQ.1)  THEN
            READ(36,1000)  COMMENT
            READ(36,2750,ERR=936)  BATMC
            DO IY=1,NY
             DO IX=1,NX
              IF(FSM(IX,IY).EQ.1.)  BATM(IX,IY,ISYS)=BATMC
             ENDDO
            ENDDO
          ELSE
            READ(36,1000)  COMMENT
            DO IY=1,NY
             READ(36,2750,ERR=936) (BATM(IX,IY,ISYS),IX=1,NX)
             DO IX=1,NX
              IF(FSM(IX,IY).NE.1.)  BATM(IX,IY,ISYS)=0.0
             ENDDO
            ENDDO
          ENDIF
        ELSE
          IF(NOATM(ISYS).EQ.1)  THEN
            READ(36)  BATMC
            DO IY=1,NY
             DO IX=1,NX
              IF(FSM(IX,IY).EQ.1.)  BATM(IX,IY,ISYS)=BATMC
             ENDDO
            ENDDO
          ELSE
            DO IY=1,NY
             READ(36) (BATM(IX,IY,ISYS),IX=1,NX)
             DO IX=1,NX
              IF(FSM(IX,IY).NE.1.)  BATM(IX,IY,ISYS)=0.0
             ENDDO
            ENDDO
          ENDIF
        ENDIF
        IF(LIST(3).EQ.1)  THEN
         IF(NOATM(ISYS).EQ.1) THEN
           WRITE(OUT,5070)  NXATMT,BATMC
 5070      FORMAT(40X,'ATMOSPHERIC LOADING AT TIME ',F6.1' EQUALS',
     .            E10.3)
         ELSE
           WRITE(OUT,5075)  NXATMT
 5075      FORMAT(40X,'ATMOSPHERIC LOADINGS AT TIME ',F6.1/3X,'IY')
           DO IY=1,NY
            WRITE(OUT,5060)  IY,(BATM(IX,IY,ISYS),IX=1,NX)
           ENDDO
         ENDIF
        ENDIF
        DO IY=1,NY
         DO IX=1,NX
           BATM(IX,IY,ISYS) = 1000.*SCALATM(ISYS)*BATM(IX,IY,ISYS) 
         ENDDO
        ENDDO

        TOTLD = 0.0
        DO 360 IY=1,NY
         DO 360 IX=1,NX
           TOTLD = TOTLD + XAZ(IX,IY)*BATM(IX,IY,ISYS)
  360   CONTINUE
        WRITE(OUT,5080)  TIME,TOTLD/1000.
 5080   FORMAT(/36X,'TOTAL ATMOSPHERIC LOADING AT TIME,',F7.3,
     .   ' =',E13.3,' KG/DAY')

      ENDIF

      GO TO 380 

  370 WRITE(OUT,5090)    ISYS 
 5090 FORMAT(///40X,'NO ATMOSPHERIC LOADINGS FOR SYSTEM',I3) 
      WRITE(OUT,1150)   SYNAME(ISYS)

  380 CONTINUE

      IF(IATMOPT.EQ.1)  RETURN
      IF(INPCHCK.EQ.0) THEN
       IF(IATMOPT.EQ.2) THEN
        IF(IATMPWLOPT.EQ.0)  THEN
         IF(IBNRYRDOPTS(5).EQ.0)  THEN
           READ(36,1100,ERR=936)  NXATMTSECS
         ELSE
           READ(36)  NXATMTSECS
         ENDIF
         NXATMTSECS=ISCALATM*NXATMTSECS
         NXATMT=NXATMTSECS/86400.
        ELSE
         CALL RCA10(SATM,BATM,NX*NY,NOATM,NXATMT,36,IATMPWLOPT,
     .              SCALATM)
        ENDIF
       ENDIF

      ELSE

  390  IF(IATMPWLOPT.EQ.0)  THEN
         IF(IBNRYRDOPTS(5).EQ.0)  THEN
           READ(36,1100,ERR=936)  NXATMTSECS
         ELSE
           READ(36)   NXATMTSECS
         ENDIF
         NXATMTSECS=ISCALATM*NXATMTSECS
         NXATMT=NXATMTSECS/86400.
       ENDIF
       CALL RCA10(SATM,BATM,MXWK,NOATM,NXATMT,36,IATMPWLOPT,SCALATM)
       BACKSPACE(36)
       IF(INPCHCK.EQ.1)  GO TO 390
       INPCHCK=1

      ENDIF

  400 RETURN

  933 IN=33
      GO TO 950
  934 IN=34
      GO TO 950
  935 IN=35
      GO TO 950
  936 IN=36
  950 CALL FMTER
      CALL EXIT 
  943 IN=43
      GO TO 955
  944 IN=44
      GO TO 955
  945 IN=45
  955 WRITE(OUT,9550)  IN
 9550 FORMAT(///20X,'ERROR READING BINARY FILE NUMBER ',I3)
      CALL EXIT 
  960 WRITE(OUT,9600)  
 9600 FORMAT(//11X,'ERROR...USER REQUESTING MORE POINT SOURCE LOADS THAN
     . DIMENSIONED FOR IN -RCACM-'/11X,'RCA TERMINATED'//)
      CALL EXIT 
  962 WRITE(OUT,9602)  
 9602 FORMAT(//11X,'ERROR...USER REQUESTING MORE NONPOINT SOURCE LOADS T
     .HAN DIMENSIONED FOR IN -RCACM-'/11X,'RCA TERMINATED'//)
      CALL EXIT 
  964 WRITE(OUT,9604)  
 9604 FORMAT(//11X,'ERROR...USER REQUESTING MORE FALL-LINE LOADS THAN DI
     .MENSIONED FOR IN -RCACM-'/11X,'RCA TERMINATED'//)
      CALL EXIT 
  965 WRITE(OUT,9605)  
 9605 FORMAT(//
     .  11X,'ERROR...REQUESTED PS LOAD NOT FOUND IN LOADING TABLE'/
     .  11X,'SYSTEM =',I3,2X,'IWK =',I5,2X,'IPS(IWK) =',I5/
     .  11X,'RCA TERMINATED'//)
      CALL EXIT 
  970 WRITE(OUT,9700)  
 9700 FORMAT(//
     .  11X,'ERROR...REQUESTED NPS LOAD NOT FOUND IN LOADING TABLE'/
     .  11X,'SYSTEM =',I3,2X,'IWK =',I5,2X,'IPS(IWK) =',I5/
     .  11X,'RCA TERMINATED'//)
      CALL EXIT 
  975 WRITE(OUT,9705)  
 9705 FORMAT(//
     .  11X,'ERROR...REQUESTED FL LOAD NOT FOUND IN LOADING TABLE'/
     .  11X,'SYSTEM =',I3,2X,'IWK =',I5,2X,'IPS(IWK) =',I5/
     .  11X,'RCA TERMINATED'//)
      CALL EXIT 
  980 WRITE(OUT,9800)  ISYS,I  
 9800 FORMAT(//11X,'LOADING ERROR IN SYSTEM ',I2/
     .  11X,'ILLEGAL SEGMENT NUMBER SPECIFIED FOR LOAD NUMBER ',I4/
     .  11X,'RCA TERMINATED'//)
      CALL EXIT 
  990 WRITE(OUT,9900)  ISYS,I,IX,IY
 9900 FORMAT(//11X,'LOADING ERROR IN SYSTEM ',I2/
     .  11X,'LOAD NUMBER ',I4,' BEING PUT INTO A NON-WATER SEGMENT (IX,
     .IY =',2I5,')'/11X,'RCA TERMINATED'//)
      CALL EXIT 
      END 
