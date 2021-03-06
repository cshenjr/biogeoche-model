      SUBROUTINE RCA09
C 
C        RCA09 IS AN OUTPUT ROUTINE FOR INTERMEDIATE DUMPS
C 
      SAVE
      INCLUDE 'RCACM' 
      INCLUDE 'NetCDFCM' 

      REAL  TOTLOADS(4)
      EQUIVALENCE
     .      (IDUMP(1,1),I1) , (IDUMP(1,2),J1) , (IDUMP(1,3),K1) 
     .   ,  (IDUMP(2,1),I2) , (IDUMP(2,2),J2) , (IDUMP(2,3),K2) 
     .   ,  (IDUMP(3,1),I3) , (IDUMP(3,2),J3) , (IDUMP(3,3),K3) 
     .   ,  (IDUMP(4,1),I4) , (IDUMP(4,2),J4) , (IDUMP(4,3),K4) 
     .   ,  (IDUMP(5,1),I5) , (IDUMP(5,2),J5) , (IDUMP(5,3),K5) 
     .   ,  (IDUMP(6,1),I6) , (IDUMP(6,2),J6) , (IDUMP(6,3),K6) 

C        INCREMENT OUTPUT RECORD NUMBER 

      PT = ITIMESECS/86400.
      IREC=IREC+1 

      IF(IREC.EQ.1)   THEN
         WRITE(OUT,2000)  PT,((IDUMP(I,J),J=1,3),I=1,6)
 2000    FORMAT('1'//6X,'TIME =',F10.4,' DAYS',5X,'<--  SEGMENTS  -->'
     .     /1X,'SYSTEM',8X,6(I5,',',I3,',',I2)/) 
      ELSE
         IF(NOSYS.GT.20) WRITE(OUT,2001)  PT
 2001    FORMAT('1'//6X,'TIME =',F10.4,' DAYS')
         IF(NOSYS.LE.20) WRITE(OUT,2002)  PT
 2002    FORMAT(//6X,'TIME =',F10.4,' DAYS')
      ENDIF

      DO 20 ISYS=1,NOSYS
  20  WRITE(OUT,1000) ISYS,SYNAME(ISYS),
     .   CARAY(I1,J1,K1,ISYS),CARAY(I2,J2,K2,ISYS),
     .   CARAY(I3,J3,K3,ISYS),CARAY(I4,J4,K4,ISYS),
     .   CARAY(I5,J5,K5,ISYS),CARAY(I6,J6,K6,ISYS)
 1000 FORMAT(1X,I3,1X,A8,1X,6E12.4) 

      WRITE(OUT,3000)  PT
 3000 FORMAT(//10X,'LOADING RATES (KG/DAY) BY SYSTEM FOR TIME =',
     .  F10.4,' DAYS'/19X,'POINT     NONPOINT'/5X,
     .  'SYSTEM       SOURCE       SOURCE    FALL-LINE  ATMOSPHERIC')
      DO 200 ISYS=1,NOSYS
       DO ILOAD=1,4
        TOTLOADS(ILOAD)=0.0
       ENDDO
       IF(NOPS(ISYS).GT.0)  THEN
        DELTPS = TIME - NXPST
        DO 145 I=1,NOPS(ISYS) 
         TOTLOADS(1) = TOTLOADS(1) + DELTPS*SPS(I,ISYS) + BPS(I,ISYS)
  145   CONTINUE
       ENDIF
       IF(NONPS(ISYS).GT.0)   THEN
        DELTNPS = TIME - NXNPST
        DO 150 I=1,NONPS(ISYS)
         TOTLOADS(2) = TOTLOADS(2) + DELTNPS*SNPS(I,ISYS) + BNPS(I,ISYS)
  150   CONTINUE
       ENDIF
       IF(NOFL(ISYS).GT.0)   THEN
        DELTFL = TIME - NXFLT
        DO 155 I=1,NOFL(ISYS)
          TOTLOADS(3) = TOTLOADS(3) + DELTFL*SFL(I,ISYS) + BFL(I,ISYS)
  155   CONTINUE
       ENDIF
       IF(NOATM(ISYS).GT.0)   THEN
        DELTATM = TIME - NXATMT
        DO 160 IY=1,NY
         IS = IXS(IY)
         IE = IXE(IY)
         IF(IS.EQ.0)   GO TO 160
         DO 157 IX=IS,IE
          IF(FSM(IX,IY).LE.0.)  GO TO 157
          TOTLOADS(4) = TOTLOADS(4) + (DELTATM*SATM(IX,IY,ISYS)
     .                + BATM(IX,IY,ISYS))*XAZ(IX,IY)
  157    CONTINUE
  160   CONTINUE
       ENDIF
       IF(TOTLOADS(1).EQ.0.0 .AND. TOTLOADS(2).EQ.0.0 .AND.
     .    TOTLOADS(3).EQ.0.0 .AND. TOTLOADS(4).EQ.0.0) GO TO 200
       WRITE(OUT,3100)  ISYS,(TOTLOADS(ILOAD)/1000.,ILOAD=1,4)
 3100  FORMAT(8X,I3,4E13.4)
  200 CONTINUE

!
!-------------------------------------------------------
! Optionally write information file 
!-------------------------------------------------------
!
      IF(NETCDFOPT.EQ.1) THEN
      ELSE
       REWIND(15)
       WRITE(15)   CARAY
      ENDIF

      RETURN
      END 
