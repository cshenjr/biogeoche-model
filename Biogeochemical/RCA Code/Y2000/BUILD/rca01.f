      SUBROUTINE RCA01
C
C        RCA01 IS THE INITILIZATION ROUTINE
C
      SAVE
      CHARACTER   TITLE*80,TYPE*15
      INTEGER*2  SYSGDP(40)
      INCLUDE  'RCACM'
      INCLUDE  'NetCDFCM'

C        INITIAL -SYSBY- SETTINGS
      DO I=1,40
       SYSGDP(I)=1
      ENDDO

C        DEFINE I/O DEVICES
      IN = 5
      OUT = 6
 
C        THE -LIST- ARRAY CONTROLS THE LISTING OF THE USERS INPUT DECK
C         LIST(#)        CONTROLS THE LISTING OF
C              1   GEOMETRY...EXCHANGES, VOLUMES AND FLOWS
C              2   BOUNDARY CONDITIONS
C              3   LOADINGS
C              4   PARAMETERS, CONSTANTS AND MISC TIME FNS
C              5   INITIAL CONDITIONS
C        LIST = 1   LIST INPUT
C        LIST = 0   DO NOT LIST INPUT
C
      READ(IN,1000)  COMMENT
 1000 FORMAT(A)
      READ(IN,1100,ERR=950)   CYCLE,LIST,IDIAGDT,INPCHCK
 1100 FORMAT(8I10)
      WRITE(OUT,1500)  NX,NY,NZ,NOSYS
 1500 FORMAT( 1H1/1X,120('*')/1X,120('*')////
     .   25X,'R C A 3 D   - -  R O W / C O L U M N   A E S O P'/
     .   25X,'                  ( 3 - D I M E N S I O N A L )'//
     .   70X,'Version 3.00  Released June 30,2004'/
     .   70X,'DEVELOPED BY HYDROQUAL INC.  MAHWAH, NJ'/
     .   70X,'THIS EXECUTABLE COMPILED FOR:  NX =',I4/
     .   97X,'    NY =',I4/97X,'    NZ =',I4/97X,' NOSYS =',I4/
     .   ////1X,120('*')/1X,120('*')//)
      IF(CYCLE.EQ.0)  THEN
         WRITE(OUT,1600)
 1600    FORMAT(
     .     15X,'INITIAL CONDITIONS WILL BE READ FROM USER INPUT DECK'//)
      ELSE
!
!---------------------------------------------------------------------
!  Optionally Read intial condition when restart
!---------------------------------------------------------------------
!
         IF(NETCDFOPT.EQ.1) THEN
         WRITE(OUT,1602)
 1602    FORMAT(
     .     15X,'INITIAL CONDITIONS WILL BE READ FROM RESTART FILE', 
     .     '(_rst.nc)'//)
         ELSE
         WRITE(OUT,1601)
 1601    FORMAT(
     .     15X,'INITIAL CONDITIONS WILL BE READ FROM RCAFIC'//)
         ENDIF
      ENDIF

      READ(IN,1000) COMMENT
   10 READ(IN,1000) TITLE
      IF(TITLE(1:3).NE.'END' .AND. TITLE(1:3).NE.'end')   THEN
         WRITE(OUT,1110)   TITLE
 1110    FORMAT(10X,A80)
         GO TO 10
      ENDIF
      WRITE(OUT,2100)
 2100 FORMAT(    ////1X,120('*')/1X,120('*')////)

C        READ STATE VARIABLE NAMES
      READ(IN,1000) COMMENT
      READ(IN,1120,ERR=950)  (SYNAME(I),I=1,NOSYS)
 1120 FORMAT(10A8)
      WRITE(OUT,2250)  (I,SYNAME(I),I=1,NOSYS)
 2250 FORMAT(5X,'THE FOLLOWING STATE-VARIABLES ARE BEING MODELED'/ 20X,
     .   'SYSTEM     VARIABLE'/20X,6('-'),5X,8('-')/(20X,I4,'.',6X,A8))

C     System Bypass Options
      READ(IN,1000) COMMENT
      READ(IN,3000,ERR=950)   (SYSBY(I),I=1,NOSYS)
 3000 FORMAT(40I2)
      WRITE(OUT,4000)    NOSYS,(SYSBY(I),I=1,NOSYS)
 4000 FORMAT(15X,'SYSTEM BYPASS OPTIONS FOR SYSTEMS 1 TO',I3,' ARE',/
     .  15X,50I2)

C        DETERMINE SYSTEM NUMBER OF LAST ACTIVE SYSTEM
      DO 20 ISYS=1,NOSYS
        SYSGDP(ISYS) = SYSBY(ISYS)
        IF(SYSBY(ISYS).EQ.0)   MXACTS = ISYS
   20 CONTINUE
!
!-------------------------------------------------------
! Optionally write information file 
!-------------------------------------------------------
!
      IF(NETCDFOPT.EQ.1) THEN
      ELSE
       WRITE(10)   NX,NY,NZ,NOSYS
       WRITE(10)   SYNAME
       WRITE(10)   SYSGDP
      ENDIF

      WRITE(OUT,4500)
 4500 FORMAT(////1X,119('*')/1X,119('*'))

C        OTHER INITIALIZATION
      TIME=0.0
      INITB = 0
      IREC=0
      IREC_sed=0  !! Yun Li, UMCES/HPL, Jun-01-2011
      rstTflag=0  !! Yun Li, UMCES/HPL, Jun-03-2011

      DO 40 ISYS=1,NOSYS
        DO 25 IZ=1,NZ
        DO 25 IY=1,NY
        DO 25 IX=1,NX
          CDARAY(IX,IY,IZ,ISYS) = 0.
          CARAY(IX,IY,IZ,ISYS) = 0.
   25   CONTINUE
        DO 30 IY=1,NY
        DO 30 IX=1,NX
          BATM(IX,IY,ISYS) = 0.
          SATM(IX,IY,ISYS) = 0.
   30   CONTINUE
        DO 35 I=1,MXBC
          BBC(I,ISYS) = 0.
          SBC(I,ISYS) = 0.
   35   CONTINUE
        DO 38 I=1,MXWK
         BPS(I,ISYS) = 0.
         SPS(I,ISYS) = 0.
         BNPS(I,ISYS) = 0.
         SNPS(I,ISYS) = 0.
         BFL(I,ISYS) = 0.
         SFL(I,ISYS) = 0.
   38   CONTINUE
   40 CONTINUE

      DO 60 IZ=1,NZ+1
      DO 60 IY=1,NY
      DO 60 IX=1,NX
        IF(IZ.LE.NZ)  THEN
          VDER(IX,IY,IZ) = 0.
          QX(IX,IY,IZ) = 0.
          QY(IX,IY,IZ) = 0.
          RX(IX,IY,IZ) = 0.
          RY(IX,IY,IZ) = 0.
        ENDIF
        QZ(IX,IY,IZ) = 0.
        RZ(IX,IY,IZ) = 0.
   60 CONTINUE

      RETURN

 950  CALL FMTER
      CALL EXIT
      END
