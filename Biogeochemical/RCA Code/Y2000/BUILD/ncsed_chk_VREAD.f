      SUBROUTINE ncsed_chk_VREAD(ISEDPRNT,idvar,Vvalue,INUM)
!
!=======================================================================
!  Check print options and scale of sediment variables                 !
!                                 YUN LI, UMCES/HPL, Feb-17-2011       !
!=======================================================================
!
      USE netcdf
      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
      INCLUDE 'SEDCM'
! 
!  Imported variable declarations.
!
      INTEGER, INTENT(in)    ::   ISEDPRNT
      INTEGER, INTENT(in)    ::   idvar       ! variable identification (see NetCDFCM)
      REAL,    INTENT(inout) ::   Vvalue(NX,NY)
      INTEGER, INTENT(in)    ::   INUM
!
!  Local variable declarations
!
      INTEGER                ::   IX, IY
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncsed_chk_VREAD.f'
!      
!-----------------------------------------------------------------------
! write information of sediment variable
!-----------------------------------------------------------------------
!
      IF(INUM.EQ.1)   WRITE(OUT,1250)
 1250 FORMAT(//52X,'SEDIMENTATION RATES (CM/YR)'/)
      IF(INUM.EQ.2)   WRITE(OUT,1260)
 1260 FORMAT(//43X,'SEDIMENT SOLID-PHASE MIXING RATES (M**2/DAY)'/)
      IF(INUM.EQ.3)   WRITE(OUT,1270)
 1270 FORMAT(//41X,'SEDIMENT DISSOLVED-PHASE MIXING RATES (M**2/DAY)'/)
      IF(INUM.EQ.4)   WRITE(OUT,1271)
 1271 FORMAT(//41X,'SEDIMENT SOLIDS CONCENTRATION - AEROBIC (KG/L)'/)
      IF(INUM.EQ.5)   WRITE(OUT,1272)
 1272 FORMAT(//41X,'SEDIMENT SOLIDS CONCENTRATION - ANAEROBIC (KG/L)'/)
      IF(INUM.EQ.6)   WRITE(OUT,1273)
 1273 FORMAT(//41X,'AEROBIC LAYER SI PARTITION COEFFICIENT (L/KG)'/)
      IF(INUM.EQ.7)   WRITE(OUT,1274)
 1274 FORMAT(//41X,'ANAEROBIC LAYER SI PARTITION COEFFICIENT (L/KG)'/)
      IF(INUM.EQ.8)   WRITE(OUT,1275)
 1275 FORMAT(//41X,'AEROBIC LAYER PO4 PARTITION COEFFICIENT (L/KG)'/)
      IF(INUM.EQ.9)   WRITE(OUT,1276)
 1276 FORMAT(//41X,'ANAEROBIC LAYER PO4 PARTITION COEFFICIENT (L/KG)'/)
      WRITE(OUT,1300)  Fscale(idvar)
 1300 FORMAT(53X,'SCALE FACTOR =',E10.3//3X,'COLUMN      ROW -->'/
     .   10X,7X,'1',11X,'2',11X,'3',11X,'4',11X,'5',11X,'6',11X,'7'
     .   11X,'8',11X,'9',10X,'10'/)
      DO 150 IY=1,NY
       IF(ISEDPRNT.GT.0) WRITE(OUT,1350)  IY,(Vvalue(IX,IY),IX=1,NX)
 1350  FORMAT(4X,I3,3X,10E12.3,/(10X,10E12.3))
       DO 130 IX=1,NX
         Vvalue(IX,IY)=Fscale(idvar)*Vvalue(IX,IY)
  130  CONTINUE
  150 CONTINUE
      RETURN

      END SUBROUTINE ncsed_chk_VREAD
