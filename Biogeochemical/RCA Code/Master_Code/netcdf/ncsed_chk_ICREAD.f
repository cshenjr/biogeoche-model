      SUBROUTINE ncsed_chk_ICREAD(ISEDPRNT,idvar,Vvalue)
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
!
!  Local variable declarations
!
      INTEGER                ::   IX, IY
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncsed_chk_ICREAD.f'
!      
!-----------------------------------------------------------------------
! write information of sediment variable
!-----------------------------------------------------------------------
!
      WRITE(OUT,1200)  TRIM(ADJUSTL(Vinfo(1,idvar))),Fscale(idvar)
 1200 FORMAT(//50X,'SEDIMENT INITIAL CONDITIONS FOR',/53X,A25/
     .   53X,'SCALE FACTOR =',E10.3//3X,'COLUMN      ROW -->'/
     .   10X,7X,'1',11X,'2',11X,'3',11X,'4',11X,'5',11X,'6',11X,'7'
     .   11X,'8',11X,'9',10X,'10'/)
      DO 150 IY=1,NY
       IF(ISEDPRNT.GT.0) WRITE(OUT,1300)  IY,(Vvalue(IX,IY),IX=1,NX)
 1300  FORMAT(4X,I3,3X,10E12.3,/(10X,10E12.3))
       DO 130 IX=1,NX
         Vvalue(IX,IY)=Fscale(idvar)*Vvalue(IX,IY)
  130  CONTINUE
  150 CONTINUE
      RETURN

      END SUBROUTINE ncsed_chk_ICREAD
