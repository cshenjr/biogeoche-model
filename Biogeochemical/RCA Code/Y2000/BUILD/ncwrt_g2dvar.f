      SUBROUTINE ncwrt_g2dvar(ncID,idvar,Vvalue)
!
!=======================================================================
!  WRITE RCA variables into NetCDF output                              !
!                                 YUN LI, UMCES/HPL, Feb-17-2011       !
!=======================================================================
!
      USE netcdf
      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
! 
!  Imported variable declarations.
!
      INTEGER,          INTENT(in)   ::   ncID        ! identification of netcdf file
      INTEGER,          INTENT(in)   ::   idvar       ! variable identification (see NetCDFCM)
      REAL,             INTENT(in)   ::   Vvalue(NX,NY)
!
!  Local variable declarations
!
      INCLUDE 'ncdim_info.h'
      INTEGER              ::   IX, IY 
      INTEGER              ::   mstart(3), mcount(3), mdimen(3)
      REAL                 ::   Vvalue_mask(NX,NY)
      CHARACTER(100)       ::   coordinates
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncbry_g2dvar.f'
!
!----------------------------------------------------------------------
! inqure dimension information
!----------------------------------------------------------------------
!
      SUBNA_h=TRIM(ADJUSTL(SUBNA))//' -> ncinq_dim.h'
      INCLUDE 'ncinq_dim.h'
!
!----------------------------------------------------------------------
! check input type and set variable dimension
!----------------------------------------------------------------------
!
      IF(Iinfo(idvar).EQ.g2dvar) THEN 
        coordinates='x_r y_r '//Vinfo(5,idvar)
        mdimen=(/NXdimID, NYdimID, TdimID/)
        mcount=(/NX,NY,1/)
        mstart=(/ 1, 1,ncIREC/)
      ELSE
        WRITE(OUT,*) TRIM(ADJUSTL(SUBNA_h))
        WRITE(OUT,*) 'ERROR: misuse of ncwrt_g2dvar.f ',
     .               TRIM(ADJUSTL(Vinfo(1,idvar))),
     .               ' is not 2D time-dependent variable'
        CALL EXIT
      ENDIF

! define/inquire variable
      IF(ncIREC.EQ.1) THEN
        SUBNA_h=TRIM(ADJUSTL(SUBNA))//' -> ncdef_var.h'
        INCLUDE 'ncdef_var.h'
      ENDIF
      status=nf90_inq_varid(ncID,TRIM(ADJUSTL(Vinfo(1,idvar))),varID)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

! set land values to default spval (for ncview display)
      DO IX=1,NX
      DO IY=1,NY
        IF(FSM(IX,IY).EQ.0.0) THEN
           Vvalue_mask(IX,IY)=spval
        ELSE
           Vvalue_mask(IX,IY)=Vvalue(IX,IY)
        ENDIF
      ENDDO
      ENDDO

! put values
      status=nf90_put_var(ncID,varID,Vvalue_mask,
     .                    start=mstart,
     .                    count=mcount)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

      END SUBROUTINE ncwrt_g2dvar
