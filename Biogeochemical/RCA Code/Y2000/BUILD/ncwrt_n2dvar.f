      SUBROUTINE ncwrt_n2dvar(ncID,idvar,Vvalue)
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
      INTEGER              ::   mstart(2), mcount(2), mdimen(2)
      CHARACTER(100)       ::   coordinates
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncbry_n2dvar.f'
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
      IF(Iinfo(idvar).EQ.n2dvar) THEN
        coordinates='x_r y_r'
        mdimen=(/NXdimID, NYdimID/)
        mcount=(/NX,NY/)
        mstart=(/ 1, 1/)
      ELSE
        WRITE(OUT,*) TRIM(ADJUSTL(SUBNA_h))
        WRITE(OUT,*) 'ERROR: misuse of ncwrt_n2dvar.f ',
     .               TRIM(ADJUSTL(Vinfo(1,idvar))),
     .               ' is not 2D time-independent variable'
        CALL EXIT
      ENDIF

! define/inquire variable
      SUBNA_h=TRIM(ADJUSTL(SUBNA))//' -> ncdef_var.h'
      INCLUDE 'ncdef_var.h'
      status=nf90_inq_varid(ncID,TRIM(ADJUSTL(Vinfo(1,idvar))),varID)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

! put values
      status=nf90_put_var(ncID,varID,Vvalue,start=mstart,count=mcount)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

      END SUBROUTINE ncwrt_n2dvar
