      SUBROUTINE ncbry_BCOPT(ncID,idvar,Vvalue)
!
!=======================================================================
!  Read RCA boundary options in NetCDF format                          !
!                                 YUN LI, UMCES/HPL, Feb-22-2011       !
!=======================================================================
!
      USE netcdf
      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
!
!  Imported variable declarations.
!
      INTEGER, INTENT(in)   ::   ncID        ! identification of netcdf file
      INTEGER, INTENT(in)   ::   idvar       ! variable identification (see NetCDFCM)
      INTEGER, INTENT(out)  ::   Vvalue
!
!  Local variable declarations
!
      INTEGER               ::   mstart(1), mcount(1), Vtem(1)
!
!----------------------------------------------------------------------
! name of this subroutine
!----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncbry_BCOPT.f'
!
!----------------------------------------------------------------------
! check input type and set variable dimension
!----------------------------------------------------------------------
!
      IF(Iinfo(idvar).EQ.n0dvar) THEN
        mstart=(/1/)
        mcount=(/1/)
      ELSE
        WRITE(OUT,*) TRIM(ADJUSTL(SUBNA))
        WRITE(OUT,*) 'ERROR: misuse of ncbry_BCOPT.f ',
     .               TRIM(ADJUSTL(Vinfo(1,idvar))),
     .               ' is not a constant'
        CALL EXIT 
      ENDIF
!
!----------------------------------------------------------------------
! define/inquire variable
!----------------------------------------------------------------------
!
      status=nf90_inq_varid(ncID,Vinfo(1,idvar),varID)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)
!
!----------------------------------------------------------------------
! read values
!----------------------------------------------------------------------
!
      status=nf90_get_var(ncID,varID,Vtem,start=mstart,count=mcount)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)
      Vvalue=Vtem(1)

      END SUBROUTINE ncbry_BCOPT
