      SUBROUTINE ncset_oceantime(ncID,TREC,oceant,not_enough_len)
!
!=======================================================================
! Create a NetCDF file ready for RCA output                            !
!                            YUN LI, UMCES/HPL Sep-30-2008             !
!=======================================================================
!
      USE netcdf
      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
! 
!  Imported and exported variable declarations.
!
      INTEGER,     INTENT(in)  ::     ncID
      INTEGER,     INTENT(in)  ::     TREC
      REAL,        INTENT(out) ::     oceant
      LOGICAL,     INTENT(out) ::     not_enough_len
!
!  Local variable declarations
!
      INTEGER         ::     totalREC
      REAL            ::     oceant_tem(1)
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncset_oceantime.f'
!
!-----------------------------------------------------------------------
! check ocean_time length
!-----------------------------------------------------------------------
!
      not_enough_len=.FALSE.
      status=nf90_inq_dimid(ncID,'ocean_time',dimID)
      CALL nccheck_status(status,'ocean_time',SUBNA)
      status=nf90_inquire_dimension(ncID,dimID,len=totalREC)
      CALL nccheck_status(status,'ocean_time',SUBNA)
      IF(totalREC.LT.TREC) THEN 
        not_enough_len=.TRUE.
        oceant=-9.9e37
        RETURN
      ENDIF
!
!-----------------------------------------------------------------------
! read ocean_time
!-----------------------------------------------------------------------
!
      status=nf90_inq_varid(ncID,'ocean_time',varID)
      CALL nccheck_status(status,'ocean_time',SUBNA)
      status=nf90_get_var(ncID,varID,oceant_tem,
     .                      start=(/TREC/),
     .                      count=(/1/))
      CALL nccheck_status(status,'ocean_time',SUBNA)
      oceant=oceant_tem(1)

      END SUBROUTINE ncset_oceantime
