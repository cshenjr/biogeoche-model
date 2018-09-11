      SUBROUTINE nccheck_status(status, vname, filename)
!
!=======================================================================
! check NetCDF status and print error message if possible              !
!                                      YUN LI, UMCES/HPL Sep-30-2008   !
!=======================================================================
!
      USE netcdf
!   
      INTEGER,   INTENT(in)    ::     status
      CHARACTER(*), INTENT(in) ::     vname
      CHARACTER(*), INTENT(in) ::     filename
!
      IF(status/=nf90_noerror) THEN 
         WRITE(*,100) TRIM(ADJUSTL(filename)),TRIM(ADJUSTL(vname)), 
     .                nf90_strerror(status) 
         CALL EXIT 
      ENDIF
100   format(3X,'==> NetCDF ERROR in ',A,' ::  ',A,1X,A)
      END SUBROUTINE nccheck_status

