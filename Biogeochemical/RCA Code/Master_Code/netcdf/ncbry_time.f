      SUBROUTINE ncbry_time(ncID)
!
!=======================================================================
!  Read RCA boundary time in NetCDF format                             !
!                                 YUN LI, UMCES/HPL, Feb-22-2011       !
!=======================================================================
!
      USE netcdf
      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
! 
!  Imported variable declarations.
!
      INTEGER, INTENT(in)  ::   ncID        ! identification of netcdf file
!
!  Local variable declarations
!
      INTEGER              ::   idvar 
      INTEGER              ::   mstart(1), mcount(1)
      REAL*8               ::   Vtem(1)
      CHARACTER            ::   Vtwarp*4
      INTEGER              ::   Vtscal
!
!----------------------------------------------------------------------
! name of this subroutine
!----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncbry_time.f'
!
!----------------------------------------------------------------------
! check input type and set variable dimension
!----------------------------------------------------------------------
!
      idvar=idbry_time 
      IF(Iinfo(idvar).EQ.t1dvar) THEN
        mstart=(/bryTflag/)
        mcount=(/1/)
      ELSE
        WRITE(OUT,*) TRIM(ADJUSTL(SUBNA))
        WRITE(OUT,*) 'ERROR: misuse of ncbry_time.f ',
     .               TRIM(ADJUSTL(Vinfo(1,idvar))),
     .               ' is not 1D time-dependent variable'
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
      status=nf90_get_var(ncID,VarID,Vtem,start=mstart,count=mcount)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)
      Vtwarp=TRIM(ADJUSTL(Vinfo(3,idbry_time)))
      Vtscal=IUNITCHECK(Vtwarp,TRIM(ADJUSTL(Vinfo(3,idbry_time))))
      NXBCTSECS=int(Vtscal*Vtem(1))
      TWARPBC='SECS'

      END SUBROUTINE ncbry_time
