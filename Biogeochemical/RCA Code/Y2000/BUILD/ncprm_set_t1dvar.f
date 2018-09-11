      SUBROUTINE ncprm_set_t1dvar(ncID,idvar_time,idvar,
     .                            nbreak,Vtwarp,Vvalue,Vtime)
!
!=======================================================================
!  Read 1D time-dependent parameters from NetCDF file                  !
!                                 YUN LI, UMCES/HPL, May-12-2011       !
!=======================================================================
!
      USE netcdf
      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
! 
!  Imported variable declarations.
!
      INTEGER,   INTENT(in)      ::   ncID        ! identification of netcdf file
      INTEGER,   INTENT(in)      ::   idvar_time  ! variable identification (see NetCDFCM)
      INTEGER,   INTENT(in)      ::   idvar       ! variable identification (see NetCDFCM)
      INTEGER,   INTENT(out)     ::   nbreak
      CHARACTER, INTENT(out)     ::   Vtwarp*4
      REAL,      INTENT(out)     ::   Vvalue(MXFUNCT)
      REAL,      INTENT(out)     ::   Vtime(MXFUNCT)
!
!  Local variable declarations
!
      INTEGER               ::   mstart(1), mcount(1)
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncprm_set_t1dvar.f'
!
!----------------------------------------------------------------------
! check input type and set variable dimension
!----------------------------------------------------------------------
!
      IF(Iinfo(idvar).EQ.t1dvar) THEN 
      ELSE
        WRITE(OUT,*) TRIM(ADJUSTL(SUBNA))
        WRITE(OUT,*) 'ERROR: misuse of ncprm_set_t1dvar.f ',
     .               TRIM(ADJUSTL(Vinfo(1,idvar))),
     .               ' is not 1D time-independent variable'
        CALL EXIT
      ENDIF

! inquire variable dimension
      status=nf90_inq_dimid(ncID,
     .       TRIM(ADJUSTL(Vinfo(1,idvar_time))),dimID)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      status=nf90_inquire_dimension(ncID,dimID,len=nbreak)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      Vvalue=-9999.
      Vtime=-9999.
      mstart=(/1/)
      mcount=(/nbreak/)

! inquire variable time
      status=nf90_inq_varid(ncID,TRIM(ADJUSTL(Vinfo(1,idvar_time))),
     .                      varID)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      status=nf90_get_var(ncID,VarID,Vtime(1:nbreak),
     .                    start=mstart,
     .                    count=mcount)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      Vtwarp = TRIM(ADJUSTL(Vinfo(3,idvar_time)))

! inquire variable
      status=nf90_inq_varid(ncID,TRIM(ADJUSTL(Vinfo(1,idvar))),varID)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

! read values
      status=nf90_get_var(ncID,varID,Vvalue(1:nbreak),
     .                    start=mstart,
     .                    count=mcount)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

      END SUBROUTINE ncprm_set_t1dvar
