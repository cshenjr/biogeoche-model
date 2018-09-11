      SUBROUTINE ncprm_set_TSS(ncID,idvar_time,idvar,itrec,
     .                         nbreak,Vtime,Vvalue)
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
      INTEGER,   INTENT(inout)   ::   idvar_time  ! variable identification (see NetCDFCM)
      INTEGER,   INTENT(in)      ::   idvar       ! variable identification (see NetCDFCM)
      INTEGER,   INTENT(in)      ::   itrec
      INTEGER,   INTENT(out)     ::   nbreak
      REAL,      INTENT(out)     ::   Vtime(7320)
      REAL,      INTENT(out)     ::   Vvalue(NX,NY)
!
!  Local variable declarations
!
      INTEGER               ::   mstart(3), mcount(3)
      CHARACTER             ::   Vtwarp*4
      INTEGER               ::   Vtscal
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncprm_set_TSS.f'
!
!----------------------------------------------------------------------
! check input type and set variable dimension
!----------------------------------------------------------------------
!
      IF(Iinfo(idvar).EQ.t2dvar) THEN 
      ELSE
        WRITE(OUT,*) TRIM(ADJUSTL(SUBNA))
        WRITE(OUT,*) 'ERROR: misuse of ncprm_set_TSS.f ',
     .               TRIM(ADJUSTL(Vinfo(1,idvar))),
     .               ' is not 2D time-dependent variable'
        CALL EXIT
      ENDIF

! inquire variable dimension
      status=nf90_inq_dimid(ncID,
     .       TRIM(ADJUSTL(Vinfo(1,idvar_time))),dimID)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      status=nf90_inquire_dimension(ncID,dimID,len=nbreak)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      mstart=(/1,1,itrec/)
      mcount=(/NX,NY,1/)

! inquire variable time
      IF(itrec.EQ.1) THEN
        status=nf90_inq_varid(ncID,TRIM(ADJUSTL(Vinfo(1,idvar_time))),
     .                        varID)
        CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
        status=nf90_get_var(ncID,VarID,Vtime(1:nbreak),
     .                      start=(/1/),
     .                      count=(/nbreak/))
        CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
        Vtwarp=TRIM(ADJUSTL(Vinfo(3,idvar_time)))
        Vtscal=IUNITCHECK(Vtwarp,TRIM(ADJUSTL(Vinfo(1,idvar_time))))
      ENDIF

! inquire variable
      status=nf90_inq_varid(ncID,TRIM(ADJUSTL(Vinfo(1,idvar))),varID)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

! read values
      status=nf90_get_var(ncID,varID,Vvalue,
     .                    start=mstart,
     .                    count=mcount)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

      END SUBROUTINE ncprm_set_TSS
