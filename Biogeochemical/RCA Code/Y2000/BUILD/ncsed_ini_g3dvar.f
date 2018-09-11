      SUBROUTINE ncsed_ini_g3dvar(ncID,idvar,Vvalue)
!
!=======================================================================
!  Read sediment parameters from NetCDF file                           !
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
      INTEGER, INTENT(in)   ::   ncID        ! identification of netcdf file
      INTEGER, INTENT(in)   ::   idvar       ! variable identification (see NetCDFCM)
      REAL,    INTENT(out)  ::   Vvalue(NX,NY,3)
!
!  Local variable declarations
!
      INTEGER               ::   mstart(4), mcount(4)
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncsed_ini_g3dvar.f'
!
!----------------------------------------------------------------------
! check input type and set variable dimension
!----------------------------------------------------------------------
!
      IF(Iinfo(idvar).EQ.g3dvar) THEN 
        mcount=(/NX,NY,3,       1/)
        mstart=(/ 1, 1,1,iniTflag/)
      ELSE
        WRITE(OUT,*) TRIM(ADJUSTL(SUBNA))
        WRITE(OUT,*) 'ERROR: misuse of ncsed_ini_g3dvar.f ',
     .               TRIM(ADJUSTL(Vinfo(1,idvar))),
     .               ' is not 3D time-dependent variable'
        CALL EXIT
      ENDIF

! inquire variable
      status=nf90_inq_varid(ncID,TRIM(ADJUSTL(Vinfo(1,idvar))),varID)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

! read values
      status=nf90_get_var(ncID,varID,Vvalue,start=mstart,count=mcount)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

      END SUBROUTINE ncsed_ini_g3dvar
