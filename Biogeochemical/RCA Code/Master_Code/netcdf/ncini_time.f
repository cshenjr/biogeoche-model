      SUBROUTINE ncini_time(ncID)
!
!=======================================================================
!  Read RCA boundary variables in NetCDF format                        !
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
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncini_time.f'
!
!----------------------------------------------------------------------
! check input type and set variable dimension
!----------------------------------------------------------------------
!
      idvar=idini_time 
      IF(Iinfo(idvar).EQ.t1dvar) THEN
        mstart=(/iniTflag/)
        mcount=(/1/)
      ELSE
        WRITE(OUT,*) TRIM(ADJUSTL(SUBNA))
        WRITE(OUT,*) 'ERROR: misuse of ncini_time.f ',
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
      Vtwarp=TRIM(ADJUSTL(Vinfo(3,idini_time)))
      Vtscal=IUNITCHECK(Vtwarp,TRIM(ADJUSTL(Vinfo(3,idini_time))))
      TZERO=Vtscal*Vtem(1)/86400.
      WRITE(OUT,'(2A,f12.5,A)') TRIM(ADJUSTL(SUBNA)),
     .                      '   Model TZERO set to ', TZERO, ' DAYS'
      WRITE(OUT,100)
 100  FORMAT(/////1X,119('*')//)

      END SUBROUTINE ncini_time
