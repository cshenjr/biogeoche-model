      SUBROUTINE ncsed_rst_t1dvar(ncID,idvar,Vvalue)
!
!=======================================================================
!  WRITE RCA variables into NetCDF output                              !
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
      INTEGER,          INTENT(in)   ::   ncID        ! identification of netcdf file
      INTEGER,          INTENT(in)   ::   idvar       ! variable identification (see NetCDFCM)
      REAL,             INTENT(in)   ::   Vvalue(1)
!
!  Local variable declarations
!
      INCLUDE 'ncdim_info.h'
      INTEGER              ::   mstart(1), mcount(1), mdimen(1)
      CHARACTER(100)       ::   coordinates
      LOGICAL              ::   create_var
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncsed_rst_t1dvar.f'
!
!----------------------------------------------------------------------
! inqure dimension information
!----------------------------------------------------------------------
!
      SUBNA_h=TRIM(ADJUSTL(SUBNA))//' -> ncinq_dim.h'
      INCLUDE 'ncinq_dim.h'
!
!----------------------------------------------------------------------
! turn on create_var switch when rstTflag=0 and temporary rstTflag=1
!----------------------------------------------------------------------
!
      create_var=.FALSE.
      IF(rstTflag.EQ.0) THEN
        create_var=.TRUE.  !! create variable of restart file
        rstTflag=1         !! write value to first record
      ENDIF
!
!----------------------------------------------------------------------
! check input type and set variable dimension
!----------------------------------------------------------------------
!
      IF(Iinfo(idvar).EQ.t1dvar) THEN 
        coordinates=Vinfo(5,idvar)
        mdimen=(/TdimID/)
        mstart=(/rstTflag/)
        mcount=(/1/)
      ELSE
        WRITE(OUT,*) TRIM(ADJUSTL(SUBNA))
        WRITE(OUT,*) 'ERROR: misuse of ncsed_rst_t1dvar.f ',
     .               TRIM(ADJUSTL(Vinfo(1,idvar))),
     .               ' is not 1D time-dependent variable'
        CALL EXIT
      ENDIF

! define/inquire variable
      IF(create_var) THEN
        SUBNA_h=TRIM(ADJUSTL(SUBNA))//' -> ncdef_var.h'
        INCLUDE 'ncdef_var.h'
      ENDIF
      status=nf90_inq_varid(ncID,TRIM(ADJUSTL(Vinfo(1,idvar))),varID)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

! put values
      status=nf90_put_var(ncID,varID,Vvalue,start=mstart,count=mcount)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)
!
!----------------------------------------------------------------------
! Change rstTflag back to zero 
!----------------------------------------------------------------------
!
      IF(create_var) rstTflag=0
      
      END SUBROUTINE ncsed_rst_t1dvar

