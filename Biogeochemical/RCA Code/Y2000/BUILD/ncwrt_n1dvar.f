      SUBROUTINE ncwrt_n1dvar(ncID,idvar,Vvalue,Vdim)
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
      INTEGER,          INTENT(in)   ::   Vdim        ! dimension of variable
      REAL,             INTENT(in)   ::   Vvalue(Vdim)
!
!  Local variable declarations
!
      INCLUDE 'ncdim_info.h'
      INTEGER              ::   mstart(1), mcount(1), mdimen(1)
      CHARACTER(100)       ::   coordinates
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncwrt_n1dvar.f'
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
      IF(Iinfo(idvar).EQ.n1dvar) THEN
      ELSE
        WRITE(OUT,*) TRIM(ADJUSTL(SUBNA_h))
        WRITE(OUT,*) 'ERROR: misuse of ncwrt_n1dvar.f ',
     .               TRIM(ADJUSTL(Vinfo(1,idvar))),
     .               ' is not 1D time-independent variable'
        CALL EXIT
      ENDIF

      mstart=(/1/)
      SELECT CASE (TRIM(ADJUSTL(Vinfo(1,idvar))))
        CASE('x_r')
         mcount=(/NX/)
         mdimen=(/NXdimID/)
         coordinates='NX'
        CASE('y_r')
         mcount=(/NY/)
         mdimen=(/NYdimID/)
         coordinates='NY'
        CASE('s_r')
         mcount=(/NZ/)
         mdimen=(/NZdimID/)
         coordinates='NZ'
        CASE('s_w')
         mcount=(/NZ+1/)
         mdimen=(/NZZdimID/)
         coordinates='NZZ'
        CASE('group')
         mcount=(/group/)
         mdimen=(/NGdimID/)
         coordinates='NG'
        CASE('system')
         mcount=(/NSYS/)
         mdimen=(/NSYSdimID/)
         coordinates='NSYS'
      END SELECT

! define/inquire variable
      SUBNA_h=TRIM(ADJUSTL(SUBNA))//' -> ncdef_var.h'
      INCLUDE 'ncdef_var.h'
      status=nf90_inq_varid(ncID,TRIM(ADJUSTL(Vinfo(1,idvar))),varID)
      CALL nccheck_status(status,TRIM(ADJUSTL(Vinfo(1,idvar))),SUBNA)

! put values
      status=nf90_put_var(ncID,varID,Vvalue,start=mstart,count=mcount)
      CALL nccheck_status(status,TRIM(ADJUSTL(Vinfo(1,idvar))),SUBNA)

      END SUBROUTINE ncwrt_n1dvar
