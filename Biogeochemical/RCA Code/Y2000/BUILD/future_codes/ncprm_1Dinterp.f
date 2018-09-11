      SUBROUTINE ncprm_1Dinterp(ncID,idvar_time,idvar,Trec,
     .                          NXprmT,Sprm,Bprm)
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
      INTEGER, INTENT(in)   ::   ncID        ! identification of netcdf file
      INTEGER, INTENT(in)   ::   idvar_time  ! variable identification (see NetCDFCM)
      INTEGER, INTENT(in)   ::   idvar       ! variable identification (see NetCDFCM)
      INTEGER, INTENT(inout)::   Trec(2)     ! two time record number in NetCDF
      REAL,    INTENT(out)  ::   NXprmT      ! Next NetCDF time for parameter
      REAL,    INTENT(out)  ::   Sprm        ! slope estimation for parameter
      REAL,    INTENT(out)  ::   Bprm        ! parameter value at one instant before
!
!  Local variable declarations
!
      INTEGER               ::   totalREC
      INTEGER               ::   tstart(1), tcount(1)
      REAL                  ::   tt(2), vv(2), Vtem(1)
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncprm_1Dinterp.f'
!
!----------------------------------------------------------------------
! check input type and set variable dimension
!----------------------------------------------------------------------
!
      IF(Iinfo(idvar).EQ.t1dvar) THEN 
        tstart=(/Trec(1)/)
        tcount=(/      2/)
      ELSE
        WRITE(OUT,*) 'ERROR: misuse of ncprm_1Dinterp.f ',
     .               TRIM(ADJUSTL(Vinfo(1,idvar))),
     .               ' is not 1D time-independent variable'
        CALL EXIT
      ENDIF
!
!----------------------------------------------------------------------
! linear interpolation between NetCDF time record to RCA TIME
!----------------------------------------------------------------------
!
! check time record number and make sure t(1)<= TIME < t(2),
! otherwise, update Trec.
      status=nf90_inq_dimid(ncID,
     .       TRIM(ADJUSTL(Vinfo(1,idvar_time))),dimID)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      status=nf90_inquire_dimension(ncID,dimID,len=totalREC)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      status=nf90_inq_varid(ncID,Vinfo(1,idvar_time),varID)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      status=nf90_get_var(ncID,VarID,tt,start=tstart,count=tcount)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      DO WHILE(tt(2).LE.TIME)
        Trec(1)=Trec(2)
        Trec(2)=Trec(2)+1 
        tstart=(/Trec(1)/)
        IF(Trec(2).GT.totalREC) THEN
          WRITE(OUT,'(A,f12.2,2A)') 'ERROR: current Model TIME (',TIME,
     .                              ') exceed the time range of ',
     .                              'params input file'
          CALL EXIT
        ENDIF
        status=nf90_get_var(ncID,VarID,tt,start=tstart,count=tcount)
        CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      ENDDO

! linear interpolation
      status=nf90_get_var(ncID,VarID,tt,start=tstart,count=tcount)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      status=nf90_inq_varid(ncID,TRIM(ADJUSTL(Vinfo(1,idvar))),varID)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)
      status=nf90_get_var(ncID,varID,vv,start=tstart,count=tcount)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)
!
!                 tt  - tt1
!       vvout = ------------- x (vv2-vv1)+vv1
!                 tt2 - tt1
!
      NXprmT=tt(2)
      Sprm=(vv(2)-vv(1))/(tt(2)-tt(1))*(vv(2)-vv(1))
      Bprm=vv(1)

      END SUBROUTINE ncprm_1Dinterp
