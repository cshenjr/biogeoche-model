      SUBROUTINE ncprm_2Dinterp(ncID,idvar_time,idvar,Trec,Vvalue)
!
!=======================================================================
!  Read 2D time-dependent parameters from NetCDF file                  !
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
      REAL,    INTENT(out)  ::   Vvalue(NX,NY)
!
!  Local variable declarations
!
      INTEGER               ::   mstart(3), mcount(3)
      INTEGER               ::   tstart(1), tcount(1)
      REAL                  ::   tt(2), vv(NX,NY,2)
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncprm_2Dinterp.f'
!
!----------------------------------------------------------------------
! check input type and set variable dimension
!----------------------------------------------------------------------
!
      IF(Iinfo(idvar).EQ.g2dvar) THEN 
        mcount=(/NX,NY,     2 /)
        mstart=(/ 1, 1,Trec(1)/)
        tstart=(/Trec(1)/)
        tcount=(/      2/)
      ELSE
        WRITE(OUT,*) 'ERROR: misuse of ncprm_2Dinterp.f ',
     .               TRIM(ADJUSTL(Vinfo(1,idvar))),
     .               ' is not 2D time-independent variable'
        CALL EXIT
      ENDIF
!
!----------------------------------------------------------------------
! linear interpolation between NetCDF time record to RCA TIME
!----------------------------------------------------------------------
!
! check time record number and make sure t(1)<= TIME < t(2),
! otherwise, update Trec.
      status=nf90_inq_varid(ncID,Vinfo(1,idvar_time),varID)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      status=nf90_get_var(ncID,VarID,tt,start=tstart,count=tcount)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      DO WHILE(tt(2).LE.TIME)
        Trec(1)=Trec(2)
        Trec(2)=Trec(2)+1 
        mstart=(/ 1, 1,Trec(1)/)
        tstart=(/Trec(1)/)
        status=nf90_get_var(ncID,VarID,tt,start=tstart,count=tcount)
        CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      ENDDO

! linear interpolation
      status=nf90_get_var(ncID,VarID,tt,start=tstart,count=tcount)
      CALL nccheck_status(status,Vinfo(1,idvar_time),SUBNA)
      status=nf90_inq_varid(ncID,TRIM(ADJUSTL(Vinfo(1,idvar))),varID)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)
      status=nf90_get_var(ncID,varID,vv,start=mstart,count=mcount)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)
!
!                 tt  - tt1
!       vvout = ------------- x (vv2-vv1)+vv1
!                 tt2 - tt1
!
      DO IX=1,NX
      DO IY=1,NY
        Vvalue(IX,IY)=(TIME-tt(1))/(tt(2)-tt(1))
     .               *(vv(IX,IY,2)-vv(IX,IY,1)) +vv(IX,IY,1)
      ENDDO
      ENDDO

      END SUBROUTINE ncprm_2Dinterp
