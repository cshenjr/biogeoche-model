      SUBROUTINE ncwrt_rstvar(ncID,idtime,idvars,Tvalue,Vvalue)
!
!=======================================================================
!  WRITE RCA variables into NetCDF restart file                        !
!                                 YUN LI, UMCES/HPL, Jun-04-2011       !
!=======================================================================
!
      USE netcdf
      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
! 
!  Imported variable declarations.
!
      INTEGER,          INTENT(in)   ::   ncID         ! identification of netcdf file
      INTEGER,          INTENT(in)   ::   idtime       ! variable identification (see NetCDFCM)
      INTEGER,          INTENT(in)   ::   idvars(NOSYS)! variable identification (see NetCDFCM)
      REAL,             INTENT(in)   ::   Tvalue(1)
      REAL,             INTENT(in)   ::   Vvalue(NX,NY,NZ,NOSYS)
!
!  Local variable declarations
!
      INCLUDE 'ncdim_info.h'
      INTEGER              ::   idvar
      INTEGER              ::   IX, IY, IZ
      INTEGER,ALLOCATABLE  ::   mstart(:), mcount(:), mdimen(:)
      REAL                 ::   Vvalue_mask(NX,NY,NZ)
      CHARACTER(100)       ::   coordinates
      LOGICAL              ::   create_var
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncwrt_rstvar.f'
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
! check input TIME and set variable dimension
!----------------------------------------------------------------------
!
      idvar=idtime
      ALLOCATE(mdimen(1))
      ALLOCATE(mstart(1))
      ALLOCATE(mcount(1))
! check time variable type
      IF(Iinfo(idvar).EQ.t1dvar) THEN
        coordinates=Vinfo(6,idvar)
        mdimen=(/TdimID/)
        mstart=(/rstTflag/)
        mcount=(/1/)
      ELSE
        WRITE(OUT,*) 'ERROR: misuse of ncwrt_t1dvar.f ',
     .               TRIM(ADJUSTL(Vinfo(1,idvar))),
     .               ' is not 1D time-dependent variable'
        CALL EXIT
      ENDIF

! define/inquire time variable
      IF(create_var) THEN
        SUBNA_h=TRIM(ADJUSTL(SUBNA))//' -> ncdef_var.h'
        INCLUDE 'ncdef_var.h'
      ENDIF
      status=nf90_inq_varid(ncID,
     .                TRIM(ADJUSTL(Vinfo(1,idvar))),varID)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

! put value
      status=nf90_put_var(ncID,varID,Tvalue,start=mstart,count=mcount)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)
      DEALLOCATE(mdimen)
      DEALLOCATE(mstart)
      DEALLOCATE(mcount)
!
!----------------------------------------------------------------------
! check input CARAY and set variable dimension
!----------------------------------------------------------------------
!
      ALLOCATE(mdimen(4))
      ALLOCATE(mstart(4))
      ALLOCATE(mcount(4))
      DO ISYS=1,NOSYS
! check system variable type
        idvar=idvars(ISYS)
        IF(Iinfo(idvar).EQ.g3dvar) THEN 
          coordinates='x_r y_r s_r  '//Vinfo(5,idvar)
          mdimen=(/NXdimID, NYdimID, NZdimID, TdimID/)
          mcount=(/NX,NY,NZ,1/)
          mstart=(/ 1, 1,1,rstTflag/)
        ELSE
          WRITE(OUT,*) TRIM(ADJUSTL(SUBNA_h))
          WRITE(OUT,*) 'ERROR: misuse of ncwrt_rstvar.f ',
     .                 TRIM(ADJUSTL(Vinfo(1,idvar))),
     .                ' is not 3D time-dependent variable'
          CALL EXIT
        ENDIF

! define/inquire variable
        IF(create_var) THEN
          SUBNA_h=TRIM(ADJUSTL(SUBNA))//' -> ncdef_var.h'
          INCLUDE 'ncdef_var.h'
        ENDIF
        status=nf90_inq_varid(ncID,
     .                 TRIM(ADJUSTL(Vinfo(1,idvar))),varID)
        CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

! set land values to default spval (for ncview display)
        DO IX=1,NX
        DO IY=1,NY
          IF(FSM(IX,IY).EQ.0.0) THEN
             Vvalue_mask(IX,IY,:)=spval
          ELSE
            DO IZ=1,NZ
             Vvalue_mask(IX,IY,IZ)=Vvalue(IX,IY,IZ,ISYS)
            ENDDO
          ENDIF
        ENDDO
        ENDDO

! put values
        status=nf90_put_var(ncID,varID,Vvalue_mask,
     .                      start=mstart,
     .                      count=mcount)
        CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)
      END DO
      DEALLOCATE(mdimen)
      DEALLOCATE(mstart)
      DEALLOCATE(mcount)
!
!----------------------------------------------------------------------
! Change rstTflag back to zero 
!----------------------------------------------------------------------
!
      IF(create_var) rstTflag=0

      END SUBROUTINE ncwrt_rstvar
