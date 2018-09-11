      SUBROUTINE ncsed_rst_g3dvar(ncID,idvar,Vvalue)
!
!=======================================================================
!  WRITE RCA sediment variables into NetCDF restart file               !
!                                 YUN LI, UMCES/HPL, Jun-05-2011       !
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
      REAL,             INTENT(in)   ::   Vvalue(NX,NY,NZ)
!
!  Local variable declarations
!
      INCLUDE 'ncdim_info.h'
      INTEGER              ::   IX, IY, IZ
      INTEGER              ::   mstart(4), mcount(4), mdimen(4)
      REAL                 ::   Vvalue_mask(NX,NY,NZ)
      CHARACTER(100)       ::   coordinates
      LOGICAL              ::   create_var
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncsed_rst_g3dvar.f'
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
      IF(Iinfo(idvar).EQ.g3dvar) THEN 
        coordinates='x_r y_r group '//Vinfo(5,idvar)
        mdimen=(/NXdimID, NYdimID, NGdimID, TdimID/)
        mcount=(/NX,NY,group,        1/)
        mstart=(/ 1, 1,    1, rstTflag/)
      ELSE
        WRITE(OUT,*) TRIM(ADJUSTL(SUBNA))
        WRITE(OUT,*) 'ERROR: misuse of ncsed_rst_g3dvar.f ',
     .               TRIM(ADJUSTL(Vinfo(1,idvar))),
     .               ' is not 3D time-dependent variable'
        CALL EXIT
      ENDIF

! define/inquire variable
      IF(create_var) THEN
        SUBNA_h=TRIM(ADJUSTL(SUBNA))//' -> ncdef_var.h'
        INCLUDE 'ncdef_var.h'
      ENDIF
      status=nf90_inq_varid(ncID,TRIM(ADJUSTL(Vinfo(1,idvar))),varID)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)

! set land values to default spval (for ncview display)
      DO IX=1,NX
      DO IY=1,NY
        IF(FSM(IX,IY).EQ.0.0) THEN
           Vvalue_mask(IX,IY,:)=spval
        ELSE
          DO IZ=1,NZ
           Vvalue_mask(IX,IY,IZ)=Vvalue(IX,IY,IZ)
          ENDDO
        ENDIF
      ENDDO
      ENDDO

! put values
      status=nf90_put_var(ncID,varID,Vvalue_mask,
     .                    start=mstart,
     .                    count=mcount)
      CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)
!
!----------------------------------------------------------------------
! Change rstTflag back to zero 
!----------------------------------------------------------------------
!
      IF(create_var) rstTflag=0

      END SUBROUTINE ncsed_rst_g3dvar
