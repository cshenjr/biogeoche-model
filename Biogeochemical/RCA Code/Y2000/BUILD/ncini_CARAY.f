      SUBROUTINE ncini_CARAY(ncID)
!
!=======================================================================
!  Read RCA initial tracer concentrations in NetCDF format             ! 
!                                 YUN LI, UMCES/HPL, Feb-22-2011       !
!=======================================================================
!
      USE netcdf
      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
! 
!  Imported variable declarations.
!
      INTEGER, INTENT(in)   ::   ncID        ! identification of netcdf file
!
!  Local variable declarations
!
      INTEGER               ::   idvar
      REAL                  ::   Vscale
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncini_CARAY.f'
!
!----------------------------------------------------------------------
! check input type
!----------------------------------------------------------------------
!
      DO ISYS=1,NOSYS
        WRITE(OUT,1400)   ISYS
 1400   FORMAT(/44X,'INITIAL CONDITIONS FOR SYSTEM',I3)
        WRITE(OUT,2150)   SYNAME(ISYS)
 2150   FORMAT(41X,11('-'),3X,A8,3X,11('-')/)

        idvar=idTvar(ISYS)
        Vvalue=Fscale(idvar)
        IF(Iinfo(idvar).EQ.g3dvar) THEN
        ELSE
          WRITE(OUT,*) TRIM(ADJUSTL(SUBNA))
          WRITE(OUT,*) 'ERROR: misuse of ncini_CARAY.f ',
     .                 TRIM(ADJUSTL(Vinfo(1,idvar))),
     .                 ' is not 3D gridded variable'
          CALL EXIT
        ENDIF
!
!-----------------------------------------------------------------------
! load initial tracer concentrations
!-----------------------------------------------------------------------
!
        IF(SYSBY(ISYS).EQ.0) THEN
          status=nf90_inq_varid(ncID,TRIM(ADJUSTL(Vinfo(1,idvar))),
     .                          VarID)
          status=nf90_get_var(ncID,VarID,CARAY(:,:,:,ISYS),
     .                        start=(/1,1,1,iniTflag/),
     .                        count=(/NX,NY,NZ,1/))
          CALL nccheck_status(status,Vinfo(1,idvar),SUBNA)
        ELSE
          CARAY(:,:,:,ISYS)=0.0
        ENDIF
      END DO

      END SUBROUTINE ncini_CARAY
