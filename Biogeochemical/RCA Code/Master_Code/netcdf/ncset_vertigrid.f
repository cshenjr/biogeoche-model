      SUBROUTINE ncset_vertigrid
!
!======================================================================
! Read the vertical grid information required by RCA                  !
!    DZ:  thickness of the sigma level (fractional #, nondimensional) !
!    DZZ: avarage depth of the grid element (same as above)           !
!                                      YUN LI, UMCES/HPL Feb-23-2009  !
!======================================================================
!
      USE netcdf
      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
!
!  Variables in ROMS NetCDF file 
!
      REAL*8            ::  Cs_r(NZ)
      REAL*8            ::  Cs_w(NZ+1)
!
!  Local variable declarations.
!
      INTEGER           ::  ncID, IZ
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncset_vertigrid.f'
!
!----------------------------------------------------------------------      
! Get ROMS vertical grid information 
!----------------------------------------------------------------------      
!
! Get Cs_w, Cs_r from ROMS history file
!
      WRITE(OUT,'(1x,a,t30,a)'),
     .      '  Vertical GRID FILE NAME ', TRIM(ADJUSTL(HISNA))
      status=nf90_open(TRIM(ADJUSTL(HISNA)),nf90_nowrite,ncID)
      CALL nccheck_status(status,HISNA,SUBNA)
      status=nf90_inq_varid(ncID,'Cs_w',varID) ! S-coordinate at W-points ROMS3.0
      status=nf90_get_var(ncID,varID,Cs_w)
      CALL nccheck_status(status,'Cs_w',SUBNA)
      status=nf90_inq_varid(ncID,'Cs_r',varID) ! S-coordinate at RHO-points ROMS3.0
      status=nf90_get_var(ncID, varID, Cs_r)
      CALL nccheck_status(status,'Cs_r',SUBNA)
      status=nf90_close(ncID)
      CALL nccheck_status(status,HISNA,SUBN)
!
!----------------------------------------------------------------------      
! transfer ROMS S-coord to RCA Sigma-coord
!  NOTE: 1) layer indices reversed from ROMS to RCA
!        2) ROMS does not include sediment layer
!           without sediment, DZ(NZ+1)=DZZ(NZ)=DZZ(NZ+1)=0
!
!                             ROMS       RCA
!  SEA SURFACE   Cs_w(21)  ---------  ---------  --|--  
!                Cs_r(20)      20         1       DZ(1)     --|--
!                Cs_w(20)  ---------  ---------  --|--      DZZ(1)
!                Cs_r(19)      19         2       DZ(2)     --|--
!                Cs_w(19)  ---------  ---------  --|--
!                            ....        ....
!                Cs_w(3)   ---------  ---------  --|--
!                Cs_r(2)        2       NZ-1      DZ(NZ-1)  --|--
!                Cs_w(2)   ---------  ---------  --|--   DZZ(NZ-1)
!                Cs_r(1)        1        NZ       DZ(NZ)    --|--
!  SEA-BED       Cs_w(1)   ---------  ---------  -----   DZZ(NZ:NZ+1) = 0.0
!    SEDIMENT                 N/A       NZ+1      DZ(NZ+1)=0.0  
!
!----------------------------------------------------------------------      
!
      DZ=0.0
      DZZ=0.0
      DO IZ=1,NZ
         DZ(IZ)=Cs_w(NZ+2-IZ)-Cs_w(NZ+1-IZ)     ! => DZ(NZ+1)==0.0
      ENDDO
      DO IZ=1,NZ-1
         DZZ(IZ)=Cs_r(NZ+1-IZ)-Cs_r(NZ-IZ)      ! => DZZ(NZ:NZ+1)=0.0
      ENDDO
      
      END SUBROUTINE ncset_vertigrid
