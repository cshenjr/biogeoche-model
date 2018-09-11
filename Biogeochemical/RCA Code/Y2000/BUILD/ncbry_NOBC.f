      SUBROUTINE ncbry_NOBC
!
!=======================================================================
!  Calculate number of boundary cells in RCA                           ! 
!                                 YUN LI, UMCES/HPL, Feb-22-2011       !
!=======================================================================
!
      USE netcdf
      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
!
!  Local variable declarations
!
      INTEGER               ::   IX, IY 
      INTEGER               ::   num_obc 
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncbry_NOBC.f'
!
!-----------------------------------------------------------------------
! calculate number of boundary cells
!-----------------------------------------------------------------------
!
      ! riverine boundary
      num_obc=nriver*NZ

      ! lateral boundary
      DO IY=2,NY-1
        IF(FSM(NX-1,IY).EQ.-2.0) num_obc=num_obc+NZ  ! east
        IF(FSM(   2,IY).EQ.-2.0) num_obc=num_obc+NZ  ! west
      ENDDO
      DO IX=2,NX-1
        IF(FSM(IX,NY-1).EQ.-2.0) num_obc=num_obc+NZ  ! north
        IF(FSM(IX,   2).EQ.-2.0) num_obc=num_obc+NZ  ! south
      ENDDO
      NOBC(ISYS)=num_obc

      IF(SYSBY(ISYS).EQ.1) NOBC(ISYS)=0

      END SUBROUTINE ncbry_NOBC
