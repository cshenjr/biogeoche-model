      SUBROUTINE nccrt_rcasys(FNAME)
!
!=======================================================================
! Create a NetCDF file ready for RCA output                            !
!                            YUN LI, UMCES/HPL Sep-30-2008             !
!=======================================================================
!
      USE netcdf
      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
! 
!  Imported variable declarations.
!
      CHARACTER(*),INTENT(in)  ::     FNAME   
!
!  Local variable declarations
!
      INCLUDE 'ncdim_info.h'
      INTEGER       ::     ncID 
      REAL          ::     Xg(NX), Yg(NY), Zg(NZ), ZZg(NZ+1),Gg(group)
      INTEGER       ::     IX, IY, IZ, IZZ, IG
      CHARACTER(8)  ::     mydate 
      CHARACTER(10) ::     mytime
      CHARACTER(19) ::     Y4MMDDhhmmss
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL nccrt_rcasys.f'
!
!-----------------------------------------------------------------------
! set up dimensions
!-----------------------------------------------------------------------
!
       DO IX=1,NX
          Xg(IX)=IX*1.0
       ENDDO
       DO IY=1,NY
          Yg(IY)=IY*1.0
       ENDDO
       DO IZ=1,NZ
          Zg(IZ)=IZ*1.0
       ENDDO
       DO IZZ=1,NZ+1
          ZZg(IZZ)=IZZ*1.0
       ENDDO
       DO IG=1,group
          Gg(IG)=IG*1.0
       ENDDO
!
!-----------------------------------------------------------------------
! create netcdf file
!-----------------------------------------------------------------------
!
! Global infomation into netcdf 
      status=nf90_create(TRIM(ADJUSTL(FNAME)),nf90_clobber,ncID)
      CALL nccheck_status(status,TRIM(ADJUSTL(FNAME)),TRIM(SUBNA))
      CALL date_and_time(mydate,mytime)
      Y4MMDDhhmmss( 1:10)=      mydate(1:4)//  ! year
     .                     '-'//mydate(5:6)//  ! month
     .                     '-'//mydate(7:8)    ! day
      Y4MMDDhhmmss(11:19)=','//mytime(1:2)//   ! hour
     .                    ':'//mytime(3:4)//   ! minute
     .                    ':'//mytime(5:6)     ! second
      status=nf90_put_att(ncID,nf90_global,    ! version info
     .       'Version',TRIM(ADJUSTL(RCA_version)))
      status=nf90_put_att(ncID,nf90_global,    ! author
     .       'Author',author)
      status=nf90_put_att(ncID,nf90_global,    ! title/case
     .       'Title',RCA_application)
      status=nf90_put_att(ncID,nf90_global,    ! date
     .       'Date', Y4MMDDhhmmss )

! define dimensions
      SUBNA_h=TRIM(ADJUSTL(SUBNA))//' -> ncdef_dim.h'
      INCLUDE 'ncdef_dim.h'

! define variables
      CALL ncwrt_n1dvar(ncID,idNX , Xg, NX)
      CALL ncwrt_n1dvar(ncID,idNY , Yg, NY)
      CALL ncwrt_n1dvar(ncID,idNZ , Zg, NZ)
      CALL ncwrt_n1dvar(ncID,idNZZ,ZZg,NZ+1)
      CALL ncwrt_n1dvar(ncID,idNG , Gg,group)

! close netcdf file
      status=nf90_close(ncID)
      CALL nccheck_status(status,FNAME,SUBNA)

! write grid information
      status=nf90_open(TRIM(ADJUSTL(FNAME)),nf90_write,ncID)
      CALL nccheck_status(status,FNAME,SUBNA)
      CALL ncwrt_n2dvar(ncID,idFSM,FSM)
      CALL ncwrt_n2dvar(ncID,idH  ,H)
      CALL ncwrt_n2dvar(ncID,idDX ,DX)
      CALL ncwrt_n2dvar(ncID,idDY ,DY)
      CALL ncwrt_n2dvar(ncID,idLAT,LAT)
      CALL ncwrt_n2dvar(ncID,idLON,LON)
      status=nf90_close(ncID)
      CALL nccheck_status(status,FNAME,SUBNA)

      END SUBROUTINE nccrt_rcasys
