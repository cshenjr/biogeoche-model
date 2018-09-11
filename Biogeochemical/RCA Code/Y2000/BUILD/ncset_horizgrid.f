       SUBROUTINE ncset_horizgrid
!
!======================================================================
! Read the horizontal grid information that required by RCA           !
!    H  : Mean water depth (m)                                        !
!    DX : Cell length or a distance across a cell in XI-direction     !
!             through its center (m)                                  !
!    DY : Cell width or a distance across a cell in ETA-direction     !
!             through its center (m)                                  !
!    FSM: Land/water mask at the center of the cell                   !
!                                      YUN LI, UMCES/HPL Sep-30-2008  !
!======================================================================
!
      USE netcdf
!
      SAVE
      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
!
! Variables in ROMS NetCDF file 
!
      REAL*8      ::  h_roms(NX-2,NY-2)
      REAL*8      ::  pm(NX-2,NY-2)
      REAL*8      ::  pn(NX-2,NY-2)
      REAL*8      ::  maskr(NX-2,NY-2)
      REAL*8      ::  latr(NX-2,NY-2)
      REAL*8      ::  lonr(NX-2,NY-2)
      REAL*8      ::  river_transport(nriver)
      REAL*8      ::  river_Xposition(nriver)
      REAL*8      ::  river_Eposition(nriver)
      REAL*8      ::  river_direction(nriver)
!
!  Local variable declarations.
!
      INTEGER     ::  ncID
      INTEGER     ::  IX, IY, IR, rI, rJ
      INTEGER     ::  I_roms, J_roms, num_river
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncset_horizgrid.f'
!
!-----------------------------------------------------------------------
! Get ROMS horizontal grid information
!-----------------------------------------------------------------------
!
! Get h, pm, pn from ROMS grid file                              
!
      WRITE(OUT,'(1x,a,t30,a)'),
     .      'Horizontal GRID FILE NAME ', TRIM(ADJUSTL(GRDNA))
      status=nf90_open(TRIM(ADJUSTL(GRDNA)),nf90_nowrite,ncID)
      CALL nccheck_status(status,GRDNA,SUBNA)
      status=nf90_inq_varid(ncID,'pm',varID)       ! 1/dx
      status=nf90_get_var(ncID,varID,pm)
      CALL nccheck_status(status,'pm',SUBNA)
      status=nf90_inq_varid(ncID,'pn',varID)       ! 1/dy
      status=nf90_get_var(ncID,varID,pn)
      CALL nccheck_status(status,'pn',SUBNA)
      status=nf90_inq_varid(ncID,'h',varID)        ! Bathymetry
      status=nf90_get_var(ncID,varID,h_roms)
      CALL nccheck_status(status,'h_roms',SUBNA)
      status=nf90_inq_varid(ncID,'mask_rho',varID) ! Mask
      status=nf90_get_var(ncID,varID,maskr)
      CALL nccheck_status(status,'mask_rho',SUBNA)
      status=nf90_inq_varid(ncID,'lat_rho',varID)  ! latitude
      status=nf90_get_var(ncID,varID,latr)
      CALL nccheck_status(status,'lat_rho',SUBNA)
      status=nf90_inq_varid(ncID,'lon_rho',varID)  ! longitude
      status=nf90_get_var(ncID,varID,lonr)
      CALL nccheck_status(status,'lon_rho',SUBNA)
      status=nf90_close(ncID)
      CALL nccheck_status(status,GRDNA,SUBNA)
!
! Get river position indices from ROMS forcing file
!
      WRITE(OUT,'(1x,a,t30,a)'),
     .      '          RIVER FILE NAME ', TRIM(ADJUSTL(RVRNA))
      status=nf90_open(TRIM(ADJUSTL(RVRNA)),nf90_nowrite,ncID)
      CALL nccheck_status(status,RVRNA,SUBNA)
      status=nf90_inq_dimid(ncID,'river',dimID)
      CALL nccheck_status(status,'river dimID',SUBNA)
      status=nf90_inquire_dimension(ncID,dimID,len=num_river)
      CALL nccheck_status(status,'river var. length',SUBNA)
      IF(num_river.NE.nriver) THEN
        WRITE(OUT,*) 'ERROR: in ', TRIM(ADJUSTL(SUBNA)),
     .               ' number of rivers does not match 
     .                 nriver in NetCDFCM' 
      ENDIF
      status=nf90_inq_varid(ncID,'river_Xposition',varID)
      status=nf90_get_var(ncID,varID,river_Xposition)
      CALL nccheck_status(status,'river_Xposition',SUBNA)
      status=nf90_inq_varid(ncID,'river_Eposition',varID)
      status=nf90_get_var(ncID,varID,river_Eposition)
      CALL nccheck_status(status,'river_Eposition',SUBNA)
      status=nf90_inq_varid(ncID,'river_direction',varID)
      status=nf90_get_var(ncID,varID,river_direction)
      CALL nccheck_status(status,'river_direction',SUBNA)
      status=nf90_inq_varid(ncID,'river_transport',varID)
      status=nf90_get_var(ncID,varID,river_transport,
     .                    start = (/     1, 1/),
     .                    count = (/nriver, 1/))
      CALL nccheck_status(status,'river_transport',SUBNA)
      status=nf90_close(ncID)
      CALL nccheck_status(status,RVRNA,SUBNA)
!
!-----------------------------------------------------------------------
! Calcualte the varible to what RCA required
!    H   :  mean water depth    [Unit:m]
!    DX  :  length of each cell [Unit:m]
!    DY  :  width of each cell  [Unit:m]
!    FSM :  land mask array     [0=land 1=water -1=river_bc -2=open_bc]
! NOTE:
! due to adding land mask aroumd ROMS domain
!    I_roms=IX-1
!    J_roms=IY-1
!-----------------------------------------------------------------------
!
! land point around domain, and mask is equal to 0.0
!
      DX = 0.0
      DY = 0.0
      H = 1.0E-10
      FSM = 0.0
      LAT = 0.0
      LON = 0.0
      DO IY=2,NY-1
      DO IX=2,NX-1
         I_roms=IX-1
         J_roms=IY-1
         DY(IX,IY)=1.0/pn(I_roms,J_roms)
         DX(IX,IY)=1.0/pm(I_roms,J_roms)
         H(IX,IY)=h_roms(I_roms,J_roms)*maskr(I_roms,J_roms)
         FSM(IX,IY)=maskr(I_roms,J_roms)
         LAT(IX,IY)=latr(I_roms,J_roms)
         LON(IX,IY)=lonr(I_roms,J_roms)
         IF(FSM(IX,IY).EQ.0.0) H(IX,IY)=1.0E-10
      ENDDO
      ENDDO
!
! mask of river boundary cell is -1.0
!
      DO IR=1,nriver
         rI = int(river_Xposition(IR))+1
         rJ = int(river_Eposition(IR))+1
         IF(river_direction(IR).EQ.0.0) rJ=rJ+1
         IF(river_direction(IR).EQ.1.0) rI=rI+1
         IF((river_direction(IR).EQ.0.0) .AND. 
     .      (river_transport(IR).GT.0.0)) rI=rI+1
         IF((river_direction(IR).EQ.1.0) .AND. 
     .      (river_transport(IR).GT.0.0)) rJ=rJ+1
         FSM(rI,rJ)=-1.0
         riverX(IR)=rI
         riverY(IR)=rJ
      ENDDO
!
! mask of open boundary is -2.0
!
      DO IY=2,NY-1
      DO IX=2,NX-1
         I_roms=IX-1
         J_roms=IY-1
         IF((IX.EQ. 2 .OR. IX.EQ.NX-1 .OR. 
     .       IY.EQ. 2 .OR. IY.EQ.NY-1)) THEN
           IF(maskr(I_roms,J_roms).EQ.1.0) FSM(IX,IY)=-2.0
         ENDIF
      ENDDO
      ENDDO

      END SUBROUTINE ncset_horizgrid
