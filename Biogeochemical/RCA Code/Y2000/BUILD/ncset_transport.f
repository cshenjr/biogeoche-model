      SUBROUTINE ncset_transport(ncID)
!
!=======================================================================
! READ ROMS HYDRODYNAMIC OUTPUT AND APPLY TRANSFERMATION               !
! NOTE:                                                                !
!      vertical indice reversed                                        !
!         NZ indices      ROMS       RCA                               !
!        ------------   --------  --------                             !
!              1         bottom    surface                             ! 
!              2           |         /|\                               !
!            .....         |          |                                !
!              NZ-1       \|/         |                                !
!              NZ        surface    bottom                             !
!                                                                      !
!                                      Yun Li, UMCES/HPL, Sep-23-2008  !
!=======================================================================
!
      USE netcdf
!
      SAVE
      INCLUDE 'RCACM'
      INCLUDE 'NetCDFCM'
!
!  Imported and exported variable declarations. 
!
      INTEGER,      INTENT(in)   ::    ncID
!
!  Variables in ROMS NetCDF file
!
      REAL*8     ::  pm(NX-2,NY-2)
      REAL*8     ::  pn(NX-2,NY-2)
      REAL*8     ::  h_roms(NX-2,NY-2)
      REAL*8     ::  masku(NX-3,NY-2)
      REAL*8     ::  maskv(NX-2,NY-3)
      REAL*8     ::  maskr(NX-2,NY-2)
      REAL*8     ::  ocean_time(1)
      REAL*8     ::  Cs_w(NZ+1)
      REAL*8     ::  zeta(NX-2,NY-2)
      REAL*8     ::  salt(NX-2,NY-2,NZ)
      REAL*8     ::  temp(NX-2,NY-2,NZ)
      REAL*8     ::  u(NX-3,NY-2,NZ)
      REAL*8     ::  v(NX-2,NY-3,NZ)
      REAL*8     ::  omega(NX-2,NY-2,NZ+1)
      REAL*8     ::  AKs(NX-2,NY-2,NZ+1)
!
!  Local variable declarations
!
      INTEGER        ::  IX, IY, IZ
      INTEGER        ::  I_roms, J_roms
      REAL           ::  D(NX,NY)
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncset_transport.f'
!
!-----------------------------------------------------------------------
! Get ROMS input
!-----------------------------------------------------------------------
!
! ROMS grid file
      status=nf90_open(TRIM(ADJUSTL(GRDNA)),nf90_nowrite,ncID_grd)
      CALL nccheck_status(status,GRDNA,SUBNA)
      status=nf90_inq_varid(ncID_grd,'pm',varID)       ! 1/dx 
      status=nf90_get_var(ncID_grd,varID,pm)
      CALL nccheck_status(status,'pm',SUBNA)
      status=nf90_inq_varid(ncID_grd,'pn',varID)       ! 1/dy
      status=nf90_get_var(ncID_grd,varID,pn)
      CALL nccheck_status(status,'pn',SUBNA)
      status=nf90_inq_varid(ncID_grd,'h',varID)        ! Bathymetry
      status=nf90_get_var(ncID_grd,varID,h_roms)
      CALL nccheck_status(status,'h',SUBNA)
      status=nf90_inq_varid(ncID_grd,'mask_u',varID)   ! u-point mask
      status=nf90_get_var(ncID_grd,varID,masku)
      CALL nccheck_status(status,'mask_u',SUBNA)
      status=nf90_inq_varid(ncID_grd,'mask_v',varID)   ! v-point mask
      status=nf90_get_var(ncID_grd,varID,maskv)
      CALL nccheck_status(status,'mask_v',SUBNA)
      status=nf90_inq_varid(ncID_grd,'mask_rho',varID) ! r-point mask
      status=nf90_get_var(ncID_grd,varID,maskr)
      CALL nccheck_status(status,'mask_rho',SUBNA)
      status=nf90_close(ncID_grd)
      CALL nccheck_status(status,GRDNA,SUBNA)

! Get u, v, zeta, omega, Cs_w, AKs from ROMS history file
      status=nf90_inq_varid(ncID,'Cs_w',varID)     ! S-level, ROMS3.0
      status=nf90_get_var(ncID, varID, Cs_w)
      CALL nccheck_status(status,'Cs_w',SUBNA)
      status=nf90_inq_varid(ncID,'ocean_time',varID)
      status=nf90_get_var(ncID,varID,ocean_time,
     .                    start=(/hydTflag/),
     .                    count=(/1/))
      CALL nccheck_status(status,'ocean_time',SUBNA)
      status=nf90_inq_varid(ncID,'zeta',varID)     ! free surface
      status=nf90_get_var(ncID, varID, zeta,
     .                    start=(/1,1,hydTflag/),
     .                    count=(/NX-2,NY-2,1/))
      CALL nccheck_status(status,'zeta',SUBNA)
      status=nf90_inq_varid(ncID,'salt',varID)     ! salinity
      status=nf90_get_var(ncID, varID, salt,
     .                    start=(/1,1,1,hydTflag/),
     .                    count=(/NX-2,NY-2,NZ,1/))
      CALL nccheck_status(status,'salt',SUBNA)
      status=nf90_inq_varid(ncID,'temp',varID)     ! temperature
      status=nf90_get_var(ncID,varID, temp,
     .                    start=(/1,1,1,hydTflag/),
     .                    count=(/NX-2,NY-2,NZ,1/))
      CALL nccheck_status(status,'temp',SUBNA)
      status=nf90_inq_varid(ncID,'u',varID)        ! u-momentum
      status=nf90_get_var(ncID,varID,u, 
     .                    start = (/1,1,1,hydTflag/),
     .                    count = (/NX-3,NY-2,NZ,1/))
      CALL nccheck_status(status,'u',SUBNA)
      status=nf90_inq_varid(ncID,'v',varID)        ! v-momentum
      status=nf90_get_var(ncID,varID,v,
     .                    start = (/1,1,1,hydTflag/),
     .                    count = (/NX-2,NY-3,NZ,1/))
      CALL nccheck_status(status,'v',SUBNA)
      status=nf90_inq_varid(ncID,'omega',varID)    ! omega-momentum
      status=nf90_get_var(ncID, varID, omega,
     .                    start = (/1,1,1,hydTflag/),
     .                    count = (/NX-2,NY-2,NZ+1,1/))
      CALL nccheck_status(status,'omega',SUBNA)
      status=nf90_inq_varid(ncID,'AKs',varID)      ! diffusivity
      status=nf90_get_var(ncID, varID, AKs,
     .                    start = (/1,1,1,hydTflag/),
     .                    count = (/NX-2,NY-2,NZ+1,1/))
      CALL nccheck_status(status,'AKs',SUBNA)
!
!-----------------------------------------------------------------------
!  Calculate the variable RCA required
!-----------------------------------------------------------------------
!
!    QX,QY,QZ  :   low pass filtered volume flow rate in the XI-1,
!                  XI-2,Vertical direction (m3/sec) 	 
!    EX,EY,EZ  :   low pass filtered horizontal eddy diffusivity in 
!                  XI-1,XI-2,vertical direction (m2/sec)
!    ETA       :   initial surface elevation (m)
!    DETA      :   time rate of change of elevation (m/sec)
!    HYDSAL    :   low pass filtered salinity (psu)
!    HYDTEMP   :   low pass filtered temperature (Deg Celsius)
!
       QX=0.0
       QY=0.0
       QZ=0.0
       EX=0.0
       EY=0.0
       EZ=0.0
       ETA=0.0
       DETA=0.0
       HYDSAL=0.0
       HYDTEMP=0.0

! Get cell information      
      DZ=0.0
      DO IZ=1,NZ
         DZ(IZ)=Cs_w(NZ+2-IZ)-Cs_w(NZ+1-IZ)   ! reversed index in ROMS
      ENDDO
      DO IY=2,NY+1
      DO IX=2,NX+1
         D(IX,IY)=(h_roms(IX,IY)+zeta(IX,IY)) ! Total depth (h+zeta)
     .            *maskr(IX,IY)
      ENDDO
      ENDDO

! Calculate the hydrodynamic variables RCA required 
! NOTE:
! due to adding land mask aroumd ROMS domain
!     I_roms=IX-1
!     J_roms=IY-1
!
      DO IY=2,NY-1
      DO IX=2,NX-1
         I_roms=IX-1 
         J_roms=IY-1

         ! Calculate the surface elevation (m)
         ETA(IX,IY)=zeta(I_roms,J_roms)
     .            *maskr(I_roms,J_roms)

         DO IZ=1,NZ
           ! salinity (psu)
           HYDSAL(IX,IY,IZ)=salt(I_roms,J_roms,NZ+1-IZ)
     .                    *maskr(I_roms,J_roms)
           ! temperature
           HYDTEMP(IX,IY,IZ)=temp(I_roms,J_roms,NZ+1-IZ)
     .                     *maskr(I_roms,J_roms)
           ! advective fluxs across cell faces (m3/sec)
           IF(IX.EQ.2) THEN
              QX(IX,IY,IZ)=0.0
           ELSE
              QX(IX,IY,IZ)
     .       =2.0/(pn(I_roms,J_roms)+pn(I_roms-1,J_roms))
     .       *0.5*( D(I_roms,J_roms)+ D(I_roms-1,J_roms))
     .       *DZ(IZ)*u(I_roms-1,J_roms,NZ+1-IZ)
     .       *masku(I_roms-1,J_roms)
           ENDIF
           IF(IY.EQ.2) THEN
              QY(IX,IY,IZ)=0.0
           ELSE
              QY(IX,IY,IZ)
     .       =2.0/(pm(I_roms,J_roms)+pm(I_roms,J_roms-1))
     .       *0.5*( D(I_roms,J_roms)+ D(I_roms,J_roms-1))
     .       *DZ(IZ)*v(I_roms,J_roms-1,NZ+1-IZ)
     .       *maskv(I_roms,J_roms-1)
           ENDIF
           ! horizontal diffusivity (m2/sec)             
           EX(IX,IY,IZ)=0.95*KH*maskr(I_roms,J_roms)
           EY(IX,IY,IZ)=0.95*kH*maskr(I_roms,J_roms)
         ENDDO

         DO IZ=1,NZ+1
            ! vertical advective fluxs across cell faces (m3/sec)
            QZ(IX,IY,IZ)=1.0/pm(I_roms,J_roms)/pn(I_roms,J_roms)
     .                   *omega(I_roms,J_roms,NZ+2-IZ)
     .                   *maskr(I_roms,J_roms)
            ! vertical diffusivity (m2/sec)
            EZ(IX,IY,IZ)=0.85*AKs(I_roms,J_roms,NZ+2-IZ)
     .                *maskr(I_roms,J_roms)
         ENDDO
      ENDDO
      ENDDO

! Calculate the surface elevation rate of change (m/s)
      DO IY=2,NY-1
      DO IX=2,NX-1
         I_roms=IX-1         
         J_roms=IY-1         
         DO IZ=1,NZ                                       
            DETA(IX,IY)=DETA(IX,IY)
     .                 +QX(IX,IY,IZ)-QX(IX+1,IY,IZ)
     .                 +QY(IX,IY,IZ)-QY(IX,IY+1,IZ)            ! m3/s 
         ENDDO
         DETA(IX,IY)=DETA(IX,IY)
     .              *pm(I_roms,J_roms)*pn(I_roms,J_roms)
     .              *maskr(I_roms,J_roms)                      ! m/s
      ENDDO
      ENDDO

      END SUBROUTINE ncset_transport
