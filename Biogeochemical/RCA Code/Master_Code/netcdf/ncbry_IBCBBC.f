      SUBROUTINE ncbry_IBCBBC(ncID)
!
!=======================================================================
!  Read RCA boundary indices (ICB) and concentrations (BBC) in NetCDF  !
!  format                                                              ! 
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
      INTEGER               ::   IX, IY, IZ, IR, idvar 
      INTEGER               ::   num_obc 
      REAL                  ::   Vvalue(NZ)
      CHARACTER(100)        ::   Vname
!
!-----------------------------------------------------------------------
! name of this subroutine
!-----------------------------------------------------------------------
!
      SUBNA=TRIM(ADJUSTL(RCANA))//' when CALL ncbry_IBCBBC.f'
!
!----------------------------------------------------------------------
! check input type
!----------------------------------------------------------------------
!
      IF(Iinfo(ideast) .EQ.t2dvar .AND.
     .   Iinfo(idwest) .EQ.t2dvar .AND.
     .   Iinfo(idnorth).EQ.t2dvar .AND.
     .   Iinfo(idsouth).EQ.t2dvar) THEN
      ELSE
        WRITE(OUT,*) TRIM(ADJUSTL(SUBNA))
        WRITE(OUT,*) 'ERROR: misuse of ncbry_IBCBBC.f ',
     .               TRIM(ADJUSTL(Vinfo(1,ideast))), ' or ',
     .               TRIM(ADJUSTL(Vinfo(1,idwest))), ' or ',
     .               TRIM(ADJUSTL(Vinfo(1,idsouth))),' or ',
     .               TRIM(ADJUSTL(Vinfo(1,idnorth))),
     .               ' is not 2D time-dependent variable'
        CALL EXIT
      ENDIF
!
!-----------------------------------------------------------------------
! get boundary indices and tracer concentrations
!-----------------------------------------------------------------------
!
      num_obc=0
      idvar=idTvar(ISYS)
!
! riverine boundary if upsteam conc. is prescribed
!
      Vname=TRIM(ADJUSTL(Vinfo(1,idriver)))
     .    //TRIM(ADJUSTL(Vinfo(1,idvar)))
      status=nf90_inq_varid(ncID,Vname,varID)
      CALL nccheck_status(status,Vname,SUBNA)
      DO IR=1,nriver
        status=nf90_get_var(ncID,varID,Vvalue,
     .                      start=(/IR, 1,bryTflag/),
     .                      count=(/ 1,NZ,1/))
        CALL nccheck_status(status,Vname,SUBNA)
        DO IZ=1,NZ
           num_obc=num_obc+1
           IBC(1,num_obc,ISYS)=riverX(IR)
           IBC(2,num_obc,ISYS)=riverY(IR)
           IBC(3,num_obc,ISYS)=IZ
           BBC(  num_obc,ISYS)=Vvalue(IZ)		   
        ENDDO
      ENDDO
!
! ocean boundary -- east
!
      DO IY=2,NY-1
        IF(FSM(NX-1,IY).EQ.-2.0) THEN
          Vname=TRIM(ADJUSTL(Vinfo(1,idvar)))
     .        //TRIM(ADJUSTL(Vinfo(1,ideast)))
          status=nf90_inq_varid(ncID,Vname,varID)
          CALL nccheck_status(status,Vname,SUBN)
          status=nf90_get_var(ncID,varID,Vvalue,
     .                        start=(/IY, 1,bryTflag/),
     .                        count=(/ 1,NZ,1/))
          CALL nccheck_status(status,Vname,SUBNA)
          DO IZ=1,NZ
             num_obc=num_obc+1
             IBC(1,num_obc,ISYS)=NX-1
             IBC(2,num_obc,ISYS)=IY
             IBC(3,num_obc,ISYS)=IZ
			 BBC(  num_obc,ISYS)=Vvalue(IZ)
          ENDDO
        ENDIF
      ENDDO
!
! ocean boundary -- west
!
      DO IY=2,NY-1
        IF(FSM(2,IY).EQ.-2.0) THEN
          Vname=TRIM(ADJUSTL(Vinfo(1,idvar)))
     .        //TRIM(ADJUSTL(Vinfo(1,idwest)))
          status=nf90_inq_varid(ncID,Vname,varID)
          CALL nccheck_status(status,Vname,SUBNA)
          status=nf90_get_var(ncID,varID,Vvalue,
     .                        start=(/IY, 1,bryTflag/),
     .                        count=(/ 1,NZ,1/))
          CALL nccheck_status(status,Vname,SUBNA)
          DO IZ=1,NZ
             num_obc=num_obc+1
             IBC(1,num_obc,ISYS)=2
             IBC(2,num_obc,ISYS)=IY
             IBC(3,num_obc,ISYS)=IZ
             BBC(  num_obc,ISYS)=Vvalue(IZ)
          ENDDO
        ENDIF
      ENDDO
!
! ocean boundary -- north
!
      DO IX=2,NX-1
        IF(FSM(IX,NY-1).EQ.-2.0) THEN
          Vname=TRIM(ADJUSTL(Vinfo(1,idvar)))
     .        //TRIM(ADJUSTL(Vinfo(1,idnorth)))
          status=nf90_inq_varid(ncID,Vname,varID)
          CALL nccheck_status(status,Vname,SUBNA)
          status=nf90_get_var(ncID,varID,Vvalue,
     .                        start=(/IX, 1,bryTflag/),
     .                        count=(/ 1,NZ,1/))
          CALL nccheck_status(status,Vname,SUBNA)
          DO IZ=1,NZ
             num_obc=num_obc+1
             IBC(1,num_obc,ISYS)=IX
             IBC(2,num_obc,ISYS)=NY-1
             IBC(3,num_obc,ISYS)=IZ
             BBC(  num_obc,ISYS)=Vvalue(IZ)
          ENDDO
        ENDIF
      ENDDO
!
! ocean boundary -- south
!
      DO IX=2,NX-1
        IF(FSM(IX,2).EQ.-2.0) THEN
          Vname=TRIM(ADJUSTL(Vinfo(1,idvar)))
     .        //TRIM(ADJUSTL(Vinfo(1,idsouth)))
          status=nf90_inq_varid(ncID,Vname,varID)
          CALL nccheck_status(status,Vname,SUBNA)
          status=nf90_get_var(ncID,varID,Vvalue,
     .                        start=(/IX, 1,bryTflag/),
     .                        count=(/ 1,NZ,1/))
          CALL nccheck_status(status,Vname,SUBNA)
          DO IZ=1,NZ
             num_obc=num_obc+1
             IBC(1,num_obc,ISYS)=IX
             IBC(2,num_obc,ISYS)=2
             IBC(3,num_obc,ISYS)=IZ
             BBC(  num_obc,ISYS)=Vvalue(IZ)
          ENDDO
        ENDIF
      ENDDO

      END SUBROUTINE ncbry_IBCBBC
