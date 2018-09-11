!======================================================================
! inquire dimension ID for RCA variables
!                              Yun Li, UMCES/HPL, Feb-18-2011
!======================================================================

      status=nf90_inq_dimid(ncID,'NX',NXdimID)
      CALL nccheck_status(status,'NX',TRIM(SUBNA_h))
      status=nf90_inq_dimid(ncID,'NY',NYdimID)
      CALL nccheck_status(status,'NY',TRIM(SUBNA_h))
      status=nf90_inq_dimid(ncID,'NZ',NZdimID)
      CALL nccheck_status(status,'NZ',TRIM(SUBNA_h))
      status=nf90_inq_dimid(ncID,'NZZ',NZZdimID)
      CALL nccheck_status(status,'NZZ',TRIM(SUBNA_h))
      status=nf90_inq_dimid(ncID,'group',NGdimID)
      CALL nccheck_status(status,'group',TRIM(SUBNA_h))
      status=nf90_inq_dimid(ncID,'NSYS',NSYSdimID)
      CALL nccheck_status(status,'NSYS',TRIM(SUBNA_h))
      status=nf90_inq_dimid(ncID,'TIME',TdimID)
      CALL nccheck_status(status,'TIME',TRIM(SUBNA_h))
      status=nf90_inq_dimid(ncID,'CONST',CONSTdimID)
      CALL nccheck_status(status,'CONST',TRIM(SUBNA_h))
