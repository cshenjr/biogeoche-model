!======================================================================
! define dimension ID for RCA variables
!                            Yun Li, UMCES/HPL, Feb-18-2011
!======================================================================

      ! define dimensions
      status=nf90_def_dim(ncID,'NX',NX,NXdimID)
      CALL nccheck_status(status,'NX',TRIM(SUBNA_h))
      status=nf90_def_dim(ncID,'NY',NY,NYdimID)
      CALL nccheck_status(status,'NY',TRIM(SUBNA_h))
      status=nf90_def_dim(ncID,'NZ',NZ,NZdimID)
      CALL nccheck_status(status,'NZ',TRIM(SUBNA_h))
      status=nf90_def_dim(ncID,'NZZ',NZ+1,NZZdimID)
      CALL nccheck_status(status,'NZZ',TRIM(SUBNA_h))
      status=nf90_def_dim(ncID,'group',group,NGdimID)
      CALL nccheck_status(status,'group',TRIM(SUBNA_h))
      status=nf90_def_dim(ncID,'NSYS',NOSYS,NSYSdimID)
      CALL nccheck_status(status,'NSYS',TRIM(SUBNA_h))
      status=nf90_def_dim(ncID,'TIME',nf90_unlimited,TdimID)
      CALL nccheck_status(status,'TIME',TRIM(SUBNA_h))
      status=nf90_def_dim(ncID,'CONST',1,CONSTdimID)
      CALL nccheck_status(status,'CONST',TRIM(SUBNA_h))

      ! end define
      status=nf90_enddef(ncID)
      CALL nccheck_status(status,'end define attr.',SUBNA_h)

