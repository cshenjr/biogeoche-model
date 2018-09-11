!======================================================================
! put all the attributes listed below to variable
!                               Yun Li, UMCES/HPL, Feb-17-2011
!======================================================================

        ! start define dimension and attributes
        status=nf90_redef(ncID)
        CALL nccheck_status(status,'nf90_redef()',SUBNA_h)

        ! dimension
        status=nf90_def_var(ncID,TRIM(ADJUSTL(Vinfo(1,idvar))),
     .         nf90_real,mdimen,varID)
        CALL nccheck_status(status,'nf90_def_var()',SUBNA_h)

        ! attributes
        status=nf90_put_att(ncID,varID,
     .                     'long_name',TRIM(ADJUSTL(Vinfo(2,idvar))))
        CALL nccheck_status(status,'long_name',SUBNA_h)
        status=nf90_put_att(ncID,varID,
     .                     'units',TRIM(ADJUSTL(Vinfo(3,idvar))))
        CALL nccheck_status(status,'unit',SUBNA_h)
        status=nf90_put_att(ncID,varID,
     .                     'field',TRIM(ADJUSTL(Vinfo(4,idvar))))
        CALL nccheck_status(status,'field',SUBNA_h)
        status=nf90_put_att(ncID,varID,
     .                     'time',TRIM(ADJUSTL(Vinfo(5,idvar))))
        CALL nccheck_status(status,'time',SUBNA_h)
        status=nf90_put_att(ncID,varID,
     .                     'coordinates',TRIM(ADJUSTL(coordinates)))
        CALL nccheck_status(status,'coordinate',SUBNA_h)
        status=nf90_put_att(ncID,varID,
     .                     '_FillValue',REAL(spval))
        CALL nccheck_status(status,'_FillValue',SUBNA_h)

        ! end define
        status=nf90_enddef(ncID)
        CALL nccheck_status(status,'nf90_enddef()',SUBNA_h)
