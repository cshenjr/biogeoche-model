      SUBROUTINE RCA13(NXCALL)
C 
C        RCA13 UPDATES THE MISC TIME FUNCTIONS
C 
      SAVE
      INCLUDE 'RCACM' 
      REAL   NT,NXCALL

      DO 300 I=1,NOFUNC 
       ITIME=ITIMF(I)
       IF(TIMEMTF(I,ITIME).LE.TIME)  THEN
        IF(ITVFPWLOPT.EQ.0)  THEN
          MFUNC(I) = 0.
          BFUNC(I) = VALMTF(I,ITIME) 
        ELSE
          MFUNC(I) = (VALMTF(I,ITIME)-VALMTF(I,ITIME+1))/
     .                (TIMEMTF(I,ITIME)-TIMEMTF(I,ITIME+1))
          BFUNC(I) = VALMTF(I,ITIME+1) 
        ENDIF
        NTSECS=ISCALTVF*TIMEMTF(I,ITIME+1)
        NT=NTSECS/86400.
        NXFUNT(I) = NT
        ITIMF(I) = ITIME+1
       ENDIF
  300 CONTINUE

      NXCALL=9999.
      DO I=1,NOFUNC
       IF(NXFUNT(I).LE.NXCALL)   NXCALL=NXFUNT(I)
      ENDDO

      RETURN
      END 
