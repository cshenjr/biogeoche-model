      FUNCTION IUNITCHECK(TWARP,NAME)
      SAVE
      CHARACTER  TWARP*4,NAME*8
      INTEGER  OUT
      COMMON  IN,OUT
      IF(TWARP.NE.'SECS' .AND. TWARP.NE.'secs' .AND.
     .   TWARP.NE.'MINS' .AND. TWARP.NE.'mins' .AND.
     .   TWARP.NE.'HRS ' .AND. TWARP.NE.'hrs ' .AND.
     .   TWARP.NE.' HRS' .AND. TWARP.NE.' hrs' .AND.
     .   TWARP.NE.'DAYS' .AND. TWARP.NE.'days')  THEN
       WRITE(OUT,9000)  NAME,TWARP
 9000  FORMAT(/5X,'THE UNITS CHOSEN FOR ',A8,1X,A4,' ARE NOT VALID'/
     .   5X,'RCA TERMINATED')
       CALL EXIT
      ENDIF
      IF(TWARP.EQ.'SECS'.OR.TWARP.EQ.'secs')  IUNITCHECK=1
      IF(TWARP.EQ.'MINS'.OR.TWARP.EQ.'mins')  IUNITCHECK=60
      IF(TWARP.EQ.'HRS '.OR.TWARP.EQ.'hrs ')  IUNITCHECK=3600
      IF(TWARP.EQ.' HRS'.OR.TWARP.EQ.' hrs')  IUNITCHECK=3600
      IF(TWARP.EQ.'DAYS'.OR.TWARP.EQ.'days')  IUNITCHECK=86400
      RETURN
      END
