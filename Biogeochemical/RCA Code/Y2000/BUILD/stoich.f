      SUBROUTINE STOICH(NALGRP,NUTR,KADNUTR,PTAM,PC1,PC2,PC3
     .   ,CNUTR1,CNUTR2,CNUTR3,ANUTRC1,ANUTRC2,ANUTRC3)

      INCLUDE   'RCACM' 

C        -STOICH- DETERMINES THE NUTRIENT TO CARBON RATIOS GIVEN ...
C  INPUTS
C     NALGRP  =  NUMBER OF ALGAL GROUPS IN MODEL APPLICATION (=1, 2, OR 3)
C     NUTR    =  TOTAL INORGANIC NUTRIENT CONCENTRATION
C     KADNUTR =  PARTITION COEFFICIENT
C     PTAM    =  PARTICULATE TOTAL ACTIVE METAL
C     PC1     =  BIOMASS CONCENTRATION OF PHYTOPLANKTON GROUP 1
C     PC2     =  BIOMASS CONCENTRATION OF PHYTOPLANKTON GROUP 2
C     PC3     =  BIOMASS CONCENTRATION OF PHYTOPLANKTON GROUP 3
C     CNUTR1  =  CARBON TO NUTRIENT COEFFICIENTS FOR PHYTOPLANKTON GROUP 1
C     CNUTR2  =  CARBON TO NUTRIENT COEFFICIENTS FOR PHYTOPLANKTON GROUP 2
C     CNUTR3  =  CARBON TO NUTRIENT COEFFICIENTS FOR PHYTOPLANKTON GROUP 3
C  OUTPUTS
C     ANUTRC1 =  NUTRIENT TO CARBON RATIO FOR PHYTOPLANKTON GROUP 1
C     ANUTRC2 =  NUTRIENT TO CARBON RATIO FOR PHYTOPLANKTON GROUP 2
C     ANUTRC3 =  NUTRIENT TO CARBON RATIO FOR PHYTOPLANKTON GROUP 3

      REAL   NUTR,KADNUTR,CNUTR1(3),CNUTR2(3),CNUTR3(3)

      IF(NALGRP.LE.0 .OR. NALGRP.GE.4) GO TO 500

      CN11=CNUTR1(1)
      CN12=CNUTR1(2)
      CN13=CNUTR1(3)
      CN21=CNUTR2(1)
      CN22=CNUTR2(2)
      CN23=CNUTR2(3)
      CN31=CNUTR3(1)
      CN32=CNUTR3(2)
      CN33=CNUTR3(3)

      XL     = 0.0
      XR     = NUTR
      XMOLD  = NUTR+1.E-30
      ALGCAR = PC1+PC2+PC3

      GO TO (100,200,300), NALGRP

  100 IF (ALGCAR.GT.0.0.AND.NUTR.GT.0.0) THEN
        DF = 1./(1.+KADNUTR*PTAM)
        DO 120 I=1,10
         FXL = XL + PC1/(CN11+CN12*EXP(-MIN(10.,CN13*XL*DF)))
     .            - NUTR
         FXR = XR + PC1/(CN11+CN12*EXP(-MIN(10.,CN13*XR*DF)))
     .            - NUTR
         XM  = (XL*FXR-XR*FXL)/(FXR-FXL)
         FXM = XM + PC1/(CN11+CN12*EXP(-MIN(10.,CN13*XM*DF)))
     .            - NUTR
         PCTDF = ABS((XM-XMOLD)/XMOLD)
         IF (PCTDF.LT.0.01.OR.FXM.EQ.0.) GO TO 130
         IF (FXM*FXL.LT.0.) THEN
           XR = XM
         ELSE
           XL = XM
         ENDIF
         XMOLD = XM
  120   CONTINUE

  130   ANUTRC1 = (NUTR-XM)/PC1
        ANUTRC1 = MAX(ANUTRC1,1./(CN11+CN12))
        ANUTRC1 = MIN(ANUTRC1,1./CN11)
 
      ELSE
 
        ANUTRC1 = 0.0

      ENDIF
!=================================================================
! add system bypass  Yun Li, Aug-18-2011
      IF(SYSBY(2).EQ.1) ANUTRC1 = 0.0
!=================================================================
      RETURN

  200 IF (ALGCAR.GT.0.0.AND.NUTR.GT.0.0) THEN
        DF = 1./(1.+KADNUTR*PTAM)
        DO 220 I=1,10
         FXL = XL + PC1/(CN11+CN12*EXP(-MIN(10.,CN13*XL*DF)))
     .            + PC2/(CN21+CN22*EXP(-MIN(10.,CN23*XL*DF)))
     .            - NUTR
         FXR = XR + PC1/(CN11+CN12*EXP(-MIN(10.,CN13*XR*DF)))
     .            + PC2/(CN21+CN22*EXP(-MIN(10.,CN23*XR*DF)))
     .            - NUTR
         XM  = (XL*FXR-XR*FXL)/(FXR-FXL)
         FXM = XM + PC1/(CN11+CN12*EXP(-MIN(10.,CN13*XM*DF)))
     .            + PC2/(CN21+CN22*EXP(-MIN(10.,CN23*XM*DF)))
     .            - NUTR
         PCTDF = ABS((XM-XMOLD)/XMOLD)
         IF (PCTDF.LT.0.01.OR.FXM.EQ.0.) GO TO 230
         IF (FXM*FXL.LT.0.) THEN
           XR = XM
         ELSE
           XL = XM
         ENDIF
         XMOLD = XM
  220   CONTINUE

  230   ANUTRC1 = (NUTR-XM
     .               -PC2/(CN21+CN22*EXP(-MIN(10.,CN23*XM*DF))))/PC1
        ANUTRC1 = MAX(ANUTRC1,1./(CN11+CN12))
        ANUTRC1 = MIN(ANUTRC1,1./CN11)
        ANUTRC2 = (NUTR-XM-PC1*ANUTRC1)/PC2
        ANUTRC2 = MAX(ANUTRC2,1./(CN21+CN22))
        ANUTRC2 = MIN(ANUTRC2,1./CN21)
 
      ELSE
 
        ANUTRC1 = 0.0
        ANUTRC2 = 0.0

      ENDIF
!=================================================================
! add system bypass  Yun Li, Aug-18-2011
      IF(SYSBY(2).EQ.1) ANUTRC1 = 0.0
      IF(SYSBY(3).EQ.1) ANUTRC2 = 0.0
!=================================================================
      RETURN

  300 IF (ALGCAR.GT.0.0.AND.NUTR.GT.0.0) THEN
        DF = 1./(1.+KADNUTR*PTAM)
        DO 320 I=1,10
         FXL = XL + PC1/(CN11+CN12*EXP(-MIN(10.,CN13*XL*DF)))
     .            + PC2/(CN21+CN22*EXP(-MIN(10.,CN23*XL*DF)))
     .            + PC3/(CN31+CN32*EXP(-MIN(10.,CN33*XL*DF)))
     .            - NUTR
         FXR = XR + PC1/(CN11+CN12*EXP(-MIN(10.,CN13*XR*DF)))
     .            + PC2/(CN21+CN22*EXP(-MIN(10.,CN23*XR*DF)))
     .            + PC3/(CN31+CN32*EXP(-MIN(10.,CN33*XR*DF)))
     .            - NUTR
         XM  = (XL*FXR-XR*FXL)/(FXR-FXL)
         FXM = XM + PC1/(CN11+CN12*EXP(-MIN(10.,CN13*XM*DF)))
     .            + PC2/(CN21+CN22*EXP(-MIN(10.,CN23*XM*DF)))
     .            + PC3/(CN31+CN32*EXP(-MIN(10.,CN33*XM*DF)))
     .            - NUTR
         PCTDF = ABS((XM-XMOLD)/XMOLD)
         IF (PCTDF.LT.0.01.OR.FXM.EQ.0.) GO TO 330
         IF (FXM*FXL.LT.0.) THEN
           XR = XM
         ELSE
           XL = XM
         ENDIF
          XMOLD = XM
  320   CONTINUE

  330   ANUTRC1 = (NUTR-XM
     .               -PC2/(CN21+CN22*EXP(-MIN(10.,CN23*XM*DF)))
     .               -PC3/(CN31+CN32*EXP(-MIN(10.,CN33*XM*DF))))/PC1
        ANUTRC1 = MAX(ANUTRC1,1./(CN11+CN12))
        ANUTRC1 = MIN(ANUTRC1,1./CN11)
        ANUTRC2 = (NUTR-XM-PC1*ANUTRC1
     .               -PC3/(CN31+CN32*EXP(-MIN(10.,CN33*XM*DF))))/PC2
        ANUTRC2 = MAX(ANUTRC2,1./(CN21+CN22))
        ANUTRC2 = MIN(ANUTRC2,1./CN21)
        ANUTRC3 = (NUTR-XM-PC1*ANUTRC1-PC2*ANUTRC2)/PC3
 
      ELSE
 
        ANUTRC1 = 0.0
        ANUTRC2 = 0.0
        ANUTRC3 = 0.0

      ENDIF
!=================================================================
! add system bypass  Yun Li, Aug-18-2011
      IF(SYSBY(2).EQ.1) ANUTRC1 = 0.0
      IF(SYSBY(3).EQ.1) ANUTRC2 = 0.0
      IF(SYSBY(4).EQ.1) ANUTRC3 = 0.0
	  ANUTRC3 = 0.0
!=================================================================
      RETURN

  500 WRITE(6,1000)  NALGRP
 1000 FORMAT(//5X,'ERROR -SUBR STOICH- CALLED WITH INVALID NUMBER OF ALG
     .AL GROUPS'/5X,'NALGRP =',I5/5X,'RCA TERMINATED')
      CALL EXIT
      END

