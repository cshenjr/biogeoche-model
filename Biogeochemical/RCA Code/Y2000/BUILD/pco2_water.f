      SUBROUTINE PCO2_WATER(TempC,Sal,TIC,TALK,PO4,SIO3,PH,PCO2)
C
C      Compute the equilibrium partial pressure of co2 in surface water
C
      IMPLICIT NONE
	  
	     REAL, INTENT(IN) :: TempC,Sal,TIC,TALK,PO4,SIO3,PH
		 REAL, INTENT(OUT) :: PCO2
		 
	     REAL  DIC,TA,TP,TSi,sqrSal,RGasConstant,TempK
		 REAL  RT,logTempK,Pbar,TB,TF,TS,TempK100,lnK0
		 REAL  K0,IonS,lnKS,KS,lnKF,KF,SWStoTOT,FREEtoTOT
		 REAL  fH,lnKBtop,lnKB,KB,lnKW,KW,lnKP1,KP1
		 REAL  lnKP2,KP2,lnKP3,KP3,lnKSi,KSi,F1,pK1,K1
		 REAL  F2,pK2,K2,Delta,b,P1atm,FugFac,pHTol,ln10
		 REAL  pHx,deltapH,H,Denom,CAlk,BAlk,OH,PhosTop
		 REAL  PhosBot,PAlk,SiAlk,Hfree,HSO4,HF,Residual
		 REAL  Slope,H1
		 
		 INTEGER  NEWTON,INEWTON
		 
		 INEWTON = 10
		 
		 DIC = TIC*0.000001
	     TA = TALK*0.000001
	     TP = PO4*0.000001
	     TSi = SIO3*0.000001
		 
		 sqrSal = sqrt(Sal)
	     RGasConstant = 83.1451
	  
	     TempK = TempC + 273.15
		 RT = RGasConstant*TempK
		 logTempK = log(TempK)
		 Pbar = 0.0/10
		 
		 TB =  0.0004157*Sal/35
		 TF = (0.000067/18.998)*(Sal/1.80655)
		 TS = (0.14/96.062)*(Sal/1.80655)
		 
		 TempK100  = TempK/100
		 lnK0 = -60.2409 + 93.4517/TempK100 + 23.3585*log(TempK100) + Sal
     .          *(0.023517 - 0.023656*TempK100 + 0.0047036*TempK100**2)
         K0   = exp(lnK0); 
		 
		 IonS = 19.924*Sal/(1000 - 1.005* Sal)
         lnKS = -4276.1/TempK + 141.328 - 23.093*logTempK +
     .         (-13856/TempK + 324.57 - 47.986*logTempK)*sqrt(IonS) +
     .         (35474/TempK - 771.54 + 114.723*logTempK)*IonS +
     .         (-2698/TempK)*sqrt(IonS)*IonS + (1776/TempK)*IonS**2
		 KS = exp(lnKS)*(1 - 0.001005*Sal)
		 
		 lnKF = 1590.2/TempK - 12.641 + 1.525*IonS**0.5
		 KF   = exp(lnKF)*(1 - 0.001005*Sal)
		 
		 SWStoTOT  = (1 + TS/KS)/(1 + TS/KS + TF/KF)
		 FREEtoTOT =  1 + TS/KS
		 
		 fH = 1.2948 - 0.002036*TempK + (0.0004607 -
     .		 0.000001475*TempK)*Sal**2
		 
		 lnKBtop = -8966.9 - 2890.53*sqrSal - 77.942*Sal +
     .		 1.728*sqrSal*Sal - 0.0996*Sal**2
		 lnKB = lnKBtop/TempK + 148.0248 + 137.1942*sqrSal +
     .		 1.62142*Sal + (-24.4344 - 25.085*sqrSal - 0.2474*
     .		 Sal)*logTempK + 0.053105*sqrSal*TempK
		 KB = exp(lnKB)/SWStoTOT
		 
		 lnKW = 148.9802 - 13847.26/TempK - 23.6521*logTempK +
     .		 (-5.977 + 118.67/TempK + 1.0495*logTempK)*
     .		 sqrSal - 0.01615*Sal
		 KW = exp(lnKW)
		 
		 lnKP1= -4576.752/TempK+115.54-18.453*logTempK+(-106.736/TempK+
     .		0.69171)*sqrSal + (-0.65643/TempK - 0.01844)*Sal 
		 KP1 = exp(lnKP1)
		 
		 lnKP2 =-8814.715/TempK+172.1033-27.927*logTempK+(-160.34/TempK+
     .		 1.3566)*sqrSal+(0.37335/TempK-0.05778)*Sal
		 KP2= exp(lnKP2)
		 
		 lnKP3=-3070.75/TempK-18.126+(17.27039/TempK+2.81197)*sqrSal+
     .		 (-44.99486/TempK-0.09984)*Sal
		 KP3 = exp(lnKP3)
		 
		 lnKSi= -8904.2/TempK + 117.4 - 19.334*logTempK + (-458.79/TempK +
     .		 3.5913)*sqrt(IonS) + (188.74/TempK - 1.5998)*IonS +
     .		 (-12.1652/TempK + 0.07871)*IonS**2
		 KSi = exp(lnKSi)*(1 - 0.001005*Sal)
		 
		 F1 = 200.1/TempK + 0.3220
		 pK1= 3404.71/TempK+0.032786*TempK-14.8435-0.071692*F1*
     .		 Sal**0.5 + 0.0021487*Sal
		 K1 = (10**(-pK1))/fH
		 
		 F2 = -129.24/TempK + 1.4381
		 pK2= 2902.39/TempK + 0.02379*TempK - 6.4980 - 0.3191
     .		 *F2*Sal**0.5 + 0.0198*Sal
		 K2= (10**(-pK2))/fH
		 
		 SWStoTOT  = (1 + TS/KS)/(1 + TS/KS + TF/KF)
		 FREEtoTOT =  1 + TS/KS
		 
		 Delta = (57.7 - 0.118*TempK)
		 b=-1636.75+12.0408*TempK-0.0327957*TempK**2+3.16528*0.00001*TempK**3
		 P1atm = 1.01325
		 FugFac = exp((b + 2*Delta)*P1atm/RT)
		 
		 pHTol       = 0.0001
		 ln10        = 2.3026
		 pHx = PH
		 deltapH     = pHTol+1
		 
!		 DO WHILE (abs(deltapH) .GT. pHTol)
         DO NEWTON = 1,INEWTON
		     H = 10**(-pHx)
		     Denom =(H*H + K1*H + K1*K2)
		     CAlk = DIC*K1*(H + 2*K2)/Denom
		     BAlk = TB*KB/(KB + H)
		     OH   = KW/H
		     PhosTop=KP1*KP2*H + 2*KP1*KP2*KP3 - H*H*H
		     PhosBot= H*H*H+KP1*H*H+KP1*KP2*H+KP1*KP2*KP3
		     PAlk   = TP*PhosTop/PhosBot
		     SiAlk  = TSi*KSi/(KSi + H)
		     FREEtoTOT = (1 + TS/KS)
		     Hfree     = H/FREEtoTOT
		     HSO4      = TS/(1 + KS/Hfree)
		     HF        = TF/(1 + KF/Hfree)
		     Residual  = TA - CAlk - BAlk - OH - PAlk - SiAlk 
     .			       + Hfree + HSO4 + HF
		     Slope     = ln10*(DIC*K1*H*(H*H + K1*K2 + 4*H*K2)/Denom/Denom 
     .			       + BAlk*H/(KB + H) + OH + H)
	         deltapH   = Residual/Slope
			 IF (abs(deltapH) .GT. 1) THEN
			     deltapH=deltapH/2
			 ENDIF
			 pHx       = pHx + deltapH
		 END DO
		 
		 H1 = 10**(-pHx)
		 PCO2 = DIC*H1*H1/(H1*H1 + K1*H1 + K1*K2)/K0
		 PCO2 = PCO2*1e6
			 
		 RETURN
			 
      END 
