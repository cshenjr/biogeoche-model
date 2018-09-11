function [pHx] = pH_cal(TempC,Sal,DIC,TA,TP,TSi)
% pH and pco2 calculatio

dt = 1e6;

DIC = DIC/dt;
TA = TA/dt;
TP = TP/dt;
TSi  = TSi/dt;

% fc = fugfac*pc
sqrSal = sqrt(Sal);
RGasConstant = 83.1451;
Pdbar = 0;

TempK    = TempC + 273.15;
RT       = RGasConstant*TempK;
logTempK = log(TempK);
Pbar     = Pdbar/10;

TB =  0.0004157*Sal/35; % in mol/kg-SW
% TB(FF) =  0.0004326.*Sal(FF)./35; % in mol/kg-SW
TF = (0.000067/18.998)*(Sal/1.80655); % in mol/kg-SW
TS = (0.14/96.062)*(Sal/1.80655); % in mol/kg-SW

TempK100  = TempK/100;
lnK0 = -60.2409 + 93.4517./TempK100 + 23.3585*log(TempK100) + Sal.*...
    (0.023517-0.023656*TempK100 + 0.0047036*TempK100.^2);
K0   = exp(lnK0);                  % this is in mol/kg-SW/atm

IonS         = 19.924*Sal./(1000 - 1.005.* Sal);
lnKS = -4276.1./TempK + 141.328 - 23.093*logTempK +...             
      (-13856./TempK + 324.57 - 47.986*logTempK).*sqrt(IonS) +...     
      (35474./TempK - 771.54 + 114.723*logTempK).*IonS +...           
      (-2698./TempK).*sqrt(IonS).*IonS + (1776./TempK).*IonS.^2; 
KS = exp(lnKS)...            % this is on the free pH scale in mol/kg-H2O
        .* (1 - 0.001005 * Sal);   % convert to mol/kg-SW
%     pKS = 647.59 ./ TempK - 6.3451 + 0.019085.*TempK - 0.5208.*sqrt(IonS);
%     KS = 10.^(-pKS)...          % this is on the free pH scale in mol/kg-H2O
%         .* (1 - 0.001005.*Sal);    % convert to mol/kg-SW
lnKF = 1590.2./TempK - 12.641 + 1.525*IonS.^0.5;
KF   = exp(lnKF)...                 % this is on the free pH scale in mol/kg-H2O
    .*(1 - 0.001005.*Sal);          % convert to mol/kg-SW

SWStoTOT  = (1 + TS./KS)./(1 + TS./KS + TF./KF);

fH = 1.2948 - 0.002036*TempK + (0.0004607 -...
        0.000001475*TempK).*Sal.^2;

lnKBtop = -8966.9 - 2890.53*sqrSal - 77.942*Sal +...
        1.728*sqrSal.*Sal - 0.0996.*Sal.^2;
lnKB = lnKBtop./TempK + 148.0248 + 137.1942.*sqrSal +...
        1.62142.*Sal + (-24.4344 - 25.085.*sqrSal - 0.2474.*...
        Sal).*logTempK + 0.053105.*sqrSal.*TempK;
KB = exp(lnKB)...    % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT; 

lnKW = 148.9802 - 13847.26./TempK - 23.6521*logTempK +...
        (-5.977 + 118.67./TempK + 1.0495*logTempK).*...
        sqrSal - 0.01615*Sal;
KW = exp(lnKW);

lnKP1 = -4576.752./TempK + 115.54 - 18.453*logTempK + (-106.736./TempK +...
        0.69171).*sqrSal + (-0.65643./TempK - 0.01844).*Sal;
KP1 = exp(lnKP1);
lnKP2 = -8814.715./TempK + 172.1033 - 27.927*logTempK + (-160.34./TempK +...
        1.3566).*sqrSal + (0.37335./TempK - 0.05778).*Sal;
KP2 = exp(lnKP2);
lnKP3 = -3070.75./TempK - 18.126 + (17.27039./TempK + 2.81197).*sqrSal +...
        (-44.99486./TempK - 0.09984).*Sal;
KP3 = exp(lnKP3);
lnKSi = -8904.2./TempK + 117.4 - 19.334*logTempK + (-458.79./TempK +...
        3.5913).*sqrt(IonS) + (188.74./TempK - 1.5998).*IonS +...
        (-12.1652./TempK + 0.07871).*IonS.^2;
KSi = exp(lnKSi)...                % this is on the SWS pH scale in mol/kg-H2O
        .*(1 - 0.001005.*Sal); 

% cai and wang's method for k1 k2
F1 = 200.1./TempK + 0.3220;
pK1 = 3404.71./TempK + 0.032786*TempK - 14.8435 - 0.071692*F1.*Sal.^0.5 + 0.0021487*Sal;
K1  = 10.^-pK1...         % this is on the NBS scale
        ./fH;                    % convert to SWS scale (uncertain at low Sal due to junction potential);
F2 = -129.24./TempK + 1.4381;
pK2 = 2902.39./TempK + 0.02379*TempK - 6.4980 - 0.3191*F2.*Sal.^0.5 + 0.0198*Sal;
K2  = 10.^-pK2...         % this is on the NBS scale
        ./fH;                    % convert to SWS scale (uncertain at low Sal due to junction potential); 
%------------------------------------
deltaV  = -25.5 + 0.1271*TempC;
    %                 'deltaV = deltaV - .151.*(Sali - 34.8); % Millero, 1979
Kappa   = (-3.08 + 0.0877*TempC)/1000;
    %                 'Kappa = Kappa  - .578.*(Sali - 34.8)/1000.; % Millero, 1979
lnK1fac = (-deltaV + 0.5*Kappa.*Pbar).*Pbar./RT;

deltaV  = -15.82 - 0.0219*TempC;
    %                  'deltaV = deltaV + .321.*(Sali - 34.8); % Millero, 1979
Kappa   = (1.13 - 0.1475*TempC)/1000;
    %                 'Kappa = Kappa - .314.*(Sali - 34.8)./1000: % Millero, 1979
lnK2fac = (-deltaV + 0.5*Kappa.*Pbar).*Pbar./RT;

deltaV  = -29.48 + 0.1622*TempC - 0.002608*TempC.^2;

Kappa   = -2.84/1000; % Millero, 1979

lnKBfac = (-deltaV + 0.5*Kappa.*Pbar).*Pbar./RT;

deltaV  = -20.02 + 0.1119*TempC - 0.001409*TempC.^2;
    %               Millero, 1992 and Millero, 1995 have:
Kappa   = (-5.13 + 0.0794*TempC)/1000; % Millero, 1983
    %               Millero, 1995 has this too, but Millero, 1992 is different.
lnKWfac = (-deltaV + 0.5*Kappa.*Pbar).*Pbar./RT;

deltaV = -9.78 - 0.009*TempC - 0.000942*TempC.^2;
Kappa = (-3.91 + 0.054*TempC)/1000;
lnKFfac = (-deltaV + 0.5*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKS:
%       This is from Millero, 1995, which is the same as Millero, 1983.
%       It is assumed that KS is on the free pH scale.
deltaV = -18.03 + 0.0466*TempC + 0.000316*TempC.^2;
Kappa = (-4.53 + 0.09*TempC)/1000;
lnKSfac = (-deltaV + 0.5*Kappa.*Pbar).*Pbar./RT;

deltaV = -14.51 + 0.1211*TempC - 0.000321*TempC.^2;
Kappa  = (-2.67 + 0.0427*TempC)/1000;
lnKP1fac = (-deltaV + 0.5*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKP2:
deltaV = -23.12 + 0.1758*TempC - 0.002647*TempC.^2;
Kappa  = (-5.15 + 0.09 *TempC)/1000;
lnKP2fac = (-deltaV + 0.5*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKP3:
deltaV = -26.57 + 0.202*TempC - 0.003042*TempC.^2;
Kappa  = (-4.08 + 0.0714*TempC)/1000;
lnKP3fac = (-deltaV + 0.5*Kappa.*Pbar).*Pbar./RT;

deltaV = -29.48 + 0.1622*TempC - 0.002608*TempC.^2;
Kappa  = -2.84/1000;
lnKSifac = (-deltaV + 0.5*Kappa.*Pbar).*Pbar./RT;

K1fac  = exp(lnK1fac);  K1  = K1 .*K1fac;
K2fac  = exp(lnK2fac);  K2  = K2 .*K2fac;
KWfac  = exp(lnKWfac);  KW  = KW .*KWfac;
KBfac  = exp(lnKBfac);  KB  = KB .*KBfac;
KFfac  = exp(lnKFfac);  KF  = KF .*KFfac;
KSfac  = exp(lnKSfac);  KS  = KS .*KSfac;
KP1fac = exp(lnKP1fac); KP1 = KP1.*KP1fac;
KP2fac = exp(lnKP2fac); KP2 = KP2.*KP2fac;
KP3fac = exp(lnKP3fac); KP3 = KP3.*KP3fac;
KSifac = exp(lnKSifac); KSi = KSi.*KSifac;

SWStoTOT  = (1 + TS./KS)./(1 + TS./KS + TF./KF);
FREEtoTOT =  1 + TS./KS;

% fugfac vapor pressure
Delta = (57.7 - 0.118*TempK);
b = -1636.75 + 12.0408*TempK - 0.0327957*TempK.^2 + 3.16528*0.00001*TempK.^3;
% For a mixture of CO2 and air at 1 atm (at low CO2 concentrations);
P1atm = 1.01325; % in bar
FugFac = exp((b + 2.*Delta).*P1atm./RT);

VPWP = exp(24.4543 - 67.4509*(100./TempK) - 4.8489*log(TempK/100));
VPCorrWP = exp(-0.000544*Sal);
VPSWWP = VPWP.*VPCorrWP;
VPFac = 1 - VPSWWP; % this assumes 1 atmosphere

%  calculate ph
%  ------------------------------------------------
TCx = DIC; TAx = TA;
K1F=K1;   K2F=K2;   KWF =KW;
KP1F=KP1; KP2F=KP2; KP3F=KP3;  TPF=TP;
TSiF=TSi; KSiF=KSi; TBF =TB;   KBF=KB;
TSF =TS;  KSF =KS;  TFF =TF;   KFF=KF;

pHGuess     = 8;       % this is the first guess
pHTol       = 0.0001;  % tolerance for iterations end
ln10        = log(10); %
pHx = pHGuess; % creates a vector holding the first guess for all samples
deltapH     = pHTol+1;

kk = 0;

for i = 1:10
% while any(abs(deltapH) > pHTol)
    kk = kk+1;
    H         = 10.^(-pHx);
    Denom     = (H.*H + K1F.*H + K1F.*K2F);
    CAlk      = TCx.*K1F.*(H + 2.*K2F)./Denom;
    BAlk      = TBF.*KBF./(KBF + H);
    OH        = KWF./H;
    PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
    PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
    PAlk      = TPF.*PhosTop./PhosBot;
    SiAlk     = TSiF.*KSiF./(KSiF + H);
    FREEtoTOT = (1 + TSF./KSF); % pH scale conversion factor
    Hfree     = H./FREEtoTOT; % for H on the total scale
    HSO4      = TSF./(1 + KSF./Hfree); % since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree); % since KF is on the free scale
    Residual  = TAx - CAlk - BAlk - OH - PAlk - SiAlk + Hfree + HSO4 + HF;
    % find Slope dTA/dpH;
    % (this is not exact, but keeps all important terms);
    Slope     = ln10.*(TCx.*K1F.*H.*(H.*H + K1F.*K2F + 4.*H.*K2F)./Denom./Denom + BAlk.*H./(KBF + H) + OH + H);
    deltapH   = Residual./Slope; % this is Newton's method
    % to keep the jump from being too big;
    while any(abs(deltapH) > 1)
        deltapH=deltapH./2;
    end
    pHx       = pHx + deltapH; % Is on the same scale as K1 and K2 were calculated...
end

end

