      SUBROUTINE ncset_ncparams
!
!=======================================================================
!  This SUBROUTINE contains all the variables associated with forcing  !
!  input and state variable output in RCA.  It does NOT include the    !
!  ROMS hydrodynamic variables .  The I/O model is very generic and    !
!  easy to change or expand.                                           !
!                                                                      !
!  as defined in file 'NetCDFCM'                                       !
!                                                                      !
!  Vinfo      Input/output variables names and attributes:             !
!               Vinfo(1,*)  => field variable name.                    !
!               Vinfo(2,*)  => long-name attribute.                    !
!               Vinfo(3,*)  => units attribute.                        !
!               Vinfo(4,*)  => field type attribute.                   !
!               Vinfo(5,*)  => associated time variable name.          !
!               Vinfo(6,*)  => identification name                     !
!               Vinfo(7,*)  => staggered C-grid variable type          !
!                       'n0dvar' => 0D time-independent variable.      !
!                       'n1dvar' => 1D time-independent variable.      !
!                       'n2dvar' => 2D time-independent variable.      !
!                       't1dvar' => 1D time-dependent variable.        !
!                       't2dvar' => 2D time-dependent variable.        !
!                       'g2dvar' => 2D gridded variable.               !
!                       'g3dvar' => 3D gridded variable.               !
!                                                                      !
!  Fscale     Scale to convert input data to model units.              !
!                                                                      !
!  Iinfo      Input/output fields integer variable grid type.          !
!                                                                      !
!                                 YUN LI, UMCES/HPL, Feb-14-2011       !
!=======================================================================
!
        SAVE
        INCLUDE 'RCACM'
        INCLUDE 'NetCDFCM'
!
!  Local variable declarations.
!
        INTEGER        :: i, idvar
        CHARACTER(100) :: Vtemp
!
!  Initialize IO information variables.
!
        DO i=1,NV
          Iinfo(i)=0
          Fscale(i)=1.0
        END DO
!
!-----------------------------------------------------------------------
!  Open input variable information file.
!-----------------------------------------------------------------------
!
        WRITE(OUT,'(1x,a,t30,a)'),
     .      '  Variable INFO FILE NAME ', TRIM(ADJUSTL(VARNA))
        OPEN (31, FILE=TRIM(ADJUSTL(VARNA)), FORM='formatted',  
     .        STATUS='old', ERR=10)
        GOTO 20
  10    WRITE(OUT,60) TRIM(ADJUSTL(VARNA))
        CALL EXIT
  20    CONTINUE
!
!-----------------------------------------------------------------------
!  Read in variable information.  Ignore blank and comment '!' input lines.
!-----------------------------------------------------------------------
!
        idvar = 0
  30    READ (31,*,ERR=40,END=50) Vtemp
!
!-----------------------------------------------------------------------
!  Load variable data into information arrays.
!-----------------------------------------------------------------------
!
        IF ((LEN_TRIM(Vtemp).gt.0).and.(Vtemp(1:1).ne.'!')) THEN
          idvar = idvar+1
          Vinfo(1,idvar)=Vtemp
          READ (31,*,ERR=30) Vinfo(2,idvar)
          READ (31,*,ERR=30) Vinfo(3,idvar)
          READ (31,*,ERR=30) Vinfo(4,idvar)
          READ (31,*,ERR=30) Vinfo(5,idvar)
          READ (31,*,ERR=30) Vinfo(6,idvar)
          READ (31,*,ERR=30) Vinfo(7,idvar)
          READ (31,*,ERR=30) Fscale(idvar)
!
!-----------------------------------------------------------------------
!  Determine variable type.
!-----------------------------------------------------------------------
!
          SELECT CASE (TRIM(ADJUSTL(Vinfo(7,idvar))))
            CASE ('n0dvar')
              Iinfo(idvar)=n0dvar
            CASE ('n1dvar')
              Iinfo(idvar)=n1dvar
            CASE ('n2dvar')
              Iinfo(idvar)=n2dvar
            CASE ('t1dvar')
              Iinfo(idvar)=t1dvar
            CASE ('t2dvar')
              Iinfo(idvar)=t2dvar
            CASE ('g2dvar')
              Iinfo(idvar)=g2dvar
            CASE ('g3dvar')
              Iinfo(idvar)=g3dvar
            CASE DEFAULT
              Iinfo(idvar)=0
          END SELECT
!
!-----------------------------------------------------------------------
!  Assign identification indices.
!-----------------------------------------------------------------------
!
          SELECT CASE (TRIM(ADJUSTL(Vinfo(6,idvar))))
!
!  model dimension 
!
            CASE ('idNX')
              idNX=idvar
            CASE ('idNY')
              idNY=idvar
            CASE ('idNZ')
              idNZ=idvar
            CASE ('idNZZ')
              idNZZ=idvar
            CASE ('idNG')
              idNG=idvar
            CASE ('idNSYS')
              idNSYS=idvar
!
!  grid info
!
            CASE ('idLAT')
              idLAT=idvar
            CASE ('idLON')
              idLON=idvar
            CASE ('idH')
              idH=idvar
            CASE ('idDX')
              idDX=idvar
            CASE ('idDY')
              idDY=idvar
            CASE ('idDZ')
              iddZ=idvar
            CASE ('idDZZ')
              idDZZ=idvar
            CASE ('idFSM')
              idFSM=idvar
!
!  initial 
!
            CASE ('idini_time')
              idini_time=idvar
!
!  boundary
!
            CASE ('idIBCOPT')
              idIBCOPT=idvar
            CASE ('idIBCPWLOPT')
              idIBCPWLOPT=idvar
            CASE ('idbry_time')
              idbry_time=idvar
            CASE ('idriver')
              idriver=idvar
            CASE ('ideast')
              ideast=idvar
            CASE ('idwest')
              idwest=idvar
            CASE ('idnorth')
              idnorth=idvar
            CASE ('idsouth')
              idsouth=idvar
!
!  10sTest model output
!
            CASE ('id10s_TIME')
              id10s_TIME=idvar
            CASE ('id10s_ETA')
              id10s_ETA=idvar
            CASE ('id10s_CONC')
              id10s_CONC=idvar

            ! Global dump averaging is performed
            CASE ('id10s_CONC_GDA')
              id10s_CONC_GDA=idvar
!
!  pathogens model output
!
            CASE ('idpath_TIME')
              idpath_TIME=idvar
            CASE ('idpath_ETA')
              idpath_ETA=idvar
            CASE ('idpath_SAL')
              idpath_SAL=idvar
            CASE ('idTRACR')
              idTRACR=idvar
            CASE ('idTCOLI')
              idTCOLI=idvar
            CASE ('idFCOLI')
              idFCOLI=idvar
            CASE ('idENTERO')
              idENTERO=idvar

            ! Global dump averaging is performed
            CASE ('idpath_SAL_GDA')
              idpath_SAL_GDA=idvar
            CASE ('idTRACR_GDA')
              idTRACR_GDA=idvar
            CASE ('idTCOLI_GDA')
              idTCOLI_GDA=idvar
            CASE ('idFCOLI_GDA')
              idFCOLI_GDA=idvar
            CASE ('idENTERO_GDA')
              idENTERO_GDA=idvar
!
!  eutrophication model input
!
            CASE ('idRPOP')
              idRPOP=idvar
            CASE ('idLPOP')
              idLPOP=idvar
            CASE ('idRDOP')
              idRDOP=idvar
            CASE ('idLDOP')
              idLDOP=idvar
            CASE ('idPO4T')
              idPO4T=idvar
            CASE ('idRPON')
              idRPON=idvar
            CASE ('idLPON')
              idLPON=idvar
            CASE ('idRDON')
              idRDON=idvar
            CASE ('idLDON')
              idLDON=idvar
            CASE ('idNH4T')
              idNH4T=idvar
            CASE ('idBSI')
              idBSI=idvar
            CASE ('idSIT')
              idSIT=idvar
            CASE ('idRPOC')
              idRPOC=idvar
            CASE ('idLPOC')
              idLPOC=idvar
            CASE ('idRDOC')
              idRDOC=idvar
            CASE ('idLDOC')
              idLDOC=idvar
            CASE ('idEXDOC')
              idEXDOC=idvar
            CASE ('idREPOC')
              idREPOC=idvar
            CASE ('idREDOC')
              idREDOC=idvar
            CASE ('idDO')
              idDO=idvar
!
!  eutrophication model parameters
!
            CASE ('idKL')
              idKL=idvar
            CASE ('idVSNET1')
              idVSNET1=idvar
            CASE ('idVSNET2')
              idVSNET2=idvar
            CASE ('idVSNET3')
              idVSNET3=idvar
            CASE ('idVSNET4')
              idVSNET4=idvar
            CASE ('idKEBS')
              idKEBS=idvar
            CASE ('idSSLDS')
              idSSLDS=idvar
            CASE ('idITOTSF')
              idITOTSF=idvar
            CASE ('idF')
              idF=idvar
            CASE ('idWIND')
              idWIND=idvar
            CASE ('idKETVF')
              idKETVF=idvar
            CASE ('idTSS')
              idTSS=idvar
            CASE ('idITOTSF_time')
              idITOTSF_time=idvar
            CASE ('idF_time')
              idF_time=idvar
            CASE ('idWIND_time')
              idWIND_time=idvar
            CASE ('idKETVF_time')
              idKETVF_time=idvar
            CASE ('idTSS_time')
              idTSS_time=idvar
!
!  eutrophication model output
!
            CASE ('ideutro_TIME')
              ideutro_TIME=idvar
            CASE ('ideutro_ETA')
              ideutro_ETA=idvar
            CASE ('ideutro_SAL')
              ideutro_SAL=idvar
            CASE ('idPHYT1')
              idPHYT1=idvar
            CASE ('idPHYT2')
              idPHYT2=idvar
            CASE ('idPHYT3')
              idPHYT3=idvar
            CASE ('idTPOP')
              idTPOP=idvar
            CASE ('idTDOP')
              idTDOP=idvar
            CASE ('idDPO4')
              idDPO4=idvar
            CASE ('idTPON')
              idTPON=idvar
            CASE ('idTDON')
              idTDON=idvar
            CASE ('idDNH4')
              idDNH4=idvar
            CASE ('idNO23')
              idNO23=idvar
            CASE ('idTPSI')
              idTPSI=idvar
            CASE ('idDSI')
              idDSI=idvar
            CASE ('idTPOC')
              idTPOC=idvar
            CASE ('idTDOC')
              idTDOC=idvar
            CASE ('idO2EQ')
              idO2EQ=idvar
            CASE ('idDOAVEG')
              idDOAVEG=idvar
            CASE ('idDOMING')
              idDOMING=idvar
            CASE ('idDOMAXG')
              idDOMAXG=idvar
            CASE ('idHYDSAL')
              idHYDSAL=idvar
            CASE ('idCHLAVEG')
              idCHLAVEG=idvar
            CASE ('idCHLMING')
              idCHLMING=idvar
            CASE ('idCHLMAXG')
              idCHLMAXG=idvar
            CASE ('idTGPP')
              idTGPP=idvar
            CASE ('idTNPP')
              idTNPP=idvar
            CASE ('idTRESP')
              idTRESP=idvar
            CASE ('idHYDTEMP')
              idHYDTEMP=idvar
            CASE ('idSKE')
              idSKE=idvar

            ! Global dump averaging is performed
            CASE ('ideutro_SAL_GDA')
              ideutro_SAL_GDA=idvar
            CASE ('idPHYT1_GDA')
              idPHYT1_GDA=idvar
            CASE ('idPHYT2_GDA')
              idPHYT2_GDA=idvar
            CASE ('idPHYT3_GDA')
              idPHYT3_GDA=idvar
            CASE ('idTPOP_GDA')
              idTPOP_GDA=idvar
            CASE ('idTDOP_GDA')
              idTDOP_GDA=idvar
            CASE ('idDPO4_GDA')
              idDPO4_GDA=idvar
            CASE ('idTPON_GDA')
              idTPON_GDA=idvar
            CASE ('idTDON_GDA')
              idTDON_GDA=idvar
            CASE ('idDNH4_GDA')
              idDNH4_GDA=idvar
            CASE ('idNO23_GDA')
              idNO23_GDA=idvar
            CASE ('idTPSI_GDA')
              idTPSI_GDA=idvar
            CASE ('idTPOC_GDA')
              idTPOC_GDA=idvar
            CASE ('idTDOC_GDA')
              idTDOC_GDA=idvar
            CASE ('idO2EQ_GDA')
              idO2EQ_GDA=idvar
            CASE ('idDOAVEG_GDA')
              idDOAVEG_GDA=idvar
            CASE ('idDOMING_GDA')
              idDOMING_GDA=idvar
            CASE ('idDOMAXG_GDA')
              idDOMAXG_GDA=idvar
            CASE ('idHYDSAL_GDA')
              idHYDSAL_GDA=idvar
            CASE ('idCHLAVEG_GDA')
              idCHLAVEG_GDA=idvar
            CASE ('idCHLMING_GDA')
              idCHLMING_GDA=idvar
            CASE ('idCHLMAXG_GDA')
              idCHLMAXG_GDA=idvar
            CASE ('idTGPP_GDA')
              idTGPP_GDA=idvar
            CASE ('idTNPP_GDA')
              idTNPP_GDA=idvar
            CASE ('idRESP_GDA')
              idRESP_GDA=idvar
            CASE ('idHYDTMP_GDA')
              idHYDTMP_GDA=idvar
            CASE ('idSKE_GDA')
              idSKE_GDA=idvar
!
!  hydrodynamic forcing diagnostics
!
            CASE ('iddia_QX')
              iddia_QX=idvar
            CASE ('iddia_QY')
              iddia_QY=idvar
            CASE ('iddia_QZ')
              iddia_QZ=idvar
            CASE ('iddia_RX')
              iddia_RX=idvar
            CASE ('iddia_RY')
              iddia_RY=idvar
            CASE ('iddia_RZ')
              iddia_RZ=idvar
!
!  eutrophication model diagnostics
!
            CASE ('iddia_ATTENL')
              iddia_ATTENL=idvar
            CASE ('iddia_GPP1')
              iddia_GPP1=idvar
            CASE ('iddia_RESPR1')
              iddia_RESPR1=idvar
            CASE ('iddia_K1C')
              iddia_K1C=idvar
            CASE ('iddia_RNUTR1')
              iddia_RNUTR1=idvar
            CASE ('iddia_RLGHT1')
              iddia_RLGHT1=idvar
            CASE ('iddia_GITMAX1')
              iddia_GITMAX1=idvar
            CASE ('iddia_K1CNSAT')
              iddia_K1CNSAT=idvar
            CASE ('iddia_RESP1')
              iddia_RESP1=idvar
            CASE ('iddia_GRAZ1')
              iddia_GRAZ1=idvar
            CASE ('iddia_ALG1SET')
              iddia_ALG1SET=idvar
            CASE ('iddia_XEMP11')
              iddia_XEMP11=idvar
            CASE ('iddia_XEMP12')
              iddia_XEMP12=idvar
            CASE ('iddia_XEMP13')
              iddia_XEMP13=idvar
            CASE ('iddia_DAVEI1')
              iddia_DAVEI1=idvar

            CASE ('iddia_GPP2')
              iddia_GPP2=idvar
            CASE ('iddia_RESPR2')
              iddia_RESPR2=idvar
            CASE ('iddia_K2C')
              iddia_K2C=idvar
            CASE ('iddia_RNUTR2')
              iddia_RNUTR2=idvar
            CASE ('iddia_RLGHT2')
              iddia_RLGHT2=idvar
            CASE ('iddia_GITMAX2')
              iddia_GITMAX2=idvar
            CASE ('iddia_K2CNSAT')
              iddia_K2CNSAT=idvar
            CASE ('iddia_RESP2')
              iddia_RESP2=idvar
            CASE ('iddia_GRAZ2')
              iddia_GRAZ2=idvar
            CASE ('iddia_ALG2SET')
              iddia_ALG2SET=idvar
            CASE ('iddia_XEMP21')
              iddia_XEMP21=idvar
            CASE ('iddia_XEMP22')
              iddia_XEMP22=idvar
            CASE ('iddia_XEMP23')
              iddia_XEMP23=idvar
            CASE ('iddia_DAVEI2')
              iddia_DAVEI2=idvar
!!!SCQ			  
            CASE ('iddia_ANH4')
              iddia_ANH4=idvar
            CASE ('iddia_ANO23')
              iddia_ANO23=idvar
            CASE ('iddia_CGR')
              iddia_CGR=idvar
            CASE ('iddia_CRE')
              iddia_CRE=idvar
            CASE ('iddia_COX')
              iddia_COX=idvar
            CASE ('iddia_CDI')
              iddia_CDI=idvar
            CASE ('iddia_ACH4')
              iddia_ACH4=idvar
            CASE ('iddia_CA')
              iddia_CA=idvar			  
            CASE ('iddia_POC2')
              iddia_POC2=idvar
            CASE ('iddia_CSOD')
              iddia_CSOD=idvar
            CASE ('iddia_ASNH4')
              iddia_ASNH4=idvar
            CASE ('iddia_ANO3')
              iddia_ANO3=idvar			  
            CASE ('iddia_SUL')
              iddia_SUL=idvar	
            CASE ('iddia_AIRSEA')
              iddia_AIRSEA=idvar			  
!     
!  sediment flux  model input -- restart and initial
!
            CASE ('idCTEMP_ini')
              idCTEMP_ini=idvar
            CASE ('idCPOP')
              idCPOP=idvar
            CASE ('idCPON')
              idCPON=idvar
            CASE ('idCPOC')
              idCPOC=idvar
            CASE ('idCPOS')
              idCPOS=idvar
            CASE ('idPO4T2TM1S')
              idPO4T2TM1S=idvar
            CASE ('idNH4T2TM1S')
              idNH4T2TM1S=idvar
            CASE ('idNO3T2TM1S')
              idNO3T2TM1S=idvar
            CASE ('idHST2TM1S')
              idHST2TM1S=idvar
            CASE ('idSIT2TM1S')
              idSIT2TM1S=idvar
            CASE ('idBNTHSTR1S')
              idBNTHSTR1S=idvar
            CASE ('idPO4T1TM1S')
              idPO4T1TM1S=idvar
            CASE ('idNH4T1TM1S')
              idNH4T1TM1S=idvar
            CASE ('idNO3T1TM1S')
              idNO3T1TM1S=idvar
            CASE ('idHST1TM1S')
              idHST1TM1S=idvar
            CASE ('idSIT1TM1S')
              idSIT1TM1S=idvar
            CASE ('idCH4T1TM1S')
              idCH4T1TM1S=idvar
            CASE ('idCH4T2TM1S')
              idCH4T2TM1S=idvar
            CASE ('idSO4T2TM1S')
              idSO4T2TM1S=idvar
            CASE ('idCH41TM1S')
              idCH41TM1S=idvar
            CASE ('idPO41TM1S')
              idPO41TM1S=idvar
            CASE ('idNH41TM1S')
              idNH41TM1S=idvar
            CASE ('idNO31TM1S')
              idNO31TM1S=idvar
            CASE ('idHS1TM1S')
              idHS1TM1S=idvar
            CASE ('idSI1TM1S')
              idSI1TM1S=idvar
            CASE ('idO20TM1S')
              idO20TM1S=idvar
            CASE ('idSODTM1S')
              idSODTM1S=idvar
            CASE ('idDD0TM1S')
              idDD0TM1S=idvar
            CASE ('idBFORMAXS')
              idBFORMAXS=idvar
            CASE ('idISWBNTHS')
              idISWBNTHS=idvar
!     
!  sediment flux  model parameters
!
            CASE ('idHSED')
              idHSED=idvar
            CASE ('idVSED')
              idVSED=idvar
            CASE ('idVPMIX')
              idVPMIX=idvar
            CASE ('idVDMIX')
              idVDMIX=idvar
            CASE ('idM1')
              idM1=idvar
            CASE ('idM2')
              idM2=idvar
            CASE ('idpie1si')
              idpie1si=idvar
            CASE ('idpie2si')
              idpie2si=idvar
            CASE ('idpie1po4m')
              idpie1po4m=idvar
            CASE ('idpie1po4n')
              idpie1po4n=idvar
!
!  sediment flux  model output
!
            CASE ('idsed_time')
              idsed_time=idvar
            CASE ('idctemp')
              idctemp=idvar
            CASE ('idpop1r')
              idpop1r=idvar
            CASE ('idpop2r')
              idpop2r=idvar
            CASE ('idpop3r')
              idpop3r=idvar
            CASE ('idpopr')
              idpopr=idvar
            CASE ('idpon1r')
              idpon1r=idvar
            CASE ('idpon2r')
              idpon2r=idvar
            CASE ('idpon3r')
              idpon3r=idvar
            CASE ('idponr')
              idponr=idvar
            CASE ('idpoc1r')
              idpoc1r=idvar
            CASE ('idpoc2r')
              idpoc2r=idvar
            CASE ('idpoc3r')
              idpoc3r=idvar
            CASE ('idpocr')
              idpocr=idvar
            CASE ('idpo4t2r')
              idpo4t2r=idvar
            CASE ('idhst2r')
              idhst2r=idvar
            CASE ('idsit2r')
              idsit2r=idvar
            CASE ('idpsiavr')
              idpsiavr=idvar
            CASE ('idflxpop')
              idflxpop=idvar
            CASE ('idflxpon')
              idflxpon=idvar
            CASE ('idflxpoc')
              idflxpoc=idvar
            CASE ('ido20')
              ido20=idvar
            CASE ('idcsod')
              idcsod=idvar
            CASE ('idsod')
              idsod=idvar
            CASE ('ids')
              ids=idvar
            CASE ('idjp')
              idjp=idvar
            CASE ('idjn')
              idjn=idvar
            CASE ('idjc')
              idjc=idvar
            CASE ('idjo2nh4')
              idjo2nh4=idvar
            CASE ('idxjco2av')
              idxjco2av=idvar
            CASE ('idxjc1av')
              idxjc1av=idvar
            CASE ('idjpo4')
              idjpo4=idvar
            CASE ('idjnh4')
              idjnh4=idvar
            CASE ('idjno3')
              idjno3=idvar
            CASE ('idjhs')
              idjhs=idvar
            CASE ('idjsi')
              idjsi=idvar
            CASE ('idjch4aq')
              idjch4aq=idvar
            CASE ('idjch4g')
              idjch4g=idvar
            CASE ('idh1')
              idh1=idvar
            CASE ('idpo40')
              idpo40=idvar
            CASE ('idpo41')
              idpo41=idvar
            CASE ('idpo42')
              idpo42=idvar
            CASE ('idpo4t2')
              idpo4t2=idvar
            CASE ('idnh40')
              idnh40=idvar
            CASE ('idnh41')
              idnh41=idvar
            CASE ('idnh42')
              idnh42=idvar
            CASE ('idnh4t2')
              idnh4t2=idvar
            CASE ('idno30')
              idno30=idvar
            CASE ('idno31')
              idno31=idvar
            CASE ('idno32')
              idno32=idvar
            CASE ('idno3t2')
              idno3t2=idvar
            CASE ('idhs1')
              idhs1=idvar
            CASE ('idhs2')
              idhs2=idvar
            CASE ('idhst2')
              idhst2=idvar
            CASE ('idsi0')
              idsi0=idvar
            CASE ('idsi1')
              idsi1=idvar
            CASE ('idsi2')
              idsi2=idvar
            CASE ('idsit2')
              idsit2=idvar
            CASE ('idbnthden')
              idbnthden=idvar
            CASE ('idch41')
              idch41=idvar
            CASE ('idch42')
              idch42=idvar
            CASE ('idch4t2')
              idch4t2=idvar
            CASE ('idso41')
              idso41=idvar
            CASE ('idso42')
              idso42=idvar
            CASE ('idso4t2')
              idso4t2=idvar
            CASE ('idflxpop1')
              idflxpop1=idvar
            CASE ('idflxpop2')
              idflxpop2=idvar
            CASE ('idflxpop3')
              idflxpop3=idvar
            CASE ('idflxpon1')
              idflxpon1=idvar
            CASE ('idflxpon2')
              idflxpon2=idvar
            CASE ('idflxpon3')
              idflxpon3=idvar
            CASE ('idflxpoc1')
              idflxpoc1=idvar
            CASE ('idflxpoc2')
              idflxpoc2=idvar
            CASE ('idflxpoc3')
              idflxpoc3=idvar
            CASE ('idflxpos')
              idflxpos=idvar
            CASE ('idflxpo4t2s')
              idflxpo4t2s=idvar
            CASE ('idflxsit2s')
              idflxsit2s=idvar

            ! Global dump averaging is performed
            CASE ('idsedG_TIME')
              idsedG_TIME=idvar
            CASE ('idCTEMP_GDA')
              idCTEMP_GDA=idvar
            CASE ('idPOP1R_GDA')
              idPOP1R_GDA=idvar
            CASE ('idPOP2R_GDA')
              idPOP2R_GDA=idvar
            CASE ('idPOP3R_GDA')
              idPOP3R_GDA=idvar
            CASE ('idPOPR_GDA')
              idPOPR_GDA=idvar
            CASE ('idPON1R_GDA')
              idPON1R_GDA=idvar
            CASE ('idPON2R_GDA')
              idPON2R_GDA=idvar
            CASE ('idPON3R_GDA')
              idPON3R_GDA=idvar
            CASE ('idPONR_GDA')
              idPONR_GDA=idvar
            CASE ('idPOC1R_GDA')
              idPOC1R_GDA=idvar
            CASE ('idPOC2R_GDA')
              idPOC2R_GDA=idvar
            CASE ('idPOC3R_GDA')
              idPOC3R_GDA=idvar
            CASE ('idPOCR_GDA')
              idPOCR_GDA=idvar
            CASE ('idPO4T2R_GDA')
              idPO4T2R_GDA=idvar
            CASE ('idHST2R_GDA')
              idHST2R_GDA=idvar
            CASE ('idSIT2R_GDA')
              idSIT2R_GDA=idvar
            CASE ('idPSIAVR_GDA')
              idPSIAVR_GDA=idvar
            CASE ('idFLXPOP_GDA')
              idFLXPOP_GDA=idvar
            CASE ('idFLXPON_GDA')
              idFLXPON_GDA=idvar
            CASE ('idFLXPOC_GDA')
              idFLXPOC_GDA=idvar
            CASE ('idO20_GDA')
              idO20_GDA=idvar
            CASE ('idCSOD_GDA')
              idCSOD_GDA=idvar
            CASE ('idSOD_GDA')
              idSOD_GDA=idvar
            CASE ('idS_GDA')
              idS_GDA=idvar
            CASE ('idJP_GDA')
              idJP_GDA=idvar
            CASE ('idJN_GDA')
              idJN_GDA=idvar
            CASE ('idJC_GDA')
              idJC_GDA=idvar
            CASE ('idJO2NH4_GDA')
              idJO2NH4_GDA=idvar
            CASE ('idXJCO2AV_GDA')
              idXJCO2AV_GDA=idvar
            CASE ('idXJC1AV_GDA')
              idXJC1AV_GDA=idvar
            CASE ('idJPO4_GDA')
              idJPO4_GDA=idvar
            CASE ('idJNH4_GDA')
              idJNH4_GDA=idvar
            CASE ('idJNO3_GDA')
              idJNO3_GDA=idvar
            CASE ('idJH2S_GDA')
              idJH2S_GDA=idvar
            CASE ('idJSI_GDA')
              idJSI_GDA=idvar
            CASE ('idJCH4AQ_GDA')
              idJCH4AQ_GDA=idvar
            CASE ('idJCH4G_GDA')
              idJCH4G_GDA=idvar
            CASE ('idH1_GDA')
              idH1_GDA=idvar
            CASE ('idPO40_GDA')
              idPO40_GDA=idvar
            CASE ('idPO41_GDA')
              idPO41_GDA=idvar
            CASE ('idPO42_GDA')
              idPO42_GDA=idvar
            CASE ('idPO4T2_GDA')
              idPO4T2_GDA=idvar
            CASE ('idNH40_GDA')
              idNH40_GDA=idvar
            CASE ('idNH41_GDA')
              idNH41_GDA=idvar
            CASE ('idNH42_GDA')
              idNH42_GDA=idvar
            CASE ('idNH4T2_GDA')
              idNH4T2_GDA=idvar
            CASE ('idNO30_GDA')
              idNO30_GDA=idvar
            CASE ('idNO31_GDA')
              idNO31_GDA=idvar
            CASE ('idNO32_GDA')
              idNO32_GDA=idvar
            CASE ('idNO3T2_GDA')
              idNO3T2_GDA=idvar
            CASE ('idHS1_GDA')
              idHS1_GDA=idvar
            CASE ('idHS2_GDA')
              idHS2_GDA=idvar
            CASE ('idHST2_GDA')
              idHST2_GDA=idvar
            CASE ('idSI0_GDA')
              idSI0_GDA=idvar
            CASE ('idSI1_GDA')
              idSI1_GDA=idvar
            CASE ('idSI2_GDA')
              idSI2_GDA=idvar
            CASE ('idSIT2_GDA')
              idSIT2_GDA=idvar
            CASE ('idJN2_GDA')
              idJN2_GDA=idvar
            CASE ('idCH41_GDA')
              idCH41_GDA=idvar
            CASE ('idCH42_GDA')
              idCH42_GDA=idvar
            CASE ('idCH4T2_GDA')
              idCH4T2_GDA=idvar
            CASE ('idSO41_GDA')
              idSO41_GDA=idvar
            CASE ('idSO42_GDA')
              idSO42_GDA=idvar
            CASE ('idSO4T2_GDA')
              idSO4T2_GDA=idvar
            CASE ('idFLXPOP1_GDA')
              idFLXPOP1_GDA=idvar
            CASE ('idFLXPOP2_GDA')
              idFLXPOP2_GDA=idvar
            CASE ('idFLXPOP3_GDA')
              idFLXPOP3_GDA=idvar
            CASE ('idFLXPON1_GDA')
              idFLXPON1_GDA=idvar
            CASE ('idFLXPON2_GDA')
              idFLXPON2_GDA=idvar
            CASE ('idFLXPON3_GDA')
              idFLXPON3_GDA=idvar
            CASE ('idFLXPOC1_GDA')
              idFLXPOC1_GDA=idvar
            CASE ('idFLXPOC2_GDA')
              idFLXPOC2_GDA=idvar
            CASE ('idFLXPOC3_GDA')
              idFLXPOC3_GDA=idvar
            CASE ('idFLXPOS_GDA')
              idFLXPOS_GDA=idvar
            CASE ('idFLXPO4T2S_GDA')
              idFLXPO4T2S_GDA=idvar
            CASE ('idFLXSIT2S_GDA')
              idFLXSIT2S_GDA=idvar
          END SELECT
       ENDIF
       GOTO 30
  40   WRITE (OUT,70) TRIM(ADJUSTL(Vinfo(1,idvar)))
       CALL EXIT
  50   CLOSE (31)

  60    FORMAT (/,' ncset_ncparams - Unable to open variable 
     .             information',
     .          ' file: ',/,15x,a,/,15x,'Default file is located in',
     .          ' Master_Code/netcdf directory.')
  70    FORMAT (/,' ncset_ncparams - Error while reading information ',
     .          'for variable: ',a)
        RETURN
      END SUBROUTINE ncset_ncparams
