C     PROGRAM RCA
C 
C 
C*******************************************************************************
C 
C 
C                  RCA  -  ROW/COLUMN VERSION OF AESOP 
C
C                                                      (RELEASE 3.0)
C 
C                                     DEVELOPED BY: HYDROQUAL INC.
C                                                   MAHWAH, N.J.
C                                     PROGRAMMERS: JAMES FITZPATRICK
C                                                  KAI-YUAN YANG
C 
C*******************************************************************************
C 
C*               Copyright (C) 1980, HydroQual Inc.                       
C*                                                                        
C*  These coded instructions, statements, and computer programs contain   
C*  unpublished proprietary information of HydroQual, Inc., and are     
C*  protected by Federal copyright law.  No part of this program and/or
C*  modifications of this program may be reproduced in any form by print
C*  photocopy, microfilm, magnetic- or optical-media or by any other
C*  means, without written permission of HydroQual, Inc.
C*	Authors: James Fitzpatrick
C*               Dominic DiToro
C*               Kai-Yuan Yang
C*               HydroQual, Inc.
C*               1200 Mac Arthur Blvd.
C*               Mahwah, NJ  07430
C*               USA
C 
C*******************************************************************************
C 
C 
C        THE FOLLOWING VARIABLES MUST BE DEFINED FOR EACH SITE SPECIFIC
C        WATER QUALITY MODEL VIA THE INCLUDE BLOCK -RCACM-
C
C     NX          -  NUMBER OF CELLS IN THE X-DIRECTION IN THE MODEL DOMAIN
C     NY          -  NUMBER OF CELLS IN THE Y-DIRECTION IN THE MODEL DOMAIN
C     NZ          -  NUMBER OF DEPTH LAYERS IN THE MODEL DOMAIN
C     NOSYS       -  NUMBER OF SYSTEMS OR WATER QUALITY VARIABLES
C                     IN THE MODEL
C
C        THE FOLLOWING VARIABLES SPECIFY THE MAXIMUM DIMENSIONS FOR
C        KEY MODEL VARIABLES
C     MXHYDFILES  -  MAXIMUM NUMBER OF LINKED HYDRODYNAMIC FILES
C     MXBC        -  MAXIMUM NUMBER OF BOUNDARY CONDITIONS PER SYSTEM
C     MXWK        -  MAXIMUM NUMBER OF FORCING FUNCTIONS OR LOADS
C     MXPARM2D    -  MAXIMUM NUMBER OF DEPTH-INDEPENDENT PARAMETERS
C     MXPARM3D    -  MAXIMUM NUMBER OF DEPTH-DEPENDENT PARAMETERS
C     MXCONS      -  MAXIMUM NUMBER OF CONSTANTS 
C     MXFUNC      -  MAXIMUM NUMBER OF MISCELLANEOUS TIME FUNCTIONS 
C     MXFUNCT     -  MAXIMUM NUMBER OF TIME BREAKS FOR THE
C                    MISCELLANEOUS TIME FUNCTIONS
C
C
C        DESCRIPTION OF COMMON VARIABLES
C 
C     IN          -  INPUT DEVICE NUMBER 
C     OUT         -  OUTPUT DEVICE NUMBER
C     MXACTS      -  NUMBER OF THE MAXIMUM (OR LAST) ACTIVE SYSTEM
C     ISYS        -  SYSTEM CURRENTLY HAVING ITS DERIVATIVES CALCULATED
C     IX,IY,IZ    -  MODEL SEGMENT CURRENTLY BEING ANALYZED
C     LIST(5)     -  OPTIONS TO LIST COMPONENTS OF MODEL INPUT DECK 
C                   = 0 - DO NOT LIST INPUT
C                   = 1 - PROVIDE USER WITH LIST OF INPUT
C                   LIST(1) - GEOMETRY
C                   LIST(2) - Boundary CONDITIONS
C                   LIST(3) - LOADS
C                   LIST(4) - PARAMETERS, CONSTANTS AND TIME FUNCTIONS
C                   LIST(5) - INITIAL CONDITIONS
C 
C 
C     IDISK       -  FLAG WHICH SIGNALS -TUNER- TO INITIATE DISK STORAGE 
C                    OF VARIABLES THAT USER WISHES TO DUMP AT END OF RUN 
C     IPRNTGSECS  -  PRINT INTERVAL (IN SECONDS) REQUESTED FOR GLOBAL
C                    DUMPS 
C     NXPRTG      -  NEXT PRINT TIME (IN SECONDS) REQUESTED FOR GLOBAL
C                    DUMPS
C     IGDOPT      -  GLOBAL DUMP TIME-AVERAGING OPTION
C                    = 0 - NO TIME-AVERAGING, I.E. INSTANTANEOUS VALUE
C                    = 1 - PERFORM TIME-AVERAGING WILL BE PERFORMED
C     IAVGGDCNTR  -  COUNTER TO KEEP TRACK OF NUMBER OF INTEGRATION STEPS GOING
C                    INTO GLOBAL DUMP AVERAGE
C     IDUMP       -  ARRAY CONTAINING THE SEGMENTS NUMBERS TO BE PRINTED AS
C                    PART OF THE INTERMEDIATE DUMPS, I.E. THE SIX SEGMENTS
C                    WHICH ARE DUMPED TO STANDOUT OUT DURING THE COURSE OF
C                    THE MODEL SIMULATION (PRINTED AT INTERVAL 'IPRNTGSECS')
C     IPRNTDSECS  -  PRINT INTERVAL (IN SECONDS) USED TO GENERATE THE
C                    'DETAILED' DUMPS FOR SELECTED SEGMENTS (VIA RCAWRIT)
C     NXPRTD      -  NEXT PRINT TIME (IN SECONDS) USED TO GENERATE THE
C                    'DETAILED' DUMPS FOR SELECTED SEGMENTS (VIA RCAWRIT)
C     IDDOPT      -  DETAILED DUMP TIME-AVERAGING OPTION
C                    = 0 - NO TIME-AVERAGING, I.E. INSTANTANEOUS VALUE
C                    = 1 - PERFORM TIME-AVERAGING WILL BE PERFORMED
C     IAVGDDCNTR  -  COUNTER TO KEEP TRACK OF NUMBER OF INTEGRATION STEPS GOING
C                    INTO DETAILED DUMP AVERAGE
C     NDMPS       -  NUMBER OF SEGMENTS SPECIFIED FOR DETAILED OUTPUT
C                    DUMPS AT PRINT INTERVAL 'IPRNTDSECS'
C     IFDMPS      -  ARRAY CONTAINING LIST OF SEGMENT NUMBERS FOR WHICH
C                    'DETAILED' DUMPS, I.E. CONCENTRATIONS AND OTHER SECONDARY 
C                    VARIABLES OF INTEREST ARE TO BE SAVED VIA RCAWRIT
C     IREC        -  COUNTER USED BY RCA FOR GLOBAL DUMP DISK FILE STORAGE 
C
C     MASSBAL     -  MASS BALANCE/FLUX BALANCE OPTION
C                    = 0 - NO BALANCES REQUESTED
C                    = 1 - PERFORM MASS/FLUX BALANCE COMPUTATIONS
C     IPRNTMBSECS -  PRINT INTERVAL (IN SECONDS) FOR OUTPUTING BALANCE
C                    COMPUTATIONS
C     NXPRTMB     -  NEXT TIME (IN SECONDS) TO OUTPUT THE BALANCE
C                    COMPUTATIONS
C     IMBOPT      -  MASS/FLUX BALANCE TIME-AVERAGING OPTION
C                    = 0 - NO TIME-AVERAGING, I.E. INSTANTANEOUS VALUE
C                    = 1 - PERFORM TIME-AVERAGING WILL BE PERFORMED
C     IAVGMBCNTR  -  COUNTER TO KEEP TRACK OF NUMBER OF INTEGRATION STEPS GOING
C                    INTO THE MASS/FLUX BALANCE AVERAGE

C
C
C     SYNAME      -  ALPHA-NUMERIC NAMES FOR EACH SYSTEM 
C     SYSBY       -  VECTOR CONTAINING BYPASS OPTION FOR EACH SYSTEM 
C     RBY         -  VECTOR CONTAINING R BYPASS INDICATOR FOR EACH SYSTEM
C     QBY         -  VECTOR CONTAINING Q BYPASS INDICATOR FOR EACH SYSTEM
C     CYCLE       -  FLAG TO DETERMINE WHETHER INITIAL CONDITIONS ARE TO BE
C                    READ FROM A USER PROVIDED INPUT FILE OR FROM A PREVIOUS
C                    RCA RUN (I.E. THE EQUIVALENT OF A HOT -START); MAY ALSO BE
C                    USED TO LINK A SERIES OF MULTI -YEAR SIMULATIONS TOGETHER
C                    = 0 - GET INITIAL CONDITIONS FROM USER SUPPLIED INPUT
C                            FILE (EQUIVALENT TO A COLD -START RUN)
C                    = 1 - GET INITIAL CONDITIONS FROM A PREVIOUS RCA RUN -
C                            THE RCAFIC FILE (EQUIVALENT TO A HOT -START RUN)
C     NEGSLN      -  VARIABLE INDICATING WHETHER OR NOT NEGATIVE SOLUTIONS 
C                    WILL BE PERMITTED 
C                    = 0 - NEGATIVE SOLUTIONS NOT PERMITTED  (FLOOR = 0.0) 
C                    = 1 - NEGATIVE SOLUTIONS PERMITTED  (FLOOR = -1.E 35) 
C     INTGRTYP    -  INTEGRATION METHOD
C                    = 1 - EXPLICIT - UPWIND
C                    = 3 - SPLIT TIME STEP - UPWIND
C                    = 4 - EXPLICIT - UPWIND/SMOLARKIEWICZ
C                    = 5 - LEAPFROG - UPWIND/SMOLARKIEWICZ
C                    = 6 - SPLIT TIME STEP - UPWIND/SMOLARKIEWICZ
C     ISMOLAR     -  SMOLARKIEWICZ CORRECTOR OPTION
C                    = 0 - ECOM3D SMOLAR_2 - SECOND ORDER ACCURATE SCHEME
C                    = 1 - ECOM3D SMOLAR_R - RECURSIVE SCHEME
C     ISMOLBCOPT  -  BOUNDARY SMOLARKIEWICZ CORRECTOR OPTION
C                    = 0 - DO NOT APPLY SMOLARKIEWICZ CORRECTIONS AT
C                          BOUNDARY CONDITIONS
C                    = 1 - APPLY SMOLARKIEWICZ CORRECTIONS AT BOUNDARY
C                          CONDITIONS
C     IDIAGDT     -  FLAG USED TO SELECT WHETHER A FULL SIMULATION IS TO BE
C                    PERFORMED OR WHETHER A DIAGNOSTIC ANALYSIS OF INTEGRATION
C                    STEPSIZES IS TO BE PERFORMED
C                    = 0 - PERFORM MODEL SIMULATION
C                    = 1 - DETERMINE CRITICAL INTEGRATION STEPSIZES
C     INPCHCK     -  FLAG USED TO SELECT WHETHER A FULL SIMULATION IS TO BE
C                    PERFORMED OR WHETHER A CHECK OF THE USER'S INPUT DECK
C                    FOR ERRORS IS TO BE PERFORMED
C                    = 0 - PERFORM MODEL SIMULATION
C                    = 1 - CHECK INPUT DECK FOR ERRORS
C     CARAY       -  CONCENTRATIONS FOR APPROPRIATE SYSTEM AND SEGMENT 
C     CDARAY      -  CONCENTRATION DERIVATIVES
C     CMAX        -  MAXIMUM ALLOWABLE CONCENTRATION FOR EACH SYSTEM 
C                    IF ANY OF THESE VALUES ARE EXCEEDED IN THEIR RESPECTIVE  
C                    SYSTEMS THE PROGRAM WILL ABORT
C     CMIN        -  MINIMUM CONCENTRATIONS, FOR EACH SYSTEM, BELOW WHICH
C                    THE ACCURACY CRITERIA IS IGNORED
C
C   /TIMEINFO/  VARIABLES ASSOCIATED WITH SIMULATION TIME AND
C               INTEGRATION HISTORY
C     INITB       -  FLAG WHICH PERMITS THE USER TO PROVIDE INITIALIZATION 
C                    OF VARIABLES IN KINETIC SUBROUTINE    - TUNER-
C     IDTSECS     -  INTEGRATION INTERVAL - SECONDS
C     ITIMESECS   -  INTERNAL CLOCK USED TO KEEP TRACK OF THE
C                    SIMULATION TIME  - SECONDS
C     IDTWQSECS   -  INTEGRATION STEPSIZE TO BE USED FOR EVALUATING THE
C                    KINETIC PORTION OF THE MASS BALANCE EQUATIONS -
C                    SECONDS
C     ITIMEWQSECS -  INTERNAL CLOCK USED TO KEEP TRACK OF KINETIC PORTION OF
C                    THE INTEGRATION PROCEDURE - SECONDS
C     TZERO       -  TIME USER WISHES TO START SIMULATION AT, OR TIME FOR
C                    WHICH THE USER WANTS A STEADY STATE SOLUTION
C     DT          -  INTEGRATION INTERVAL - DAYS
C     TIME        -  INTERNAL CLOCK USED TO KEEP TRACK OF THE
C                    SIMULATION TIME - DAYS
C     DTWQ        -  INTEGRATION STEPSIZE TO BE USED FOR EVALUATING THE
C                    KINETIC PORTION OF THE MASS BALANCE EQUATIONS - DAYS
C     TIMEWQ      -  INTERNAL CLOCK USED TO KEEP TRACK OF KINETIC PORTION OF
C                    THE INTEGRATION PROCEDURE - DAYS
C     ISCALT      -  TIME WARP FACTOR ... TIME CONVERSION FACTOR - INTEGER 
C     SCALT       -  TIME WARP FACTOR ... TIME CONVERSION FACTOR - REAL 
C     TEND        -  ENDING TIME OF SIMULATION (OR TIME TO CHANGE 
C                    TO NEXT INTEGRATION INTERVAL) - DAYS
C
C
C   /INTHIST/  VARIABLES ASSOCIATED WITH INTEGRATION HISTORY TIMESTEPS
C     NSTEP       -  NUMBER OF INTEGRATION STEPSIZES/TIMES TO BE USED
C                    FOR DESCRIBING THE INTEGRATION HISTORY
C     ISTEP       -  VECTOR CONTAINING THE INTEGRATION STEPSIZES - DT (SECS)
C     TBRK        -  VECTOR CONTAINING THE CORRESPONDING TIMES FOR THE VECTOR
C                    -ISTEP-. ISTEP(I) WILL BE USED UNTIL TIME TBRK(I), THEN
C                    SWITCH TO ISTEP(I+1) UNTIL TBRK(I+1), ETC.
C
C
C   /FILENAMES/  VARIABLES ASSOCIATED WITH READING MODEL INPUTS
C     IHYDFILE    -  COUNTER TO KEEP TRACK OF WHICH HYDRODYNAMIC FILE
C                    TO OPEN AND READ
C     HYDFILNA    -  VECTOR CONTAINING NAMES OF THE HYDRODYNAMIC FILES
C     DIFFILNA    -  VECTOR CONTAINING NAMES OF THE DIFFUSER FILES
C     BCFILNA     -  NAME OF THE FILE CONTAINING BOUNDARY CONDITIONS
C     PSFILNA     -  NAME OF THE FILE CONTAINING POINT SOURCE LOADS
C     NPSFILNA    -  NAME OF THE FILE CONTAINING NONPOINT SOURCE LOADS
C     FLFILNA     -  NAME OF THE FILE CONTAINING FALL-LINE LOADS
C     ATMFILNA    -  NAME OF THE FILE CONTAINING ATMOSPHERIC LOADS
C     PCFILNA     -  NAME OF THE FILE CONTAINING PARAMETERS, CONSTANTS,
C                    MISCELLANEOUS TIME FUNCTIONS AND FILE NAMES FOR
C                    MISCELLANEOUS INPUT FILES REQUIRED BY A USER'S
C                    KINETIC SUBROUTINE
C     KINFILNA    -  NAMES OF THE MISCELLANEOUS INPUT FILES TO BE USED
C                    FOR A USER'S KINETIC SUBROUTINE (EX. INPUT REQ'D
C                    WHEN USING THE SEDIMENT NUTRIENT FLUX MODEL AS PART
C                    OF A EUTROPHICATION ANALYSIS)
C     ICFILNA     -  NAME OF THE FILE CONTAINING INITIAL CONDITIONS
C     IBNRYRDOPTS -  VECTOR SPECIFYING WHETHER ASCII OR BINARY FILES
C                    WILL BE READ FOR BOUNDARY CONDITIONS, LOADS, AND
C                    INITIAL CONDITIONS
C                    = 0 - ASCII FILE
C                    = 1 - BINARY FILE
C 
C 
C   /UNITS/  VARIABLES ASSOCIATED WITH "SCALING" OF TIME FOR READING
C            INPUT VARIABLES
C     TWARP       -  INTEGRATION TIME
C     TWARPP      -  PRINT INTERVALS
C     TWARPBC     -  BOUNDARY CONDITION
C     TWARPPS     -  POINT SOURCE LOADS
C     TWARPNPS    -  NONPOINT SOURCE LOADS
C     TWARPFL     -  FALL-LINE LOADS
C     TWARPATM    -  ATMOSPHERIC LOADS
C     TWARPTVF    -  MISCELLANEOUS TIME FUNCTIONS
C     ISCALBC     -  BOUNDARY CONDITION
C     ISCALPS     -  POINT SOURCE LOADS
C     ISCALNPS    -  NONPOINT SOURCE LOADS
C     ISCALFL     -  FALL-LINE LOADS
C     ISCALATM    -  ATMOSPHERIC LOADS
C     ISCALTVF    -  MISCELLANEOUS TIME FUNCTIONS
C 
C 
C   /VOLUMES/  MODEL VOLUMES AND VARIABLES ASSOCIATED WITH READING
C              HYDRODYNAMIC FILES
C     IHYDDTSECS  -  INTERVAL (IN SECONDS) BETWEEN HYDRODYNAMIC UPDATES
C                    (I.E., THE PERIOD USED TO AVERAGE TRANSPORT
C                    COMPUTATIONS IN THE HYDRODYNAMIC MODEL)
C     HYDDT       -  INTERVAL (IN DAYS) BETWEEN HYDRODYNAMIC UPDATES
C     NXHYDTSECS  -  TIME (IN SECONDS) AT WHICH THE NEXT READ OF THE
C                    HYDRODYNAMIC FILE IS TO OCCUR
C     NXHYDTSECS  -  TIME (IN DAYS) AT WHICH THE NEXT READ OF THE
C                    HYDRODYNAMIC FILE IS TO OCCUR
C     HYDCYCOPT   -  FLAG TO SPECIFY WHETHER THE HYDRODYNAMIC INPUT FILE
C                    IS TO CYCLED (I.E., RE-READ)
C                    = 1 - HYDRODYNAMIC FILE IS TO BE CYCLED (SAME AS
C                          ASSUMMING A PERIODIC STEADY STATE FOR HYDRO)
C                    = 2 - NO CYCLING
C     BVOL        -  VECTOR CONTAINING VOLUMES (MCM)
C     VDER        -  VOLUME DERIVATIVES (MCM/DAY)

C   /PARCONFNS/  VARIABLES ASSOCIATED WITH PARAMETERS, CONSTANTS, AND
C                MISCELLANEOUS TIME FUNCTIONS
C     PARAM2D     -  2-D ARRAY CONTAINING PARAMETERS FOR EACH LAYER
C     PARAM3D     -  3-D ARRAY CONTAINING PARAMETERS FOR EACH SEGMENT
C     NOCONS      -  NUMBER OF CONSTANTS 
C     CONST       -  A VECTOR OF CONSTANTS COMMON ACROSS ALL SEGMENTS
C     NOFUNC      -  NUMBER OF MISCELLANEOUS TIME FUNCTIONS 
C     ITVFPWLOPT  -  FLAG SPECIFYING WHETHER THE TIME FUNCTIONS WILL BE
C                    EVALUATED AS STEP FUNCTIONS OR PIECEWISE LINEAR
C                    FUNCTIONS
C                    = 0 - STEP FUNCTIONS
C                    = 1 - PIECEWISE LINEAR INTERPOLATION
C     NXFUNTSECS  -  NEXT TIME (IN SECONDS) TO UPDATE TIME FUNCTION
C     NXFUNT      -  NEXT TIME (IN DAYS) TO UPDATE TIME FUNCTION
C     NXCALL13T   -  NEXT TIME (IN DAYS) TO CALL RCA13 - SUBROUTINE
C                    THAT UPDATES THE TIME FUNCTIONS
C     ITIMF       -  COUNTER THAT TRACKS THE TIME BREAKS FOR THE
C                    MISCELLANEOUS TIME FUNCTIONS
C     MFUNC       -  VECTOR CONTAINING 'SLOPE' OF THE KINETIC TIME FUNCTIONS
C     BFUNC       -  VECTOR CONTAINING THE INTERCEPT OF THE KINETIC TIME
C                    FUNCTIONS 
C     VALMTF      -  VALUES OF THE MISCELLANEOUS TIME FUNCTIONS
C     TIMEMTF     -  TIME BREAKS FOR THE MISCELLANEOUS TIME FUNCTIONS
C
C 
C   /BOUNDCOND/  VARIABLES ASSOCIATED WITH THE BOUNDARY CONDITIONS
C     IBCOPT      -  BC OPTION
C                    = 1 - SIGMA LEVEL - CONSTANT
C                    = 2 - SIGMA LEVEL - TIME-VARYING
C                    = 3 - STANDARD LEVEL - CONSTANT
C                    = 4 - STANDARD LEVEL - TIME-VARYING
C     IBCPWLOPT   -  PIECEWISE LINEAR BOUNDARY CONCENTRATION OPTION
C                    = 0 - STEP FUNCTION APPROXIMATION FOR BCS
C                    = 1 - PIECEWISE LINEAR APPROXIMATION FOR BCS
C     NOBC        -  NUMBER OF B.C.'S FOR EACH SYSTEM
C     NXBCTSECS   -  NEXT TIME (IN SECONDS) TO UPDATE BOUNDARY CONDITIONS
C     NXBCT       -  NEXT TIME (IN DAYS) TO UPDATE BOUNDARY CONDITIONS
C     SCALBC      -  SCALE FACTORS FOR BOUNDARY CONDITIONS
C     IBC         -  VECTOR CONTAINING THE (I,J,K) LOCATIONS OF THE
C                    BOUNDARY CONDITIONS
C     BBC         -  VECTOR CONTAINING BOUNDARY CONDITIONS (MG/L)
C     SBC         -  VECTOR CONTAINING 'SLOPE' OF BOUNDARY CONDITIONS
C                    PER SYSTEM
C     SLDEPTH     -  STANDARD LEVEL DEPTHS (IF BCS TO BE READ USING
C                    STANDARD LEVEL OPTION)
C     NOBCSL      -  NUMBER OF STANDARD LEVELS FOR EACH BC
C
C
C   /FORCFNS/  VARIABLES ASSOCIATED WITH THE FORCING FUNCTIONS (I.E., THE LOADS)
C     IPSOPT      -  PS LOADING OPTION (1==>CONSTANT;2==>TIME-VARYING)
C     IPSPWLOPT   -  FLAG TO SPECIFY WHETHER TIME VARIABLE POINT SOURCE
C                    LOADS ARE STEP-FUNCTIONS OR PIECEWISE LINEAR
C     NOPS        -  NUMBER OF PS LOADS FOR EACH SYSTEM
C     IPSTABL     -  INDEX TABLE OF POINT SOURCE LOADS
C     ZFRACPS     -  VERTICAL DISTRIBUTION FOR EACH POINT SOURCE LOAD
C                    CONTAINED IN -IPSTABL-
C     NXPSTSECS   -  NEXT TIME (IN SECONDS) TO UPDATE POINT SOURCE LOADS
C     NXPST       -  NEXT TIME (IN DAYS) TO UPDATE POINT SOURCE LOADS
C     SCALPS      -  SCALE FACTORS FOR POINT SOURCE LOADS
C     IPS         -  VECTOR CONTAINING (I,J,K) LOCATIONS OF THE POINT
C                    SOURCE LOADS
C     BPS         -  VECTOR CONTAINING POINT SOURCE LOADS (GM/DAY)
C     SPS         -  VECTOR CONTAINING 'SLOPE' OF POINT SOURCE LOADS
C     INPSOPT     -  NPS OPTION (1==>CONSTANT;2==>TIME-VARYING)
C     INPSPWLOPT  -  FLAG TO SPECIFY WHETHER TIME VARIABLE NONPOINT SOURCE
C                    LOADS ARE STEP-FUNCTIONS OR PIECEWISE LINEAR
C     NONPS       -  NUMBER OF NPS LOADS FOR EACH SYSTEM
C     INPSTABL    -  INDEX TABLE OF NONPOINT SOURCE LOADS
C     ZFRACNPS    -  VERTICAL DISTRIBUTION FOR EACH NONPOINT SOURCE LOAD
C                    CONTAINED IN -INPSTABL-
C     NXNPSTSECS  -  NEXT TIME (IN SECONDS) TO UPDATE NONPOINT SOURCE LOADS
C     NXNPST      -  NEXT TIME (IN DAYS) TO UPDATE NONPOINT SOURCE LOADS
C     SCALNPS     -  SCALE FACTORS FOR NONPOINT SOURCE LOADS
C     INPS        -  VECTOR CONTAINING (I,J,K) LOCATIONS OF THE NONPOINT
C                    SOURCE LOADS
C     BNPS        -  VECTOR CONTAINING NONPOINT SOURCE LOADS (GM/DAY)
C     SNPS        -  VECTOR CONTAINING 'SLOPE' OF NONPOINT SOURCE LOADS
C     IFLPWLOPT   -  FLAG TO SPECIFY WHETHER TIME VARIABLE FALL-LINE
C                    LOADS ARE STEP-FUNCTIONS OR PIECEWISE LINEAR
C     NOFL        -  NUMBER OF FL LOADS FOR EACH SYSTEM
C     IFLTABL     -  INDEX TABLE OF FALL-LINE LOADS
C     ZFRACFL     -  VERTICAL DISTRIBUTION FOR EACH FALL-LINE LOAD
C                    CONTAINED IN -IFLTABL-
C     NXFLTSECS   -  NEXT TIME (IN SECONDS) TO UPDATE FALL-LINE LOADS
C     NXFLT       -  NEXT TIME (IN DAYS) TO UPDATE FALL-LINE LOADS
C     SCALFL      -  SCALE FACTORS FOR FALL-LINE LOADS
C     IFL         -  VECTOR CONTAINING (I,J,K) LOCATIONS OF THE FALL-LINE
C                    LOADS
C     BFL         -  VECTOR CONTAINING FALL-LINE LOADS (GM/DAY)
C     SFL         -  VECTOR CONTAINING 'SLOPE' OF FALL-LINE LOADS
C     IFLOPT      -  FALL-LINE OPTION (1==>CONSTANT;2==>TIME-VARYING)
C     IATMOPT     -  ATMOSPHERIC OPTION (1==>CONSTANT;2==>TIME-VARYING)
C     IATMPWLOPT  -  FLAG TO SPECIFY WHETHER TIME VARIABLE ATMOSPHERIC
C                    LOADS ARE STEP-FUNCTIONS OR PIECEWISE LINEAR
C     NOATM       -  NUMBER OF LOADS FOR EACH SYSTEM
C     NXATMTSECS   -  NEXT TIME (IN SECONDS) TO UPDATE ATMOSPHERIC LOADS
C     NXATMT       -  NEXT TIME (IN DAYS) TO UPDATE ATMOSPHERIC LOADS
C     SCALATM     -  SCALE FACTORS FOR ATMOSPHERIC LOADS
C     BATM        -  VECTOR CONTAINING ATMOSPHERIC LOADS (GM/DAY)
C     SBATM       -  VECTOR CONTAINING 'SLOPE' OF ATMOSPHERIC LOADS
C 
C
C   /TRANSCOEFFS/  VARIABLES ASSOCIATED WITH THE TRANSPORT COEFFICIENTS
C     IXS,IXE     -  STARTING AND ENDING LOCATIONS OF WATER CELLS FOR
C                    EACH "ROW" IN THE WATER QUALITY MODEL
C     QX,QY,QZ    -  ARRAYS CONTAINING THE ADVECTIVE FLOWS (MCM/DAY)
C     RX,RY,RZ    -  ARRAYS CONTAINING THE BULK EXCHANGE COEFFICIENTS (MCM/DAY)
C     SCALRX      -  SCALE FACTOR FOR HORIZONTAL EXCHANGE FIELD
C     SCALRY      -  SCALE FACTOR FOR LATERAL EXCHANGE FIELD
C     SCALRZ      -  SCALE FACTOR FOR VERTICAL EXCHANGE FIELD
C     AVECT       -  VECTOR CONTAINING THE OFF-DIAGONAL TERMS OF THE
C                    TRANSPORT MATRIX
C     DIAG        -  VECTOR CONTAINING THE DIAGONAL TERMS OF THE TRANSPORT
C                    MATRIX
C     EX,EY,EZ    -  ARRAYS CONTAINING THE MIXING COEFFICIENTS (MCM/DAY)
C     IDIFFOPT    -  FLAG SPECIFYING WHETHER DIFFUSER DISCHARGE FILE
C                    BEING READ
C     IDD,JDD     -  (I,J) LOCATIONS OF THE DIFFUSERS
C     VDDIST      -  VERTICAL DISTRIBUTION (FRACTIONS) OF DIFFUSER FLOWS
C                    THROUGH THE WATER COLUMN
C     QDIFF       -  DIFFUSER DISCHARGE FLOW RATES (MCM/DAY)
C 
C 
C   /GEOM/  VARIABLES ASSOCIATED WITH THE LENGTHS, WIDTHS, AND DEPTHS OF
C           THE WATER QUALITY MODEL
C     ICOLLOPT    -  FLAG WHICH SPECIFIES WHETHER A GRID-COLLAPSED
C                    VERSION OF THE HYDRODYNAMICS BEING PROVIDED TO RCA
C     LNDWTROPT   -  FLAG WHICH SPECIFIES WHETHER FULL GRID HYDRODYNAMIC
C                    TRANSPORT BEING READ OR WET-GRID ONLY
C     IWTRCNT     -  THE NUMBER OF WET-GRID CELLS
C     IWTRNDX     -  THE -I- LOCATIONS OF THE WET CELLS
C     JWTRNDX     -  THE -J- LOCATIONS OF THE WET CELLS
C     FSM         -  LAND MASK ARRAY
C     H           -  MEAN WATER DEPTH ARRAY (M)
C     DETA        -  WATER ELEVATION DERIVATIVE ARRAY (M/day)
C     ETA         -  WATER ELEVATION (RELATIVE TO MEAN WATER -H-) ARRAY (M)
C     HBAR        -  WATER DEPTH (H+ETA) ARRAY (M)
C     DZ          -  THICKNESS OF EACH SIGMA-LEVEL
C     DZZ         -  AVERAGE THICKNESS BETWEEN LEVELS
C     ZZ          -  SIGMA-LEVEL DEPTH OF A VERTICAL CELL
C     DZ          -  "LENGTHS" OF EACH WATER CELL (M)
C     DY          -  "WIDTHS" OF EACH WATER CELL (M)
C     XAX         -  CROSS-SECTIONAL AREA (ON X-INTERFACE) OF EACH WATER
C                    CELL (M^2)
C     XAY         -  CROSS-SECTIONAL AREA (ON Y-INTERFACE) OF EACH WATER
C                    CELL (M^2)
C     XAZ         -  SURFACE AREA OF EACH WATER CELL (M^2)
C 
C 
C   /SALTEMP/  SALINITY AND TEMPERATURE FROM THE HYDRO MODEL
C     HYDSAL      -  HYDRODYNAMIC MODEL COMPUTED SALINITY
C     HYDTEMP     -  HYDRODYNAMIC MODEL COMPUTED TEMPERATURE
C
C     ITIMX       -  COUNTER USED TO PICK PROPER 'B' VALUES FOR....
C     ITIMHYD     -  FOR HYDRODYNAMIC FIELDS
C     ITIMF       -  FOR MISC TIME FUNCTIONS 

      SAVE
      INCLUDE  'RCACM'
      REAL  ADDMASS(NX,NY,NZ,NOSYS)
  
C               RCA FILE STRUCTURE
C    FILE        
C     NO.  TYPE     NAME            USAGE 
C    ----  ----    ------      -------------------------
C
C      5   ASCII               STANDARD INPUT
C      6   ASCII               STANDARD OUTPUT
C
C     10   BINARY  RCAF10      PRINT TIME HISTORY FOR 'GLOBAL' DUMPS
C                              (+ NX, NY, NZ, NOSYS AS REC 1)
C                              (+ SYSTEM NAMES AS REC 2)
C                              (+ SYSBY AS REC 3)
C                              (+ FSM LAND MASK AS REC 4)
C     11   BINARY  RCAF11      CONCENTRATIONS FOR EACH MODEL SEGMENT
C     12   BINARY  RCAF12      PRINT TIME HISTORY FOR MORE DETAILED
C                              SEGMENT DUMPS
C                              (+ NOSEG DUMPED AS REC 1)
C                              (+ EQUIV TABLE OF SEGS SAVED AS REC 2)
C     13   BINARY  RCAF13      MORE DETAILED DUMPS FOR SELECTED SEGMENTS
C     14   BINARY  RCAF14      SEDIMENT LAYER DUMP FILES (IF USER MODEL
C                              INCLUDES A SEDIMENT COMPONENT)
C     15   BINARY  RCAFIC      WATER COLUMN ICs FOR MULTI-YEAR RUNS
C     16   BINARY  RCAFICSED   INITIAL CONDITIONS FOR MULTI-YEAR RUNS
C     17   BINARY  RCAFMB      WATER COLUMN MASS BALANCE/FLUX BALANCE DUMPS
C     18   BINARY  RCAFMBSED   SEDIMENT MASS BALANCE/FLUX BALANCE DUMPS
C     19   BINARY  RCAFADDM    MASS ADDED DUMPS (I.E., THE MASS ADDED TO
C                              EACH SYSTEM/SEGMENT WHEN NEGATIVE
C                              CONCENTRATIONS ENCOUTERED DURING THE
C                              INTEGRATION PROCEDURE)
C
C   20-29    *        *        CAN BE USED FOR SPECIFIC APPLICATIONS AS
C                              OUTPUT FILES AS NEEDED BY DEFINING IN THE
C                              USER-SPECIFIED KINETIC SUBROUTINE -TUNER-
C
C     30   BINARY  gcm_geom    GEOMETRY AS SUPPLIED FROM HYDRODYNAMIC MODEL
C     30   BINARY  wet_grid    (IX,IY) LOCATIONS OF THE WATER SEGMENTS
C                              IN THE HYDRODYNAMIC MODEL
C     30   BINARY  HYDFILNA    TRANSPORT TERMS AS SUPPLIED FROM HYDRODYNAMIC
C                              MODEL (THIS IS THE ECOM-3D gcm_trans FILE)
C     31   BINARY  DIFFILNA    DIFFUSER (OUTFALL) DISCHARGE VOLUMES AS SUPPLIED
C                              FROM HYDRODYNAMIC MODEL (THIS IS THE ECOM-3D
C                              gcm_qdiff FILE)
C     32   A-or-B  BCFILNA     USER SELECTED FILE FOR BOUNDARY CONDITIONS
C     33   A-or-B  PSFILNA     USER SELECTED FILE FOR POINT SOURCE LOADS
C     34   A-or-B  NPSFILNA    USER SELECTED FILE FOR NONPOINT SOURCE LOADS
C     35   A-or-B  FLFILNA     USER SELECTED FILE FOR FALL -LINE LOADS
C     36   A-or-B  ATMFILNA    USER SELECTED FILE FOR ATMOSPHERIC LOADS
C     37   A-or-B  PCFILNA     USER SELECTED FILE FOR PARAMETERS, CONSTANTS
C                              AND MISCELLANEOUS TIME FUNCTIONS
C     38   A-or-B  ICFILNA     USER SELECTED FILE FOR INITIAL CONDITIONS
C
C   39-49    *        *        CAN BE USED FOR SPECIFIC APPLICATIONS AS
C                              INPUT FILES AS NEEDED BY DEFINING IN THE
C                              USER  - SPECIFIED KINETIC SUBROUTINE -TUNER-
C
C      A-or-B = ASCII or BINARY - as specified by user selected option
C      *  AS DEFINED BY THE USER
C 

      OPEN(UNIT=10,FILE='RCAF10',FORM='UNFORMATTED')
      OPEN(UNIT=11,FILE='RCAF11',FORM='UNFORMATTED')
      OPEN(UNIT=15,FILE='RCAFIC',FORM='UNFORMATTED')

C  Open files for sediment dumps and initial conditions

      OPEN(UNIT=14,FILE='RCAF14',FORM='UNFORMATTED')
      OPEN(UNIT=16,FILE='RCAFICSED',FORM='UNFORMATTED')

C  Open and close/delete "mass added" file
      OPEN(19,FILE='RCAFADDM',FORM='UNFORMATTED')
      CLOSE(19,STATUS='DELETE')

C 
C     Calling sequence for RCA subroutines 
C 
C     Call input routines
      CALL RCA01
      CALL RCA02
      CALL RCA03
      CALL RCA04
      CALL RCA05
      CALL RCA06
      CALL RCA07
      IF(INPCHCK.EQ.1)  THEN
        WRITE(OUT,1000)
 1000   FORMAT(///
     .     25X,'CHECK OF RCA INPUT COMPLETED --- PROGRAM TERMINATED'//)
        CALL EXIT
      ENDIF

C     Call subroutines to solve state equations
C       Time-variable
      IF(INTGRTYP.EQ.1) THEN
         CALL RCAEXP                 ! Explicit - Upwind
      ELSEIF(INTGRTYP.EQ.3) THEN
         CALL RCASPLT                ! Split time step - Upwind
      ELSEIF(INTGRTYP.EQ.4) THEN
         CALL RCAEXPS                ! Explicit - Upwind, Smolarkiewicz
      ELSEIF(INTGRTYP.EQ.5) THEN
         CALL RCALFS                 ! Leapfrog - Upwind, Smolarkiewicz
      ELSEIF(INTGRTYP.EQ.6) THEN
         CALL RCASPLTS               ! Split time step - Upwind, Smolarkiewicz
      ENDIF

C     Check to see if mass added to any of the state-variables,
C        if so, inform user
      OPEN(19,FILE='RCAFADDM',FORM='UNFORMATTED',STATUS='OLD',ERR=100)
C        RCAFADDM file exists, therefore mass added
      WRITE(OUT,2000)
 2000 FORMAT(///20X,
     .  'DURING THIS SIMULATION RUN MASS WAS ADDED AS FOLLOWS...')
      READ(19)  ADDMASS
      DO 200 ISYS=1,NOSYS
        WRITE(OUT,2100)  ISYS
 2100   FORMAT(//20X,'SYSTEM',I3/
     .        13X,'IX   IY   IZ',3X,'MASS ADDED',9X,'EQUIVALENT CONC')
        TOTMASS=0.
        DO 190 IZ=1,NZ
         DO 190 IY=1,NY
          DO 190 IX=1,NX
           IF(ADDMASS(IX,IY,IZ,ISYS).EQ.0.) GO TO 190
           TOTMASS = TOTMASS+ADDMASS(IX,IY,IZ,ISYS)
           CONCMASS = ADDMASS(IX,IY,IZ,ISYS)/BVOL(IX,IY,IZ)
           WRITE(OUT,2200)  IX,IY,IZ,ADDMASS(IX,IY,IZ,ISYS)/1000.
     .                     ,CONCMASS
 2200      FORMAT(10X,3I5,E13.5,' KG',5X,E13.5,' MG/L')
  190   CONTINUE
        IF(TOTMASS.GT.0.)  WRITE(OUT,2300) ISYS,TOTMASS/1000.
 2300   FORMAT(20X,'TOTAL MASS ADDED FOR SYSTEM',I5,' =',E13.5,' KG')
  200 CONTINUE

C     Call output routines
  100 CALL RCA12

C     Close files
      CLOSE (10,STATUS='KEEP')
      CLOSE (11,STATUS='KEEP')
      CLOSE (13,STATUS='KEEP')

      CALL EXIT 
      END 