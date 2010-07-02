************************************************************************
*
*     General include file for COMBOS program
*
*     Olivier Schneider, CERN/PPE-ALE
*     December 6, 1996
*
*     Version 1.00, December  6, 1996:
*     -- original release
*     Version 1.10, December 10, 1996:
*     -- (undocumented changes)
*     -- (undocumented changes)
*     Version 1.20, December 13, 1996:
*     -- (undocumented changes)
*     Version 1.30, January 27, 1997:
*     -- added array KSYSN to handle synonyms
*     Version 1.40, February 6, 1997:
*     -- parameter MPAR increased from 30 to 40
*     -- parameter MCONT increased from 40 to 60
*     -- added array KLUMP to handle lumps
*     Version 1.50, February 12, 1997:
*     -- parameter MSTAT increased from 4 to 8
*     Version 2.00, February 20, 1997:
*     -- renamed LUNIT to LUNITS (to avoid conflict with LUNIT in master.inc)
*     -- new parameter FCONT = number of the first systematic contribution
*     -- new parameters for the switches (MSW, LADJU, LSYMM, etc ...)
*     -- new array XON (replacing, XNOAD, XNOCO_STAT, XNOCO_SYST, XNOSY, XSYMM)
*     -- new variable NDEST and array DEST (desired step values for combination)
*     -- renamed all variables associated with the combination routines; added
*         set of new variables associated with the termination routines; added
*         variables to hold logical unit numbers (LCALLC and LCALLT)
*     -- combination routine NEW_ROUTINE renamed to NEW_COMBINATION
*     -- new combination routines: COMBY_LIFETIME, DUMP_MASTER_INC 
*     -- new termination routines: SET_DMS_LIMIT, AVG_AMPLITUDE, NEW_TERMINATION
*     Version 2.10, February 27, 1997:
*     -- new EXTRAPOLATION switch (parameter LEXTR)
*     -- new termination routines KUMAC_NOSTEP and KUMAC_DMS_LIMIT
*     Version 2.11, February 27, 1997:
*     -- parameter MMETH increased from 12 to 20
*     -- parameter MNAM increased from MCONT+2 to MCONT+20
*     Version 2.20, April 3, 1997:
*     -- parameter MANAL increased from 20 to 30
*     -- parameter MEXP increased from 6 to 8
*     -- parameter MCONT increased from 60 to 80
*     -- new combination routine CHI2_SYM_CIRC_G
*     Version 2.30, May 12, 1997:
*     -- SYMMETRIZATION switch renamed SYMM_COMB (parameter LSYMM renamed LSYMC)
*     -- new SYMM_ADJU switch (parameter LSYMA)
*     -- add 4th logical unit for kumac output
*     -- termination routine KUMAC_NOSTEP removed
*     -- new combination routines BLUE and SINGLE
*     Version 2.32, Jul 16, 1997:
*     -- increase parameter MCONT from 80 to 100
*     -- new combination routine BLUE_ASYM
*     Version 2.40, October 14, 1997:
*     -- add common for MINOS errors steering (HCJS)
*     -- new combination routine CHI2_ASYM_MIN (HCJS)
*     Version 2.50, March 9, 1998:
*     -- add arrays NDISP,KDISP,DISP,XSYMD and parameters MDISP,SUPERSEDED 
*        to handle new DISPLAY logical line
*     Version 2.51, July 7, 1999:
*     -- parameter MMETH increased from 20 to 30
*     -- parameter MNAM increased from MCONT+20 to MCONT+30
*     Version 2.52, July 12, 1999:
*     -- new termination routine FINAL_FRACTIONS
*     -- introduction of the concept of preparation routines; 
*        new preparation routines ADJUST_SIGASTAT and NEW_PREPARATION
*     Version 2.53, November 27, 1999:
*     -- preparation routine  ADJUST_SIGASTAT renamed to ADJU_A_FOR_FBS
*DR   Version ???, June 1, 1999:
*     -- add parameter NQUANT,MQUANT number of quantities to be combined
*     change KCONT declaration
*     Version 2.97, July 21, 2000
*     -- parameter MCONT increased from MAX0(FCONT,MPAR+60) to
*        MAX0(FCONT,MPAR+100), i.e. from 100 to 140
*     -- parameter MQUANT increased from 10 to 12, i.e. MSTEP=MQUANT**2
*        is increased from 100 to 144.
*     Version 2.98, September 7, 2000
*     -- new preparation routine USER_DMS_LIMIT
*     Version 2.99, March 6, 2001
*     -- parameter MANAL increased from 30 to 50
*     -- parameter MEXP increased from 8 to 12
*     Version 3.02, July 21, 2002:
*     -- parameter MNAM increased from MCONT+30 to MCONT+50
*     Version 3.05, December 12, 2002:
*     -- adapt to Linux
*     -- parameter MNAM increased from MCONT+50 to MCONT+60
*     Version 3.10, January 21, 2003:
*     -- parameter MMETH increased from 30 to 50
*     -- parameter MPAR	increased from 40 to 50
*     Version 3.31, December 9, 2007:
*     -- adapt Scientific Linux SLC 4
*     -- parameter MNAM increased from MCONT+60 to MCONT+120
*
************************************************************************
*
      INTEGER NUNIT
      PARAMETER(NUNIT=4)       ! number of logical units for program output
      INTEGER LUNITS(NUNIT)
      INTEGER LUNLOG           ! logical unit for "log file" output  
      EQUIVALENCE(LUNLOG,LUNITS(1))
      INTEGER LUNERR           ! logical unit for warning and error messages
      EQUIVALENCE(LUNERR,LUNITS(2))
      INTEGER LUNOUT           ! logical unit for output cards
      EQUIVALENCE(LUNOUT,LUNITS(3))
      INTEGER LUNKUM           ! logical unit for kumac output
      EQUIVALENCE(LUNKUM,LUNITS(4))
      COMMON/COUNIT/LUNITS
*
*     Parameters
*     ----------
*
      INTEGER MWORDS
      PARAMETER(MWORDS=600)      ! max. number of words on a logical line
*DR (warning MQUANT is copied into MQUAN in master.inc)
      INTEGER MQUANT
      PARAMETER(MQUANT=50)     ! max. number of quantities to be averaged
      INTEGER MANAL
      PARAMETER(MANAL=200)      ! max. number of analyses
      INTEGER MEXP
      PARAMETER(MEXP=100)        ! max. number of experiments
      INTEGER MMETH
      PARAMETER(MMETH=100)      ! max. number of methods
      INTEGER MSTAT
      PARAMETER(MSTAT=100)       ! max. number of statuses
      INTEGER MSTEP
      PARAMETER(MSTEP=MQUANT**2)! max. number of steps 
*                                 (don't change this; increase MQUANT if 
*                                  a larger value of MSTEP is needed !)
      INTEGER MPAR
      PARAMETER(MPAR=150)        ! max. number of parameters
      INTEGER MCONT,FCONT
      PARAMETER(FCONT=4)        ! first systematic contribution
      PARAMETER(MCONT=MPAR+100) ! WARNING: MCONT must be at least equal to FCONT !!
*                               ! max. number of stat. and syst. contributions
      INTEGER MDISP
      PARAMETER(MDISP=3)        ! max. number of stat. and syst. contributions
*                               !  on a DISLAY logical line
      CHARACTER*(*) SUPERSEDED
      PARAMETER(SUPERSEDED='SUPERSEDED')
      INTEGER MNAM    
      PARAMETER(MNAM=MCONT+240) ! maximum number of names
      INTEGER MLINK
      PARAMETER(MLINK=50)       ! maximum number of links to other analyses
      INTEGER MCALLP
      PARAMETER(MCALLP=3)       ! maximum number of preparation routines to call
      INTEGER MCALLC
*OSOS increase MCALLC to accomodate CHI2_N_SYM
      PARAMETER(MCALLC=14)      ! maximum number of combination routines to call
      INTEGER MCALLT
      PARAMETER(MCALLT=5)       ! maximum number of termination routines to call
*
      INTEGER MSW               ! number of switches
      INTEGER LADJU,LSYMA,LSYMC,LSYST,LCOSY,
     &        LCOST,LINTE,LEXTR,LQUST,LLUMP,LSYRE,LSYAB,LSTSQ,LSTAB,
     &        LMXCR ! switches
      PARAMETER(LADJU=1,LSYMA=2,LSYMC=3,LSYST=4,LCOSY=5,
     &          LCOST=6,LINTE=7,LEXTR=8,LQUST=9,LLUMP=10,
     &          LSYRE=11,LSYAB=12,LSTSQ=13,LSTAB=14,LMXCR=15,MSW=15)
*
*     General information and names
*     -----------------------------
*
      INTEGER NEXP              ! number of experiments
      CHARACTER*20 CHEXP(MEXP)  ! names of the experiments
      INTEGER NMETH             ! number of methods
      CHARACTER*20 CHMETH(MMETH)! names of the methods
      INTEGER NSTAT             ! number of statuses
      CHARACTER*20 CHSTAT(MSTAT)! names of statuses
      INTEGER NNAM              ! total number of parameter and variable names 
      CHARACTER*20 CHNAM(MNAM)  ! names of parameters and variables
*
      COMMON/COGENE/NEXP,CHEXP,NMETH,CHMETH,NSTAT,CHSTAT,NNAM,CHNAM
*
      CHARACTER*16 CHCALP(MCALLP)! names of preparation routines to call
      CHARACTER*16 CHCALC(MCALLC)! names of combination routines to call
      CHARACTER*16 CHCALT(MCALLT)! names of termination routines to call
      DATA CHCALP/ ! List of preparation routines:
     &            'ADJU_A_FOR_FBS','USER_DMS_SENS','NEW_PREPARATION'/
*OSOS add CHI2_N_SYM
      DATA CHCALC/ ! List of combination routines:
     &            'BLUE','BLUE_ASYM','CHI2_SYM','CHI2_N_SYM',
     &            'CHI2_SYM_MIN',
     &            'CHI2_SYM_CIRC','CHI2_SYM_CIRC_G',
     &            'CHI2_SYM_INDEP','CHI2_ASYM_MIN',
     &            'COMBY_LIFETIME','DUMP_MASTER_INC',
     &            'SINGLE','CHI2_SYM_RVK','NEW_COMBINATION'/
      DATA CHCALT/ ! List of termination routines:
     &            'SET_DMS_LIMIT','KUMAC_DMS_LIMIT',
     &            'AVG_AMPLITUDE','FINAL_FRACTIONS',
     &            'NEW_TERMINATION'/ 
*
*     Results from the various analyses
*     ---------------------------------
*
      INTEGER NANAL             ! number of analyses
      INTEGER KEXP(MANAL)       ! experiment number
      INTEGER KMETH(MANAL)      ! method number
      INTEGER KSTAT(MANAL)      ! status number
      CHARACTER*64 CHREF(MANAL) ! reference
      INTEGER NSTEP(MANAL)      ! number of steps
      INTEGER KSTEP(MANAL)      ! number of step variable
      REAL STEP(MSTEP,MANAL)    ! values of step variable
      INTEGER NCONT(MANAL)      ! number of stat. + syst. uncertainties
*DR      INTEGER KCONT(0:MCONT,MANAL) ! indexes in array CHNAM
*DR now KCONT(-MQUANT+1,MANAL) to KCONT(0,MANAL) is the list of quantities to be fit
      INTEGER KCONT(-MQUANT+1:MCONT,MANAL) ! indexes in array CHNAM
      REAL CONT(MSTEP,-2*MCONT:2*MCONT,MANAL) ! measurements and their uncertainties
      LOGICAL XSYMC(0:MCONT,MANAL)! is uncertainty symmetric
      INTEGER NPAR(MANAL)       ! number of parameters
      INTEGER KPAR(MPAR,MANAL)  ! parameter numbers
      REAL PAR(-1:+1,MPAR,MANAL)! values of parameters
      INTEGER PCONT(MPAR,MANAL) ! syst. uncertainties affected by the parameters
*                               ! (or index of preparation routine if negative)
      INTEGER NLINK(MANAL)      ! number of links to other analyses 
      INTEGER LINK(0:MLINK,MANAL) ! links to analyses
      REAL RHOSTA(0:MLINK,MANAL)! statistical correlations 
      INTEGER NDISP(MANAL)      ! number of special data to be displayed
      INTEGER KDISP(0:MDISP,MANAL) ! indexes in array CHNAM
      REAL DISP(MSTEP,-MDISP:MDISP,MANAL) ! measurements and their uncertainties,
*                                         !   as they should be displayed on plots
      LOGICAL XSYMD(0:MDISP,MANAL)! is uncertainty symmetric
*
      COMMON/COANAL/NANAL,KEXP,KMETH,KSTAT,
     &              NSTEP,KSTEP,STEP,
     &              NCONT,KCONT,CONT,XSYMC, 
     &              NPAR ,KPAR ,PAR ,PCONT,
     &              NLINK,LINK ,RHOSTA,
     &              CHREF,
     &              NDISP,KDISP,DISP,XSYMD
*
*     Instructions for combination
*     ----------------------------
*
      INTEGER NQUANT         ! number of quantities to be averaged
      LOGICAL XCOMB(0:MANAL) ! do we combine this analysis ?
      LOGICAL XON(MSW)       ! is switch ON ?
      INTEGER NDEST          ! number of desired steps
      INTEGER DEST(MSTEP)    ! desired steps
      INTEGER KSYN(0:MNAM)   ! index for synonyms
      INTEGER KLUMP(0:MCONT,MANAL) ! index for lumps
      INTEGER NCALLP         ! number of preparation routines to call
      INTEGER KCALLP(MCALLP) ! preparation routines to call
      INTEGER NCALLC         ! number of combination routines to call
      INTEGER KCALLC(MCALLC) ! combination routines to call
      INTEGER NCALLT         ! number of combination routines to call
      INTEGER KCALLT(MCALLT) ! termination routines to call
      INTEGER LCALLP(MCALLC) ! logical units for printout (preparation routines)
      INTEGER LCALLC(MCALLC) ! logical units for printout (combination routines)
      INTEGER LCALLT(MCALLT) ! logical units for printout (termination routines)
*OSOS add NQUANT
      COMMON/CODRIV/NQUANT,XCOMB,XON,NDEST,DEST,KSYN,KLUMP,
     &              NCALLP,NCALLC,NCALLT,KCALLP,KCALLC,KCALLT,
     &              LCALLP,LCALLC,LCALLT
*
*****************************************************************************
*
* If IANAL, ICONT and ISTEP are defined as:
*
*  IANAL = analysis number (from 1 to NANAL)
*  ICONT = systematic uncertainty number (from FCONT to NCONT(IANAL))
*  ISTEP = step number (from 1 to NSTEP(IANAL))
*  IPAR  = parameter number (from 1 to NPAR(IANAL))
*
* then:
*
*  CHEXP(KEXP(IANAL)         = name of experiment
*  CHMETH(KEXP(IANAL))       = name of method
*  CHSTAT(KSTAT(IANAL))      = name of status
*  CHNAM(KCONT(0,IANAL))     = name of measured quantity
*  CHNAM(KCONT(1,IANAL))     = name of statistical uncertainty
*  CHNAM(KCONT(2,IANAL))     = name of total systematic uncertainty
*  CHNAM(KCONT(3,IANAL))     = name of total uncertainty
*  CHNAM(KCONT(ICONT,IANAL)) = name of systematic uncertainty
*  CONT(ISTEP, 0,IANAL)      = central value of measurement
*  CONT(ISTEP,+1,IANAL)      = positive statistical uncertainty (> 0)
*  CONT(ISTEP,-1,IANAL)      = negative statistical uncertainty (< 0)
*  CONT(ISTEP,+2,IANAL)      = total positive systematic uncertainty (> 0)
*  CONT(ISTEP,-2,IANAL)      = total negative systematic uncertainty (< 0)
*  CONT(ISTEP,+3,IANAL)      = total positive uncertainty (> 0)
*  CONT(ISTEP,-3,IANAL)      = total negative uncertainty (< 0)
*  CONT(ISTEP,+ICONT,IANAL)  = systematic uncertainty for positive excursion
*  CONT(ISTEP,-ICONT,IANAL)  = systematic uncertainty for negative excursion
*  CONT(ISTEP,+MCONT+ICONT,IANAL)  = power for error scaling
*  CONT(ISTEP,-MCONT-ICONT,IANAL)  = power for error scaling
*  PAR( 0,IPAR,IANAL)        = central value of parameter
*  PAR(+1,IPAR,IANAL)        = positive excursion of parameter
*  PAR(-1,IPAR,IANAL)        = negative excursion of parameter
*  PCONT(IPAR,IANAL)         = value of ICONT associated with parameter, or 0
*
*  XCOMB(0)    =.TRUE.    --> a combination of analyses is required
*  XCOMB(IANAL)=.TRUE.    --> analysis IANAL has to be included in combination 
*  KSYN(I)=KSYN(J)        --> CHNAM(I) and CHNAM(J) have to be treated as 
*                              synonyms in the combination 
*  KLUMP(ICONT,IANAL).NE.KCONT(ICONT,IANAL)
*                         --> systematic uncertainty CHNAM(KCONT(ICONT,IANAL))
*                               belongs to lump CHNAM(KLUMP(ICONT,IANAL))
*
*****************************************************************************
*
*     Information for performing MINOS error calculation
*
      LOGICAL XMINOS
      CHARACTER*16 CHMINOS(MWORDS)
      INTEGER NMINOS
      COMMON /CMINOS/CHMINOS,NMINOS,XMINOS
* CHMINOS -- List of parameters for which MINOS errors are to be calculated;
*            if one parameter is 'ALL' calculate them for all parameters
* NMINOS  -- Number of entries in the array CHMINOS
* XMINOS  -- Flag to indicate if MINOS errors to be calculated
*****************************************************************************
      INTEGER NSPAR, MSPAR	
      PARAMETER (MSPAR = 1000)
      CHARACTER*16 CHSPAR(MSPAR)
      DOUBLE PRECISION SPAR1(MSPAR), SPAR2(MSPAR)
      COMMON /SPAR/ NSPAR,CHSPAR,SPAR1,SPAR2
*****************************************************************************
