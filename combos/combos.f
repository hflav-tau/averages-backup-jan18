************************************************************************
*
*     COMBOS:  program to combine LEP B oscillation results
*     =======             COMbine     B OScillation
*                         COM         B OS
*                          COM       B OS
*                           COM     B OS
*                            COM  B OS
*                             COMB OS
*                             COMBOS
*
*     Olivier Schneider, CERN/PPE-ALE
*
*********
*
*     Version 1.00, December  6, 1996:
*     -- original release
*     Version 1.10, December 10, 1996:
*     -- now able to handle statistical correlations between analyses
*     -- now able to handle analyses that supersede other analyses
*     -- added 4 keywords: SUPERSEDE, STAT_CORR_WITH, NOSYSTEMATICS, SYMMETRIZE
*     -- changed characteristics of NOCORRELATIONS keyword
*     Version 1.20, December 13, 1996:
*     -- read name of the input file(s) from the command line
*     -- allow combinations involving a single analysis
*     Version 1.21, January 17, 1997:
*     -- bug fixed in parameter inheritance (combined analysis could inherit
*        parameters from an inappropriate analysis)
*     Version 1.30, January 27, 1997:
*     -- added new keyword: SYNONYMS
*     -- handle synonyms (array KSYN)
*     Version 1.40, February 6, 1997:
*     -- remove dependence on ICNTH in routine CHECK_CHAR
*     -- added new keyword: LUMP
*     -- handle lumps (array KLUMP, routines STORE_LUMP and STORE_COMB_GET)
*     Version 1.50, February 12, 1997:
*     -- bug fixed in the decoding of the DATA logical lines (new entry
*        STORE_DATA_RESET)
*     -- bug fixed in the decoding of STEP and MEAS logical lines (these logical
*        lines cannot contain a common name anymore)
*     -- handle new format of the STEP logical line (numerical values allowed)
*     Version 2.00, February 20, 1997:
*     -- removed 4 obsolete keywords: NOADJUSTMENTS, NOCORRELATIONS,
*                                     NOSYSTEMATICS, SYMMETRIZE
*     -- added new keyword: SWITCH
*     -- handle switches (array XON); new switches are LUMPS, QUOTE_STAT, and 
*         INTERPOLATION; old switches (i.e. new switches replacing the 4 
*         4 obsolete keywords) are ADJUSTMENTS, SYMMETRIZATION, SYSTEMATICS,
*         SYST_CORR, and STAT_CORR
*     -- handle new format of the MEASUREMENTS logical line;
*         format is now the same for individual analyses or combined analyses
*     -- handle new format of the CALL logical line; both combination and 
*         termination routines can be specified on the CALL logical line; 
*         the logical unit for the output of these routines is now also 
*         specified directly on the CALL logical line
*     -- bug fixed in decoding of CALL logical line: each name after the CALL
*         keyword is now allowed only once, or not allowed if it conflicts 
*         with the default value of a missing name
*     -- the specification of the total systematic uncertainties is not
*         required anymore, since these are not used in the combinations
*     -- the computation of the total systematic uncertainties and the total 
*         uncertainties for an individual analysis can be requested on the 
*         MEASUREMENTS logical line by explicitely specifying a name for these
*         quantitites
*     -- new routine CHECK_STEPS to perform all the nevessary checks on the
*         step values (can handle both an individual analysis and a combined 
*         analysis); if no step values are provided for a combined analyses
*         (on the STEP logical line), all the steps values used in the 
*         individual analyses included in the combination are used; 
*         step values in a combined analysis that would require extrapolation
*         in at least one of the individual analyses included in the
*         combination are rejected
*     -- suppressed warning messge in case that only one single analysis is
*         included in a combination
*     -- increased dimension of 5th argument in call to MASTER to receive
*         the statistical and systematic uncertainties on the combined results
*         (in addition to the total uncertaintities)
*     Version 2.10, February 27, 1997:
*     -- new EXTRAPOLATION switch (parameter LEXTR)
*     -- routine ANALYSIS_NAME has now an extra argument VERBOSE
*     -- step rejection procedure in CHECK_STEPS modified: only step values in 
*        a combined analysis that would require extrapolation in all the
*        individual analyses included in the combination are rejected
*     Version 2.12, March 8, 1997:
*     -- implement missing overwrite protection in routine STORE_STEP
*        (program could write more than MSTEP=100 words in array DEST)
*     Version 2.30, May 12, 1997:
*     -- print name of operating system (environmental variable OS)
*        in output files in addition to version number, date and time
*     -- SYMMETRIZATION switch renamed SYMM_COMB (parameter LSYMM renamed LSYMC)
*     -- new SYMM_ADJU switch (parameter LSYMA)
*     -- add 4th logical unit for kumac output
*     Version 2.40, October 24, 1997:
*     -- add logic and routine STORE_MINO for minos error option (HCJS)
*     Version 2.50, March 9, 1998:
*     -- move XMINOS=.FALSE. initialization from routine COMBOS to 
*        routine STORE_BEGI
*     -- handle new DISPLAY logical line
*     -- fix bug in COMB_ANAL (inheritance of step name)
*     Version 2.52, July 12, 1999:
*     -- introduce concept of preparation routine
*     -- parameter name in an analysis to be combined must be the name 
*        of a systematic contribution or of a preparation routine 
*     Version 3.00, June 1, 1999 (David Rousseau):
*     -- modifications to be able to average two correlated quantities
*     Version 3.04, December 5, 2002 (Fabrizio Parodi):
*     -- Linux version
*        - Main program in C (combos_wrapper.c)
*        - CHARACTER*(*) in concatenation: solved using constant CHARACTER*256
*        - EQUIVALENCE   eliminated by simple copy
*     Version 3.10, January 21, 2003:
*     -- rename main program, now called cmain.c (E. Barberio)
*     Version 3.20, February 5, 2003:
*     -- first release for HFAG (to be put in CVS)   
*     Version 3.22, February 10, 2003:
*     -- call CHI2_SYM instead of CHI2_SYM_CIRC_G or CHI2_SYM_CIRC to get 
*        "stat only" result
*     Version 3.25, May 20, 2005:
*     -- correct small bug in the computation of SYTOT in WRITE_KUMAC
*     Version 3.30, July 3, 2005:
*     -- put in modifications by Bob Kowalewski (routine CHI2_SYM_RVK 
*        and related input format and switches)
*     Version 3.31, December 9, 2007:
*     -- adapt to Scientific Linux SLC 4
*     -- parameter MNAM increased from MCONT+60 to MCONT+120
*     Version 3.32, January 17, 2010
*     -- patches (from Alberto Lusiani) to compile with gfortran
*        on Scientific Linux 5 
*     -- new printout format for the results of the combination routines
*        allowing for NQUAN=2 (in new subroutine PRINT_QUAN)
*     Version 3.33, March 20, 2010  
*     -- improved printout in DUMP_MASTER_INC and in CHI2_N_SYM
*        (provided by Swagato)
*
************************************************************************
*
      SUBROUTINE COMBOS(LEN,INPUT)
*     ========================
*
      IMPLICIT NONE
      INCLUDE 'combos.inc'
*
*     Argument
*
      CHARACTER*256 INPUT  
      INTEGER       LEN
*
*     Local variables
*
      INTEGER NWORDS,IERR
      CHARACTER*20 KEY,WORDS(MWORDS)
*
*     Open input file
*
      IF(INPUT(1:LEN).NE.' ') THEN
        OPEN(UNIT=5,FILE=INPUT(1:LEN),STATUS='OLD',
     &       FORM='FORMATTED',ACCESS='SEQUENTIAL',IOSTAT=IERR)
        IF(IERR.NE.0) THEN 
          WRITE(*,*)'Cannot open file "'//INPUT(1:LEN)//'"'  
          RETURN
        ENDIF
        CALL SET_INPUT_NAME(LEN,INPUT)
      ENDIF
*
*     Initialize common blocks
*
      NEXP=0
      NMETH=0
      NSTAT=0
      NNAM=0
      NANAL=0
*
*     Set default values for the output logical units 
*
      CALL STORE_OUTP_DEFAULTS
*
*     Get next keyword from input files
*
    1 CALL GET_NEXT_KEY(KEY,MWORDS,NWORDS,WORDS)
*
*     Take action depending on keyword
*
      IF(KEY.EQ.' ') THEN 
        CALL COMBOS_ERROR_SUMMARY
        RETURN
      ELSE IF(KEY.EQ.'OUTP') THEN 
        CALL STORE_OUTP(NWORDS,WORDS)
      ELSE IF(KEY.EQ.'BEGI') THEN 
        CALL STORE_BEGI(NWORDS,WORDS)
      ELSE IF(KEY.EQ.'PARA') THEN 
        CALL STORE_PARA(NWORDS,WORDS)
      ELSE IF(KEY.EQ.'STEP') THEN 
        CALL STORE_STEP(NWORDS,WORDS)
      ELSE IF(KEY.EQ.'MEAS') THEN 
        CALL STORE_MEAS(NWORDS,WORDS)
      ELSE IF(KEY.EQ.'DISP') THEN 
        CALL STORE_DISP(NWORDS,WORDS)
      ELSE IF(KEY.EQ.'DATA') THEN 
        CALL STORE_DATA(NWORDS,WORDS)
      ELSE IF(KEY.EQ.'SUPE') THEN 
        CALL STORE_LINK('SUPERSEDE',NWORDS,WORDS)
      ELSE IF(KEY.EQ.'STAT') THEN 
        CALL STORE_LINK('STAT_CORR_WITH',NWORDS,WORDS)
      ELSE IF(KEY.EQ.'COMB') THEN 
        CALL STORE_COMB(NWORDS,WORDS)
      ELSE IF(KEY.EQ.'LUMP') THEN 
        CALL STORE_LUMP(NWORDS,WORDS)
      ELSE IF(KEY.EQ.'CALL') THEN 
        CALL STORE_CALL(NWORDS,WORDS)
      ELSE IF(KEY.EQ.'SYNO') THEN 
        CALL STORE_SYNO(NWORDS,WORDS)
      ELSE IF(KEY.EQ.'SWIT') THEN 
        CALL STORE_SWIT(NWORDS,WORDS)
      ELSE IF(KEY.EQ.'MINO') THEN 
        CALL STORE_MINO(NWORDS,WORDS)
      ELSE IF(KEY.EQ.'END') THEN 
        IF(XCOMB(0)) THEN
          CALL COMB_ANAL(NWORDS,WORDS,IERR)
        ELSE
          CALL CHECK_ANAL(NWORDS,WORDS,IERR)
        ENDIF
        IF(LUNOUT.GT.0.AND.IERR.EQ.0) CALL WRITE_CARDS(NANAL)
        IF(XCOMB(0).OR.IERR.NE.0) NANAL=NANAL-1 ! discard this analysis
      ELSE 
        PRINT * ,'Unrecognized key word: ',KEY
        STOP 99
      ENDIF
      GOTO 1
      END
*
************************************************************************
*
      SUBROUTINE COMBOS_ERROR(ISEVER,TEXT,ACTION)
*     ===========================================
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER ISEVER
      CHARACTER*(*) TEXT,ACTION
*
*     Extra arguments of entry
*
      CHARACTER*(*) RFILE
      INTEGER RLINE
*
*     Externals
*
      INTEGER LENOCC
*
*     Local variables
*
      INTEGER WDEPTH,WUNIT,WLINE,WWORD,IS
      CHARACTER*1 Q
      CHARACTER*80 WFILE,NAME
      CHARACTER*7 SEVER(-1:1)
      DATA SEVER/'FATAL','WARNING','ERROR'/
      INTEGER COUNT(-1:1)
      SAVE    COUNT
      DATA    COUNT/3*0/
*
      INCLUDE 'combos.inc'
*
      IS=MIN0(MAX0(ISEVER,-1),1)
      COUNT(IS)=COUNT(IS)+1
      IF(LUNERR.LE.0) RETURN ! error messages suppressed
*
*     Where are we in the input files ?
*
      CALL GET_WHERE(WDEPTH,WUNIT,WFILE,WLINE,WWORD)
*
*     Get name of current analysis
*
      CALL ANALYSIS_NAME(.TRUE.,NANAL,NAME)
*
      Q=' '
      IF(NAME.EQ.' ') Q=':'
      WRITE(LUNERR,8881) SEVER(IS),WWORD,WLINE,WFILE(:LENOCC(WFILE)),Q
      IF(NAME  .NE.' ') WRITE(LUNERR,8882) NAME  (:LENOCC(NAME  ))
    8 IF(TEXT  .NE.' ') WRITE(LUNERR,8883) TEXT  (:LENOCC(TEXT  ))
      IF(ACTION.NE.' ') WRITE(LUNERR,8884) ACTION(:LENOCC(ACTION))
 8881 FORMAT(1X,A7,' at (or before) word',I4,' of line',I4,
     &       ' in file ',2A)
 8882 FORMAT(T10,'concerning ',A,':')
 8883 FORMAT(T10,'"',A,'"')
 8884 FORMAT(T10,'---> ',A)
      IF(ISEVER.LT.0) STOP 1
      RETURN
*
      ENTRY READER_ERROR(RFILE,RLINE,ISEVER,TEXT,ACTION)
*     ==================================================
*   
      IS=MIN0(MAX0(ISEVER,-1),1)
      COUNT(IS)=COUNT(IS)+1
      IF(LUNERR.LE.0) RETURN ! error messages suppressed
      WRITE(LUNERR,8885) SEVER(IS),RLINE,RFILE(:LENOCC(RFILE))
 8885 FORMAT(1X,A7,' at line',I4,' in file ',A,':')
      GOTO 8
*
      ENTRY COMBOS_ERROR_SUMMARY
*     ==========================
*   
      WRITE(LUNLOG,'(1X)')
      IF(LUNERR.GT.0.AND.LUNERR.NE.LUNLOG) WRITE(LUNERR,'(1X)')
      DO IS=-1,1
        IF(COUNT(IS).GT.0) THEN 
          WRITE(LUNLOG,8889) SEVER(IS),COUNT(IS)
          IF(LUNERR.GT.0.AND.LUNERR.NE.LUNLOG) 
     &    WRITE(LUNERR,8889) SEVER(IS),COUNT(IS)
 8889     FORMAT(1X,'Total number of ',A,
     &           ' messages generated by COMBOS:',I7)
        ENDIF
      ENDDO
      END
*
************************************************************************
*
      SUBROUTINE ANALYSIS_NAME(VERBOSE,IANAL,NAME)
*     ============================================
*
      IMPLICIT NONE
*
*     Arguments
*
      LOGICAL VERBOSE
      INTEGER IANAL
      CHARACTER*(*) NAME
*
*     Externals
*
      INTEGER LENOCC
*
      INCLUDE 'combos.inc'
*
      IF(IANAL.GT.0.AND.IANAL.LE.NANAL.AND.IANAL.LE.MANAL) THEN
        NAME=CHEXP (KEXP (IANAL))(:LENOCC(CHEXP (KEXP (IANAL))))//' '//
     &       CHMETH(KMETH(IANAL))(:LENOCC(CHMETH(KMETH(IANAL))))//' '//
     &       CHSTAT(KSTAT(IANAL))(:LENOCC(CHSTAT(KSTAT(IANAL))))
        IF(VERBOSE) THEN
          NAME='analysis '//NAME
          IF(IANAL.EQ.NANAL.AND.XCOMB(0)) NAME='combined '//NAME
        ENDIF
      ELSE
        NAME=' '
      ENDIF
      END
*
************************************************************************
*
      SUBROUTINE STORE_OUTP(NWORDS,WORDS)
*     ===================================
*
*     Routine called for the "OUTPUT" keyword: 
*     logical units for program output are read in and stored
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORDS(*)
*
*     Externals
*
      INTEGER LENOCC
*
*     Program version 
*
      REAL VERSION_NUMBER
      PARAMETER(VERSION_NUMBER=3.33)  ! program version number
      CHARACTER*(*) VERSION_DATE
      PARAMETER(VERSION_DATE='Mar 20, 2010') ! date of version
*
*     Local variable
*
      CHARACTER*8 CHDAY,CHTIME
      CHARACTER*32 CHOS
      INTEGER I,NVAL,LUN
      INTEGER OLDLOG,OLDERR,OLDOUT,OLDKUM
      REAL VAL
*
      INCLUDE 'combos.inc'
*
*     Save old logical units
*
      OLDLOG=LUNLOG
      OLDERR=LUNERR
      OLDOUT=LUNOUT
      OLDKUM=LUNKUM
*
*     Get new logical units
*
      I=0
    1 CALL GET_NEXT_VALUE(1,NVAL,VAL)
      IF(NVAL.GT.0) THEN 
        LUN=NINT(VAL)
        I=I+1
        IF(LUN.EQ.5.OR.LUN.GT.89) THEN
          CALL COMBOS_ERROR(1,'Invalid logical unit',
     &                        'logical unit not changed')
        ELSE IF(LUN.NE.0.AND.I.LE.NUNIT) THEN 
          LUNITS(I)=LUN
        ENDIF
        GOTO 1
      ENDIF
      IF(LUNLOG.LE.0) LUNLOG=6
      GOTO 2
*
      ENTRY STORE_OUTP_DEFAULTS
*     =========================
*
      OLDLOG=-1
      OLDERR=-1
      OLDOUT=-1
      OLDKUM=-1
      LUNLOG=6 ! default logical unit for log file output
      LUNERR=6 ! default logical unit for warning and error messages
      LUNOUT=0 ! default logical unit for output cards
      LUNKUM=0 ! default logical unit for kumac output
*
*     Print version number in all output files
*
    2 CALL DATIMH(CHDAY,CHTIME) 
      CALL GETENVF('OS',CHOS)
      IF(LUNLOG.NE.OLDLOG)
     &   WRITE(LUNLOG,1000) VERSION_NUMBER,VERSION_DATE,
     &                      CHOS(:LENOCC(CHOS)),CHDAY,CHTIME
      IF(LUNERR.NE.OLDERR.AND.LUNERR.NE.LUNLOG.AND.LUNERR.GT.0)
     &   WRITE(LUNERR,1001) VERSION_NUMBER,CHOS(:LENOCC(CHOS)),
     &                      CHDAY,CHTIME
      IF(LUNOUT.NE.OLDOUT.AND.LUNOUT.GT.0)
     &   WRITE(LUNOUT,1002) VERSION_NUMBER,CHOS(:LENOCC(CHOS)),
     &                      CHDAY,CHTIME
      IF(LUNKUM.NE.OLDKUM.AND.LUNKUM.GT.0) THEN
         WRITE(LUNKUM,1002) VERSION_NUMBER,CHOS(:LENOCC(CHOS)),
     &                      CHDAY,CHTIME
         CALL WRITE_KUMAC_INIT
      ENDIF
 1000 FORMAT(1X,'COMBOS version ',F5.2,' (',A,') run on ',A,
     &                                             ' on ',A,' at ',A,/)
* 8000 FORMAT(T18,'LEP B oscillation results combination program',
*     &     /,T18,'     COMBOS version ',F5.2,'  (',A,')',
*     &     /,T18,'         run on ',A,' at ',A,/)
* 9000 FORMAT(1X,'COMBOS version',F5.2,' (',A,
*     &  '): LEP B oscillation results combination program',/)
 1001 FORMAT(1X,'Messages generated by COMBOS version ',
     &  F5.2,' (',A,') on ',A,' at ',A,':',/)
 1002 FORMAT('* File written by COMBOS version ',F5.2,' (',A,')',
     &          ' on ',A,' at ',A,/
     &       '**************************************',8('****'))
      END
*
************************************************************************
*
      SUBROUTINE STORE_BEGI(NWORDS,WORDS)
*     ===================================
*
*     Routine called for the "BEGIN" keyword: 
*     start of a new analysis definition, the experiment name, method, 
*     status and reference (which caracterize the analysis)
*     are read in and stored
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORDS(*)
*
*     Externals
*  
      INTEGER LENOCC
*
*     Local variable
*
      INTEGER I,IANAL,ICONT,IWORD,DUMMY,NEWANAL,IPOS
      CHARACTER*20 STRING
      LOGICAL XVALUE
      REAL VALUE
*
      INCLUDE 'combos.inc'
*
*     Reset synonyms and lumps
*
      DO I=0,NNAM
        KSYN(I)=I
      ENDDO
      DO IANAL=1,NANAL
        DO ICONT=0,NCONT(IANAL)
          KLUMP(ICONT,IANAL)=KCONT(ICONT,IANAL)
        ENDDO
      ENDDO
*
*     Reset STORE_DATA
*
      CALL STORE_DATA_RESET
*
*     Create new analysis
*
      NANAL=NANAL+1
      IF(NANAL.GT.MANAL) CALL COMBOS_ERROR(-1,
     & 'Too many analyses','please increase parameter MANAL')
      NQUANT=0
      NSTEP(NANAL)=0
      KSTEP(NANAL)=0
      NPAR(NANAL)=0
      NCONT(NANAL)=FCONT-1
      NDISP(NANAL)=-1
*DR     CALL UZERO(KCONT(0,NANAL),1,MCONT+1)
      CALL UZERO(KCONT(1-MQUANT,NANAL),1,MCONT+MQUANT)
      DO I=0,MANAL
        XCOMB(I)=.FALSE.
      ENDDO
      XMINOS = .FALSE.
      CALL RESET_SWIT
      NDEST=0
      NCALLP=0
      NCALLC=0
      NCALLT=0
*
*     Set NANAL temporarily to zero
*     (to prevent COMBOS_ERROR from printing the "analysis name")
*
      NEWANAL=NANAL
      NANAL=0
*
*     Experiment name
*
      IWORD=1
      I=MIN0(MAX0(NWORDS-IWORD+1,0),1)
      CALL CHOOSE_NAME(I,WORDS(IWORD),'experiment','E',NEXP+1)
      CALL CHECK_CHAR(WORDS(IWORD),MEXP,NEXP,CHEXP,'experiment','MEXP',
     &               .TRUE.,KEXP(NEWANAL),DUMMY)
*
*     Method name
*
      IWORD=2
      IF(NWORDS.LT.IWORD) WORDS(IWORD)=' '
      CALL CHECK_CHAR(WORDS(IWORD),MMETH,NMETH,CHMETH,'method','MMETH',
     &                .TRUE.,KMETH(NEWANAL),DUMMY)
*
*     Status
*
      IWORD=3
      IF(NWORDS.LT.IWORD) WORDS(IWORD)=' '
      CALL CHECK_CHAR(WORDS(IWORD),MSTAT,NSTAT,CHSTAT,'status','MSTAT',
     &                .TRUE.,KSTAT(NEWANAL),DUMMY)
*
*     Set NANAL back to its normal value
*
      NANAL=NEWANAL
*
*     Reference
*
      CHREF(NANAL)=' '
    9 CALL GET_NEXT_WORD(IPOS,STRING,XVALUE,VALUE)
      IF(XVALUE.OR.IPOS.GT.1) THEN
        CHREF(NANAL)(LENOCC(CHREF(NANAL))+2:)=STRING
        GOTO 9
      ENDIF
      CHREF(NANAL)=CHREF(NANAL)(2:)
      CALL SET_SAME_WORD
      DO IWORD=NWORDS,4,-1
        CHREF(NANAL)=
     &          WORDS(IWORD)(:LENOCC(WORDS(IWORD)))//' '//CHREF(NANAL)
      ENDDO
*
*     Link to itself
*
      NLINK(NANAL)=0
      LINK(0,NANAL)=1000000*KEXP(NANAL)+1000*KMETH(NANAL)+KSTAT(NANAL)
      RHOSTA(0,NANAL)=1. ! statistical correlation
      END
*
************************************************************************
*
      SUBROUTINE STORE_PARA(NWORDS,WORD)
*     ==================================
*
*     Routine called for the "PARAMETER" keyword: 
*     the parameter name, value, positive and negetive excusrions are 
*     read in and stored
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORD
*
*     Externals
*
      INTEGER LENOCC
*
*     Local variable
*
      INTEGER NVAL,I,DUMMY,ISAME
      REAL VAL(3)
*
      INCLUDE 'combos.inc'
*
*     Create new parameter
*
      IF(NPAR(NANAL).EQ.MPAR) CALL COMBOS_ERROR(-1,
     & 'Too many parameters','please increase parameter MPAR')
      NPAR(NANAL)=NPAR(NANAL)+1
      CALL CHOOSE_NAME(NWORDS,WORD,'parameter','P',NPAR(NANAL))
      CALL CHECK_NAME(WORD,.TRUE.,KPAR(NPAR(NANAL),NANAL),DUMMY)
      CALL GET_NEXT_VALUE(3,NVAL,VAL)
      IF(NVAL.EQ.0) CALL COMBOS_ERROR(0,'No value for parameter '//
     & CHNAM(KPAR(NPAR(NANAL),NANAL))(:LENOCC(
     & CHNAM(KPAR(NPAR(NANAL),NANAL)))),'Value assumed to be 0')
      PAR( 0,NPAR(NANAL),NANAL)=VAL(1) ! central value
      PAR(+1,NPAR(NANAL),NANAL)=AMAX1(VAL(2),VAL(3)) ! positive excursion
      PAR(+1,NPAR(NANAL),NANAL)=AMAX1(PAR(+1,NPAR(NANAL),NANAL),0.)
      PAR(-1,NPAR(NANAL),NANAL)=AMIN1(VAL(2),VAL(3)) ! negative excursion
      PAR(-1,NPAR(NANAL),NANAL)=AMIN1(PAR(-1,NPAR(NANAL),NANAL),0.)
*
*     Do not allow two parameters with same name
*
      ISAME=0
      DO I=1,NPAR(NANAL)-1
        IF(KPAR(I,NANAL).EQ.KPAR(NPAR(NANAL),NANAL)) THEN
          CALL COMBOS_ERROR(1,'Another parameter has the same name',
     &                        'this parameter ignored')
          ISAME=ISAME+1
        ENDIF
      ENDDO
      IF(ISAME.NE.0) NPAR(NANAL)=NPAR(NANAL)-1
      END
*
************************************************************************
*
      SUBROUTINE STORE_STEP(NWORDS,WORD)
*     ==================================
*
*     Routine called for the "STEP" keyword: 
*     the name of the step variable is read in and stored
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORD
      CHARACTER*256 XWORD
*
*     Externals
*
      INTEGER LENOCC
*
*     Local variable
*
      INTEGER DUMMY,I
      INTEGER NVAL
*
      INCLUDE 'combos.inc'
      INTEGER MVAL
      PARAMETER(MVAL=MSTEP+1)
      REAL VAL(MVAL)
      CHARACTER*4 CH4
*
*     Get numerical values
*
      NVAL=0
    1 NVAL=NVAL+1
      CALL GET_NEXT_VALUE(1,I,VAL(MIN0(NVAL,MVAL)))
      IF(I.GT.0) GOTO 1
      NVAL=NVAL-1
      IF(NVAL.GT.MSTEP) THEN
        NVAL=MSTEP
        WRITE(CH4,'(I4)') NVAL
        CALL COMBOS_ERROR(0,'Too many values specified on STEP keyword',
     &    'only first'//CH4//' values used; extra values ignored')
      ENDIF
*
*     Check that this is the first STEP logical line
*
      IF(KSTEP(NANAL).GT.0) THEN 
        CALL COMBOS_ERROR(0,'More than one STEP keyword',
     &                      'this STEP logical line ignored')
        RETURN
      ENDIF
*
*     Check step name
*
      CALL CHOOSE_NAME(NWORDS,WORD,'step','STEP',0)
      CALL CHECK_NAME(WORD,.TRUE.,KSTEP(NANAL),DUMMY)
      DO I=0,NCONT(NANAL)
        IF(KCONT(I,NANAL).EQ.KSTEP(NANAL)) THEN 
          XWORD = WORD
          CALL COMBOS_ERROR(1,'Name '//XWORD(:LENOCC(WORD))//
     &    ' already appeared on a DATA or MEASUREMENT logical line',
     &    'this STEP logical line ignored')
          KSTEP(NANAL)=0
          RETURN
        ENDIF
      ENDDO
*
*     Load step values if provided
*
      CALL UZERO(STEP(1,NANAL),1,MSTEP)
      NSTEP(NANAL)=NVAL
      DO I=1,NSTEP(NANAL)
        STEP(I,NANAL)=VAL(I)
      ENDDO
      IF(NVAL.GT.0) CALL STORE_DATA_STEP_LOADED
      END
*
************************************************************************
*
      SUBROUTINE STORE_MEAS(NWORDS,WORDS)
*     ===================================
*
*     Routine called for the "MEASUREMENT" keyword: 
*     the name of the measured quantity, as well as the names of the 
*     statistical, total systematic, and total uncertainties 
*     are read in and stored
*
*OSOS start
*     - count number of MEASUREMENT keywords (variable NQUANT)
*     - only the following names from the MEASUREMENT logical lines are saved:
*          - the name of the measurement of each logical line
*          - the names of the statistical, total systematic, and total 
*            uncertainties of the first logical line only (the corresponding 
*            names from the next logical lines of the same analysis 
*            being lost ...
*OSOS end
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORDS(*)
*
*     Externals
*
      INTEGER LENOCC
*
*     Local variable
*
*OSOS start
      INTEGER IWORD,I,DUMMY,ICONT,JCONT,IKCONT
*OSOS end
      CHARACTER*16 CH16
*
      INCLUDE 'combos.inc'
      CHARACTER*11 DNAM(FCONT,0:1)
      DATA DNAM/'measurement','stat. error','syst. error','total error',
     &          'MEAS','STAT','SYST','TOTAL_ERROR'/ ! default names 
*
*     Check that there is a measurement name and statistical error 
*
      IF(NWORDS.LE.0) THEN 
        CALL COMBOS_ERROR(1,'No measurement name specified',
     &                      'this MEASUREMENT logical line ignored')
        RETURN
      ENDIF
      IF(NWORDS.EQ.1.AND..NOT.XCOMB(0)) THEN
        CALL COMBOS_ERROR(1,'No name specified for statistical error',
     &                      'this MEASUREMENT logical line ignored')
        RETURN
      ENDIF
*
*OSOS start
      IF(NQUANT.GE.MQUANT) CALL COMBOS_ERROR(0,
     & 'Too many MEASUREMENT keywords',
     & 'this logical line will overwrite the previous definition')
*OSOS        IF(KCONT(0,NANAL).GT.0) CALL COMBOS_ERROR(0,
*OSOS     & 'More than one MEASUREMENT keyword',
*OSOS     & 'this logical line will overwrite the previous definition')
*OSOS end
      IF(NCONT(NANAL).GT.FCONT-1) CALL COMBOS_ERROR(1,
     & 'MEASUREMENT keyword after DATA keyword(s)',
     & 'previous data (except step values) will be ignored')
      NCONT(NANAL)=FCONT-1
      DO ICONT=0,NCONT(NANAL)
        CALL UZERO(CONT(1,+ICONT,NANAL),1,MSTEP)
        CALL UZERO(CONT(1,-ICONT,NANAL),1,MSTEP)
crvk
        CALL UZERO(CONT(1,+MCONT+ICONT,NANAL),1,MSTEP)
        CALL UZERO(CONT(1,-MCONT-ICONT,NANAL),1,MSTEP)
        IWORD=ICONT+1
        IF(IWORD.LE.NWORDS) THEN ! name specified
          CH16=WORDS(IWORD)
          I=1
        ELSE IF(XCOMB(0)) THEN ! name not specified, combined analysis
          I=1
          CH16=DNAM(IWORD,1) ! to suppress message from CHOOSE_NAME
        ELSE IF(IWORD.GT.2) THEN ! syst or tot error name not specified,
*                                  individual analysis
          I=-1 ! to not call CHOOSE_NAME
          CH16=' '
          KCONT(ICONT,NANAL)=0
        ELSE ! stat error name not specified, individual analysis
          I=0
          CH16=' '
        ENDIF
        IF(I.GE.0) THEN
          CALL CHOOSE_NAME(I,CH16,DNAM(IWORD,0),DNAM(IWORD,1),0)
*OSOS start
          IF(ICONT.EQ.0) THEN
            NQUANT=MIN(NQUANT+1,MQUANT)
            IKCONT=1-NQUANT
          ELSE
            IKCONT=ICONT
          ENDIF
*OSOS     PRINT * ,'NANAL,IKCONT = ',NANAL,IKCONT
          IF(NQUANT.EQ.1.OR.ICONT.EQ.0) THEN
            CALL CHECK_NAME(CH16,.TRUE.,KCONT(IKCONT,NANAL),DUMMY)
*OSOS       PRINT * ,'KCONT(IKCONT,NANAL) = ',KCONT(IKCONT,NANAL),
*OSOS& CHNAM(KCONT(IKCONT,NANAL))
            IF(KCONT(IKCONT,NANAL).EQ.KSTEP(NANAL)) THEN
              CALL COMBOS_ERROR(1,'Name '//CH16(:LENOCC(CH16))//
     &         ' already appeared on a STEP logical line',
     &         'this MEASUREMENT logical line ignored')
              IF(NQUANT.LE.1) THEN
                CALL UZERO(KCONT(0,NANAL),1,NCONT(NANAL)+1)
              ELSE
                KCONT(1-NQUANT,NANAL)=0
              ENDIF
              NQUANT=NQUANT-1
              RETURN
            ENDIF
*OSOS end
          ENDIF
*OSOS     PRINT * ,' ' 
*OSOS     PRINT * ,'NANAL = ',NANAL
          DO JCONT=0,ICONT-1
*OSOS       PRINT * ,'JCONT,ICONT = ',JCONT,ICONT
*OSOS       PRINT * ,'KCONT(JCONT,NANAL),KCONT(ICONT,NANAL) = ',
*OSOS&                KCONT(JCONT,NANAL),KCONT(ICONT,NANAL),
*OSOS&        CHNAM(KCONT(ICONT,NANAL)),CHNAM(KCONT(JCONT,NANAL))
            IF(KCONT(ICONT,NANAL).EQ.KCONT(JCONT,NANAL)) THEN
              IF(IWORD.LE.NWORDS) THEN
                CALL COMBOS_ERROR(1,'Name '//CH16(:LENOCC(CH16))//
     &           ' appears more than once on MEASUREMENT logical line',
     &           'this MEASUREMENT logical line ignored')
              ELSE
                CALL COMBOS_ERROR(1,'Invalid use of name '//
     &           CH16(:LENOCC(CH16)),
     &           'this MEASUREMENT logical line ignored')
              ENDIF
*OSOS start
              IF(NQUANT.LE.1) THEN
                CALL UZERO(KCONT(0,NANAL),1,NCONT(NANAL)+1)
              ELSE
                KCONT(1-NQUANT,NANAL)=0
              ENDIF
              NQUANT=NQUANT-1
*OSOS end
              RETURN
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      END
*
************************************************************************
*
      SUBROUTINE STORE_DISP(NWORDS,WORDS)
*     ===================================
*
*     Routine called for the "DISPLAY" keyword: 
*     the name of the measured quantity, as well as the names of the 
*     statistical, total systematic, and total uncertainties 
*     to be displayed on plots are read in and stored
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORDS(*)
*
*     Externals
*
      INTEGER LENOCC
*
*     Local variable
*
      INTEGER IWORD,DUMMY,IDISP
      CHARACTER*16 CH16
*
      INCLUDE 'combos.inc'
*
*     Check that there is a measurement name and statistical error 
*
      IF(NWORDS.LE.0) THEN 
        CALL COMBOS_ERROR(1,'No measurement name specified',
     &                      'this DISPLAY logical line ignored')
        RETURN
      ENDIF
*
      IF(KDISP(0,NANAL).GT.0) CALL COMBOS_ERROR(0,
     & 'More than one DISPLAY keyword',
     & 'this logical line will overwrite the previous definition')
      IF(NCONT(NANAL).GT.FCONT-1) CALL COMBOS_ERROR(1,
     & 'DISPLAY keyword after DATA keyword(s)',
     & 'previous data (except step values) will be ignored')
      NCONT(NANAL)=FCONT-1
      NDISP(NANAL)=-1
      DO IDISP=0,MDISP
        KDISP(IDISP,NANAL)=0
        CALL UZERO(DISP(1,+IDISP,NANAL),1,MSTEP)
        CALL UZERO(DISP(1,-IDISP,NANAL),1,MSTEP)
      ENDDO
      DO IWORD=1,NWORDS
        CH16=WORDS(IWORD)
        IDISP=IWORD-1
        IF(IDISP.LE.MDISP) THEN
          IF(CH16.EQ.SUPERSEDED) THEN
            DO IDISP=0,MDISP
              KDISP(IDISP,NANAL)=-1
            ENDDO
            NDISP(NANAL)=0
          ELSE IF(KDISP(IDISP,NANAL).GE.0) THEN
            CALL CHOOSE_NAME(1,CH16,'display','DISPLAY',0)
            CALL CHECK_NAME(CH16,.TRUE.,KDISP(IDISP,NANAL),DUMMY)
            NDISP(NANAL)=IDISP
          ENDIF  
        ELSE
          CALL COMBOS_ERROR(1,
     &     'Too many words on DISPLAY logical line',
     &     'word '//CH16(:LENOCC(CH16))//' ignored') 
        ENDIF
      ENDDO
      END
*
************************************************************************
*
      SUBROUTINE STORE_DATA(NWORDS,WORDS)
*     ===================================
*
*     Routine called for the "DATA" keyword: 
*     the names and values of steps, measurements, statistical and 
*     systematic uncertainties are read in and stored
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORDS(*)
      CHARACTER*256 XWORDS
      INTEGER       LENOCC
*
*     Local variable
*
      INTEGER I,NVAL,N,K,KK
      CHARACTER*4 CH4
*
      INCLUDE 'combos.inc'
      INTEGER INDX(MWORDS),LSIGN(MWORDS),ICONT(MWORDS),IDISP(MWORDS)
      REAL VAL(MWORDS)
      LOGICAL LOAD_STEP
      LOGICAL STEP_LOADED,LOADED(0:MCONT,-1:+1),LOADIS(0:MDISP,-1:+1)
      SAVE    STEP_LOADED,LOADED               ,LOADIS
      INTEGER IANAL
      SAVE    IANAL
      DATA    IANAL/-1/
*RVK look for & store in standalone common
*DR look for % store in standalone common
*DR awful solution not to change all calling sequences
*DR indicate if error given in percentage. 
      REAL VALCOR
      LOGICAL YAMPER(MWORDS),YDOLLAR(MWORDS),YPERCENT(MWORDS)
      LOGICAL XPERCENT,XAMPER,XDOLLAR
      COMMON /PERCENT/XPERCENT,XAMPER,XDOLLAR
      XPERCENT=.FALSE.
      XAMPER=.FALSE.
      XDOLLAR=.FALSE.
*
*     Detect new analysis
*
      N=MCONT+1
      IF(NCONT(NANAL).EQ.FCONT-1) N=FCONT
      IF(NANAL.NE.IANAL) THEN
        IANAL=NANAL
        N=0
      ENDIF
      IF(N.LE.MCONT) THEN
        DO KK=-1,1
          DO K=N,MCONT
            LOADED(K,KK)=.FALSE.
          ENDDO
          DO K=0,MDISP
            LOADIS(K,KK)=.FALSE.
          ENDDO
        ENDDO
      ENDIF
*
*     Get list of items to load
*
      LOAD_STEP=.FALSE.
      DO I=1,NWORDS
        CALL CHECK_NAME(WORDS(I),.FALSE.,INDX(I),LSIGN(I))
        YAMPER(I) = XAMPER
        YDOLLAR(I) = XDOLLAR
*RVK-20050501 Need to store this per input word
        YPERCENT(I) = XPERCENT
        IF(LSIGN(I).NE.0.AND.
     &     (INDX(I).EQ.KCONT(0,NANAL).OR.INDX(I).EQ.KSTEP(NANAL))) THEN
          XWORDS = WORDS(I)
          CALL COMBOS_ERROR(0,'Invalid name '
     &     //XWORDS(:LENOCC(WORDS(I))),
     &     'name is taken as '//CHNAM(INDX(I)))
          WORDS(I)=CHNAM(INDX(I))
          LSIGN(I)=0
        ELSE IF(LSIGN(I).EQ.0.AND.
     &     (INDX(I).NE.KCONT(0,NANAL).AND.
     &      INDX(I).NE.KCONT(1,NANAL).AND.
     &      INDX(I).NE.KCONT(2,NANAL).AND.
     &      INDX(I).NE.KSTEP(NANAL))) THEN
****          CALL COMBOS_ERROR(0,'"Unsigned" name '//WORDS(I),
****     &     'assume symmetric errors for '//WORDS(I))
        ENDIF
        ICONT(I)=-1
        IDISP(I)=-1
        IF(INDX(I).NE.KSTEP(NANAL)) THEN
          DO K=0,NCONT(NANAL)
            IF(INDX(I).EQ.KCONT(K,NANAL)) ICONT(I)=K
          ENDDO
          DO K=0,NDISP(NANAL)
            IF(INDX(I).EQ.KDISP(K,NANAL)) IDISP(I)=K
          ENDDO
          IF(ICONT(I).EQ.-1.AND.IDISP(I).EQ.-1) THEN 
            IF(NCONT(NANAL).EQ.MCONT) CALL COMBOS_ERROR(-1,
     &       'Too many contributions','please increase parameter MCONT')
            NCONT(NANAL)=NCONT(NANAL)+1
            ICONT(I)=NCONT(NANAL)
            KCONT(ICONT(I),NANAL)=INDX(I)
          ENDIF
          IF(ICONT(I).NE.-1) THEN
            IF(LOADED(ICONT(I),LSIGN(I))) THEN
              XWORDS = WORDS(I)
              CALL COMBOS_ERROR(0,
     &        'Data for '//XWORDS(:LENOCC(WORDS(I)))//' already loaded',
     &        'previous data will be overwritten')
            ENDIF
            LOADED(ICONT(I),0)=.TRUE.
          ENDIF
          IF(IDISP(I).NE.-1) THEN
            IF(LOADIS(IDISP(I),LSIGN(I))) THEN
              XWORDS = WORDS(I)
              CALL COMBOS_ERROR(0,
     &        'Display data for '//XWORDS(:LENOCC(WORDS(I)))
     &        //' already loaded',
     &        'previous display data will be overwritten')
            ENDIF
            LOADIS(IDISP(I),0)=.TRUE.
          ENDIF
          IF(LSIGN(I).EQ.0) THEN
            DO KK=-1,1,2
              IF(ICONT(I).NE.-1) THEN
                CALL UZERO(CONT(1,KK*ICONT(I),NANAL),1,MSTEP)
                LOADED(ICONT(I),KK)=.TRUE.
              ENDIF
              IF(IDISP(I).NE.-1) THEN
                CALL UZERO(DISP(1,KK*IDISP(I),NANAL),1,MSTEP)
                LOADIS(IDISP(I),KK)=.TRUE.
              ENDIF
            ENDDO
          ELSE
            IF(ICONT(I).NE.-1) THEN
              CALL UZERO(CONT(1,LSIGN(I)*ICONT(I),NANAL),1,MSTEP)
              LOADED(ICONT(I),LSIGN(I))=.TRUE.
            ENDIF
            IF(IDISP(I).NE.-1) THEN
              CALL UZERO(DISP(1,LSIGN(I)*IDISP(I),NANAL),1,MSTEP)
              LOADIS(IDISP(I),LSIGN(I))=.TRUE.
            ENDIF
          ENDIF
        ELSE
          LOAD_STEP=.TRUE.
        ENDIF
      ENDDO
*
*     Load values
*
      N=0
    1 CALL GET_NEXT_VALUE(NWORDS,NVAL,VAL)
      IF(NVAL.EQ.0) THEN
        IF(N.EQ.0) THEN 
          CALL COMBOS_ERROR(0,'No numerical values provided',
     &                        'this logical line ignored')
          RETURN
        ENDIF
        IF(NSTEP(NANAL).EQ.0) THEN 
          NSTEP(NANAL)=N
        ELSE IF(N.NE.NSTEP(NANAL)) THEN
          NSTEP(NANAL)=MIN0(NSTEP(NANAL),N)
          WRITE(CH4,'(I4)',IOSTAT=K) NSTEP(NANAL)
          CALL COMBOS_ERROR(1,
     &         'Number of steps does not match expectation',
     &         'assume only '//CH4//' steps')
        ENDIF
        STEP_LOADED=STEP_LOADED.OR.LOAD_STEP
        RETURN
      ENDIF
      IF(N.EQ.MSTEP) CALL COMBOS_ERROR(-1,'Too many steps',
     &               'please increase parameter MSTEP')
      N=N+1
      DO I=1,NWORDS
        IF(ICONT(I).EQ.1.OR.ICONT(I).EQ.2.OR.ICONT(I).EQ.3.OR.
     &     IDISP(I).EQ.1.OR.IDISP(I).EQ.2.OR.IDISP(I).EQ.3) THEN ! stat. or tot. (syst.) uncertainty
          IF(FLOAT(LSIGN(I))*VAL(I).LT.0..OR.
     &       (LSIGN(I).EQ.0.AND.VAL(I).LT.0.)) THEN
            XWORDS = WORDS(I)
            CALL COMBOS_ERROR(1,
     &      XWORDS(:LENOCC(WORDS(I)))//' uncertainty specified with '
     &      //'the wrong sign',
     &      'change sign of value of '//XWORDS(:LENOCC(WORDS(I))))
            VAL(I)=-VAL(I)
          ENDIF
        ENDIF
        IF(INDX(I).EQ.KSTEP(NANAL)) THEN 
          IF(STEP_LOADED.AND.STEP(N,NANAL).NE.VAL(I))
     &      CALL COMBOS_ERROR(-1,'Unexpected step value',
     &                           'please correct input file')
          STEP(N,NANAL)=VAL(I)
        ELSE

          IF (YDOLLAR(I)) THEN
             CONT(N,+MCONT+ICONT(I),NANAL) = 0.5
             CONT(N,-MCONT-ICONT(I),NANAL) = 0.5
          ENDIF
          IF (YAMPER(I)) THEN
             CONT(N,+MCONT+ICONT(I),NANAL) = 1
             CONT(N,-MCONT-ICONT(I),NANAL) = 1
          ENDIF

*DR treat value given as percentage
          VALCOR=VAL(I)
          IF (YPERCENT(I)) THEN
            IF (CONT(N,0,NANAL).LE.0.) THEN
              CALL COMBOS_ERROR(0,
     &            'percentage error on neg or null value',
     &                           'certainly wrong')
            ENDIF
            VALCOR=VAL(I)/100.*CONT(N,0,NANAL)
          ENDIF
*DR end
*DR then replace VAL(I) by VALCOR
          IF(ICONT(I).NE.-1) THEN
            IF(INDX(I).EQ.KCONT(0,NANAL)) THEN 
              IF(ICONT(I).NE.0) PRINT * ,'problem: ',ICONT(I)
              CONT(N,0,NANAL)=VALCOR
            ELSE IF(LSIGN(I).EQ.0) THEN
              CONT(N,+ICONT(I),NANAL)=+VALCOR
              CONT(N,-ICONT(I),NANAL)=-VALCOR
            ELSE
              CONT(N,LSIGN(I)*ICONT(I),NANAL)=VALCOR
            ENDIF
          ENDIF
          IF(IDISP(I).NE.-1) THEN
            IF(INDX(I).EQ.KDISP(0,NANAL)) THEN
              IF(IDISP(I).NE.0) PRINT * ,'problem: ',IDISP(I)
              DISP(N,0,NANAL)=VALCOR
            ELSE IF(LSIGN(I).EQ.0) THEN
              DISP(N,+IDISP(I),NANAL)=+VALCOR
              DISP(N,-IDISP(I),NANAL)=-VALCOR
            ELSE
              DISP(N,LSIGN(I)*IDISP(I),NANAL)=VALCOR
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      GOTO 1
*
      ENTRY STORE_DATA_RESET
      IANAL=-1
      STEP_LOADED=.FALSE.
      RETURN
*
      ENTRY STORE_DATA_STEP_LOADED
      STEP_LOADED=.TRUE.
      END
*
************************************************************************
*
      SUBROUTINE STORE_LINK(KEYWORD,NWORDS,WORDS)
*     ===========================================
*
*     Routine called for the KEYWORD keyword ('SUPERSEDE' or 'STAT_CORR_WITH'): 
*     links to other analyses are read in and stored
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORDS(*),KEYWORD
      CHARACTER*256 XWORDS,XKEYWORD
*
*     Externals
*
      INTEGER LENOCC
*
*     Local variable
*
      INTEGER IWORD,KEXPER,KMETHO,KSTATU,DUMMY,NEWLINK,I,NRHO
      REAL RHO
      CHARACTER*16 CH16
      CHARACTER*80 CH80
*
      INCLUDE 'combos.inc'
*
*     Check number of arguments
*
      IF(KEYWORD.EQ.'STAT_CORR_WITH') THEN
        DO IWORD=4,NWORDS
          XKEYWORD = KEYWORD
          XWORDS   = WORDS(IWORD)
          CALL COMBOS_ERROR(1,'Too many words on '
     &    //XKEYWORD(:LENOCC(KEYWORD))//
     &    ' logical line','word '//XWORDS(:LENOCC(WORDS(IWORD)))
     &    //' ignored')
        ENDDO
        NWORDS=MIN0(NWORDS,3)
      ENDIF
*
*     Loop on new links
*
      IWORD=0
    1 IF(NWORDS-IWORD.LE.0) RETURN
*
*     Get analysis to link
*
      IWORD=IWORD+1
      IF(IWORD.GT.NWORDS) WORDS(IWORD)=' '
      CALL CHECK_CHAR(WORDS(IWORD),MEXP,NEXP,CHEXP,'experiment','MEXP',
     &                .TRUE.,KEXPER,DUMMY)
      IWORD=IWORD+1
      IF(IWORD.GT.NWORDS) WORDS(IWORD)=' '
      CALL CHECK_CHAR(WORDS(IWORD),MMETH,NMETH,CHMETH,'method','MMETH',
     &                .TRUE.,KMETHO,DUMMY)
      IWORD=IWORD+1
      IF(IWORD.GT.NWORDS) WORDS(IWORD)=' '
      CALL CHECK_CHAR(WORDS(IWORD),MSTAT,NSTAT,CHSTAT,'status','MSTAT',
     &                .TRUE.,KSTATU,DUMMY)
      NEWLINK=1000000*KEXPER+1000*KMETHO+KSTATU
      CH80=CHEXP (KEXPER)(:LENOCC(CHEXP (KEXPER)))//' '//
     &     CHMETH(KMETHO)(:LENOCC(CHMETH(KMETHO)))//' '//
     &     CHSTAT(KSTATU)(:LENOCC(CHSTAT(KSTATU)))
*
*     Find statistical correlation coefficient
*
      IF(KEYWORD.EQ.'SUPERSEDE') THEN 
        RHO=9999. ! Note: a correlation of 9999 means that analysis NANAL
*                         supersedes analysis KEXPER,KMETHO,KSTATU
      ELSE
        CALL GET_NEXT_VALUE(1,NRHO,RHO)
        IF(NRHO.LE.0) THEN 
          CALL COMBOS_ERROR(1,'No statistical correlation '//
     &     'coefficient specified for '//CH80(:LENOCC(CH80)),
     &     'correlation set to 0')
          RHO=0.
        ELSE IF(RHO.LT.-1.) THEN
          CALL COMBOS_ERROR(1,'Invalid statistical correlation '//
     &     'coefficient specified for '//CH80(:LENOCC(CH80)),
     &     'correlation set to -1')
          RHO=-1.
C        ELSE IF(RHO.LT.0.) THEN
C          CALL COMBOS_ERROR(0,'Negative correlation '//
C     &     'coefficient specified for '//CH80(:LENOCC(CH80)),
C     &     'probably OK if 2D fit is done')
        ELSE IF(RHO.GT.1.) THEN
          CALL COMBOS_ERROR(1,'Invalid statistical correlation '//
     &     'coefficient specified for '//CH80(:LENOCC(CH80)),
     &     'correlation set to 1')
          RHO=1.
        ENDIF
      ENDIF
*
*     Check that analysis is not already linked
*
      DO I=0,NLINK(NANAL)
        IF(NEWLINK.EQ.LINK(I,NANAL)) THEN
          CH80=CHEXP (KEXPER)(:LENOCC(CHEXP (KEXPER)))//' '//
     &         CHMETH(KMETHO)(:LENOCC(CHMETH(KMETHO)))//' '//
     &         CHSTAT(KSTATU)(:LENOCC(CHSTAT(KSTATU)))
          IF(I.EQ.0) THEN
            CH16='BEGIN'
          ELSE IF(ABS(RHOSTA(I,NANAL)).GT.1.) THEN
            CH16='SUPERSEDE'
          ELSE
            CH16='STAT_CORR_WITH'
          ENDIF
          CALL COMBOS_ERROR(1,'Analysis '//CH80(:LENOCC(CH80))//
     &        ' already specified after '//CH16(:LENOCC(CH16))//
     &        ' keyword','this logical line ignored')
           GOTO 1
        ENDIF
      ENDDO
*
*     Create new link
*
      IF(NLINK(NANAL).EQ.MLINK) CALL COMBOS_ERROR(-1,
     & 'Too many links','please increase parameter MLINK')
      NLINK(NANAL)=NLINK(NANAL)+1
      LINK(NLINK(NANAL),NANAL)=NEWLINK
      RHOSTA(NLINK(NANAL),NANAL)=RHO
      GOTO 1
      END
*
************************************************************************
*
      SUBROUTINE STORE_COMB(NWORDS,WORDS)
*     ===================================
*
*     Routine called for the "COMB" keyword: 
*     the names of the analyses to combine are read in and stored
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORDS(*)
      CHARACTER*256 XWORDS
      INTEGER LENOCC
*
*     Arguments of entry STORE_COMB_GET
*
      INTEGER IANAL
*
*     Externals
*
      INTEGER ICNTH
*
*     Local variable
*
      INTEGER IWORD,I,ICOUNT
      LOGICAL MARK
      INTEGER LASTANAL
      SAVE    LASTANAL
      DATA    LASTANAL/0/
*
      INCLUDE 'combos.inc'
*
      XCOMB(0)=.TRUE.
*
      LASTANAL=0
      IWORD=0
    1 IF(NWORDS-IWORD.LE.0) RETURN
*
*     Check value of arguments
*
      IF(IWORD+1.GT.NWORDS) WORDS(IWORD+1)=' '
      IF(WORDS(IWORD+1).NE.'*'.AND.ICNTH(WORDS(IWORD+1),CHEXP,NEXP)
     &.EQ.0) THEN
         XWORDS = WORDS(IWORD+1)
         CALL COMBOS_ERROR(1,
     &   'Invalid experiment name specified after COMB keyword: '//
     &   XWORDS(:LENOCC(WORDS(IWORD+1))),' ')
      ENDIF
      IF(IWORD+2.GT.NWORDS) WORDS(IWORD+2)=' '
      IF(WORDS(IWORD+2).NE.'*'.AND.ICNTH(WORDS(IWORD+2),CHMETH,NMETH)
     &.EQ.0) THEN
        XWORDS = WORDS(IWORD+2)
        CALL COMBOS_ERROR(1,
     &  'Invalid method name specified after COMB keyword: '//
     &  XWORDS(:LENOCC(WORDS(IWORD+2))),' ')
      ENDIF
      IF(IWORD+3.GT.NWORDS) WORDS(IWORD+3)=' '
      IF(WORDS(IWORD+3).NE.'*'.AND. ICNTH(WORDS(IWORD+3),CHSTAT,NSTAT)
     &.EQ.0) THEN
        XWORDS = WORDS(IWORD+3)
        CALL COMBOS_ERROR(1,
     &  'Invalid status name specified after COMB keyword: '//
     &  XWORDS(:LENOCC(WORDS(IWORD+3))),' ')
      ENDIF
*        
*     Mark analyses to be combined
*
      LASTANAL=0
      ICOUNT=0
      DO I=1,NANAL-1
        MARK=
     &(WORDS(IWORD+1).EQ.'*'.OR.WORDS(IWORD+1).EQ.CHEXP (KEXP (I))).AND.
     &(WORDS(IWORD+2).EQ.'*'.OR.WORDS(IWORD+2).EQ.CHMETH(KMETH(I))).AND.
     &(WORDS(IWORD+3).EQ.'*'.OR.WORDS(IWORD+3).EQ.CHSTAT(KSTAT(I)))
        XCOMB(I)=XCOMB(I).OR.MARK
        IF(MARK) THEN
          IF(LASTANAL.EQ.0) THEN 
            LASTANAL=I
          ELSE
            LASTANAL=-1
          ENDIF
        ENDIF
      ENDDO
      IWORD=IWORD+3
      GOTO 1
*
      ENTRY STORE_COMB_GET(IANAL)
      IANAL=LASTANAL
      END
*
************************************************************************
*
      SUBROUTINE STORE_LUMP(NWORDS,WORDS)
*     ===================================
*
*     Routine called for the "LUMP" keyword
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORDS(*)
      CHARACTER XWORDS
*
*     Externals
*
      INTEGER LENOCC
*
*     Local variables
*
      INTEGER IANAL,ILUMP,DUMMY,I,ICONT,INAM,ICOUNT,NLUMP
      CHARACTER*80 NAME
*
      INCLUDE 'combos.inc'
*
*     Check that there is an unambiguous analysis
*     associated to LUMP logical line
*
      CALL STORE_COMB_GET(IANAL)
      IF(IANAL.LE.0.OR.IANAL.GE.NANAL) THEN 
        IF(IANAL.EQ.0) THEN
          CALL COMBOS_ERROR(1,
     &      'Keyword LUMP does not apply to any analysis',
     &      'this logical line ignored')
        ELSE
          CALL COMBOS_ERROR(1,
     &      'Keyword LUMP does apply to several analyses',
     &      'this ambiguous logical line ignored')
        ENDIF
        RETURN
      ENDIF
*
*     Check number of tokens on logical line
*
      IF(NWORDS.LT.1) THEN
        CALL COMBOS_ERROR(1,'Missing lump name after LUMP keyword',
     &                      'this logical line ignored')
        RETURN
      ENDIF
*
*     Check lump name
*
      IF(WORDS(1)(11:).NE.' ') THEN
        XWORDS = WORDS(1)
        CALL COMBOS_ERROR(1,
     &      'Lump name too long: '//XWORDS(:LENOCC(WORDS(1))),
     &      'this LUMP logical line ignored')
        RETURN
      ENDIF
      WORDS(1)='lump '//WORDS(1)
      CALL CHECK_NAME(WORDS(1),.TRUE.,ILUMP,DUMMY)
*
*     Check contributions included in lump
*
      NLUMP=0
      DO I=2,NWORDS
        CALL CHECK_NAME(WORDS(I),.TRUE.,INAM,DUMMY)
        ICOUNT=0
        DO ICONT=FCONT,NCONT(IANAL)
          IF(KSYN(KCONT(ICONT,IANAL)).EQ.KSYN(INAM)) THEN
            ICOUNT=ICOUNT+1
            IF(KLUMP(ICONT,IANAL).NE.KCONT(ICONT,IANAL).AND.
     &         KLUMP(ICONT,IANAL).NE.ILUMP) THEN
              CALL ANALYSIS_NAME(.TRUE.,IANAL,NAME)
              CALL COMBOS_ERROR(1,'Contribution '//
     & CHNAM(INAM)(:LENOCC(CHNAM(INAM)))//' in '//
     & NAME(:LENOCC(NAME))//' already included in '//
     & CHNAM(KLUMP(ICONT,IANAL))(:LENOCC(CHNAM(KLUMP(ICONT,IANAL)))),
     & 'this contribution not included in '//
     & CHNAM(ILUMP)(:LENOCC(CHNAM(ILUMP))))
            ELSE
              KLUMP(ICONT,IANAL)=ILUMP
              NLUMP=NLUMP+1
            ENDIF
          ENDIF
        ENDDO
        IF(ICOUNT.LE.0) THEN 
          CALL ANALYSIS_NAME(.TRUE.,IANAL,NAME)
          CALL COMBOS_ERROR(1,
     & 'No contribution '//CHNAM(INAM)(:LENOCC(CHNAM(INAM)))//' in '//
     & NAME(:LENOCC(NAME)),'this contribution not included in '//
     & CHNAM(ILUMP)(:LENOCC(CHNAM(ILUMP))))
        ENDIF
        IF(ICOUNT.GT.1) THEN 
          CALL ANALYSIS_NAME(.TRUE.,IANAL,NAME)
          CALL COMBOS_ERROR(-1,'Several contributions '//
     &     CHNAM(INAM)(:LENOCC(CHNAM(INAM)))//' in '//
     &     NAME(:LENOCC(NAME)),'please report this problem')
        ENDIF
      ENDDO
      IF(NLUMP.LE.0) THEN 
        CALL ANALYSIS_NAME(.TRUE.,IANAL,NAME)
        CALL COMBOS_ERROR(0,
     &   'Empty '//CHNAM(ILUMP)(:LENOCC(CHNAM(ILUMP)))//' for '//
     &   NAME(:LENOCC(NAME)),'this lump ignored')
      ENDIF
      END
*
************************************************************************
*
      SUBROUTINE STORE_CALL(NWORDS,WORDS)
*     ===================================
*
*     Routine called for the "CALL" keyword: 
*     the names of the routines to be called (to perform the combination) 
*     are read in and stored
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORDS(*)
      CHARACTER*256 XWORDS
*
*     Externals
*
      INTEGER ICNTH
      INTEGER LENOCC
*
*     Local variable
*
      INTEGER IWORD,ICALL
      INTEGER FCALLP,GCALLP,FCALLC,GCALLC,FCALLT,GCALLT
      INTEGER I,J,NVAL,LUNIT
      REAL VAL
*
      INCLUDE 'combos.inc'
      INTEGER MCALLP1,MCALLC1,MCALL
      PARAMETER(MCALLP1=MCALLP+1,MCALLC1=MCALLC+MCALLP1,
     &          MCALL=MCALLP+MCALLC+MCALLT)
      CHARACTER*16 CHCALL(MCALL) ! names of routines to call
*
*      EQUIVALENCE (CHCALL(1),CHCALP(1))
*      EQUIVALENCE (CHCALL(MCALLP1),CHCALC(1))
*      EQUIVALENCE (CHCALL(MCALLC1),CHCALT(1))
*
      DO I=1,MCALLP
        CHCALL(I)=CHCALP(I)
      ENDDO
      DO I=1,MCALLC
        CHCALL(I+MCALLP)=CHCALC(I)
      ENDDO
      DO I=1,MCALLT
        CHCALL(I+MCALLP+MCALLC)=CHCALT(I)
      ENDDO

      XCOMB(0)=.TRUE.
*
*     Get logical unit for printout
*
      CALL GET_NEXT_VALUE(1,NVAL,VAL)
      IF(NVAL.EQ.0) THEN 
        LUNIT=LUNLOG
      ELSE
        LUNIT=NINT(VAL)
      ENDIF
      IF(IABS(LUNIT).EQ.5.OR.IABS(LUNIT).GT.89) THEN
        CALL COMBOS_ERROR(1,'Invalid logical unit',
     &                      'this logical unit ignored')
        LUNIT=LUNLOG
      ENDIF
*
      DO IWORD=1,NWORDS
*
*       Check value of arguments
*
        IF(WORDS(IWORD).NE.'*') THEN
          FCALLP=0
          GCALLP=-1
          FCALLC=0
          GCALLC=-1
          FCALLT=0
          GCALLT=-1
          I=ICNTH(WORDS(IWORD),CHCALL,MCALL)
          IF(I.EQ.0) THEN
            XWORDS = WORDS(IWORD)
            CALL COMBOS_ERROR(1,
     &       'Invalid routine name specified after CALL keyword: '//
     &       XWORDS(:LENOCC(WORDS(IWORD))),'this name ignored')
          ELSE IF(I.LE.MCALLP) THEN
            FCALLP=I
            GCALLP=I
          ELSE IF(I.LE.MCALLP+MCALLC) THEN
            FCALLC=I-MCALLP
            GCALLC=I-MCALLP
          ELSE IF(I.LE.MCALL) THEN
            FCALLT=I-MCALLP-MCALLC
            GCALLT=I-MCALLP-MCALLC
          ENDIF
        ELSE
          FCALLP=1
          GCALLP=MCALLP
          FCALLC=1
          GCALLC=MCALLC
          FCALLT=1
          GCALLT=MCALLT
        ENDIF
*
*       Add preparation routine(s) to the list of routines to be called
*
        DO 1 ICALL=FCALLP,GCALLP
          J=0
          DO I=1,NCALLP
            IF(KCALLP(I).EQ.ICALL) J=I
          ENDDO
          IF(J.EQ.0) THEN
            NCALLP=NCALLP+1
            IF(NCALLP.GT.MCALLP) CALL COMBOS_ERROR(-1,
     &       'Too many names of preparation routines',
     &       'please increase parameter MCALLP')
            J=NCALLP
          ENDIF
          KCALLP(J)=ICALL
          LCALLP(J)=LUNIT
    1   CONTINUE
*
*       Add combination routine(s) to the list of routines to be called
*
        DO 2 ICALL=FCALLC,GCALLC
          J=0
          DO I=1,NCALLC
            IF(KCALLC(I).EQ.ICALL) J=I
          ENDDO
          IF(J.EQ.0) THEN
            NCALLC=NCALLC+1
            IF(NCALLC.GT.MCALLC) CALL COMBOS_ERROR(-1,
     &       'Too many names of combination routines',
     &       'please increase parameter MCALLC')
            J=NCALLC
          ENDIF
          KCALLC(J)=ICALL
          LCALLC(J)=LUNIT
    2   CONTINUE
*
*       Add termination routine(s) to the list of routines to be called
*
        DO 3 ICALL=FCALLT,GCALLT
          J=0
          DO I=1,NCALLT
            IF(KCALLT(I).EQ.ICALL) J=I
          ENDDO
          IF(J.EQ.0) THEN
            NCALLT=NCALLT+1
            IF(NCALLT.GT.MCALLT) CALL COMBOS_ERROR(-1,
     &       'Too many names of termination routines',
     &       'please increase parameter MCALLT')
            J=NCALLT
          ENDIF
          KCALLT(J)=ICALL
          LCALLT(J)=LUNIT
    3   CONTINUE
      ENDDO
      END
*
************************************************************************
*
      SUBROUTINE STORE_SYNO(NWORDS,WORDS)
*     ===================================
*
*     Routine called for the "SYNO" keyword
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORDS(*)
*
*     Local variables
*
      INTEGER I,J,INAM,DUMMY,ISYN1,ISYN
*
      INCLUDE 'combos.inc'
*
      XCOMB(0)=.TRUE.
*
      DO I=1,NWORDS
        CALL CHECK_NAME(WORDS(I),.TRUE.,INAM,DUMMY)
        ISYN=KSYN(INAM)
        IF(I.EQ.1) THEN
          ISYN1=ISYN
        ELSE
          DO J=1,NNAM
            IF(KSYN(J).EQ.ISYN) KSYN(J)=ISYN1
          ENDDO
        ENDIF
      ENDDO
      END
*
************************************************************************
*
      SUBROUTINE STORE_SWIT(NWORDS,WORDS)
*     ===================================
*
*     Routine called for the "SWIT" keyword; switch various options 
*     "ON" or "OFF".
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORDS(*)
      CHARACTER*256 XWORDS
*
*     Externals
*
      INTEGER ICNTH,LENOCC
*
*     Local variables
*
      INTEGER NVAL,I,NSW,J,K
      REAL VAL
      LOGICAL ALL
*
      INCLUDE 'combos.inc'
      INTEGER MSW3
      PARAMETER(MSW3=MSW+3)
      INTEGER ISW(MSW3)
      CHARACTER*16 CHSW(MSW3)
      LOGICAL SWDEF(MSW3)
*                                      Switch name      Default value
      DATA CHSW(LADJU),SWDEF(LADJU)/'ADJU*STMENTS    ',.TRUE./
      DATA CHSW(LSYMA),SWDEF(LSYMA)/'SYMM_ADJU*      ',.FALSE./
      DATA CHSW(LSYMC),SWDEF(LSYMC)/'SYMM_COMB*      ',.FALSE./
      DATA CHSW(LSYST),SWDEF(LSYST)/'SYST*EMATICS    ',.TRUE./
      DATA CHSW(LCOSY),SWDEF(LCOSY)/'SYST_CORR*      ',.TRUE./
      DATA CHSW(LCOST),SWDEF(LCOST)/'STAT_CORR*      ',.TRUE./
      DATA CHSW(LINTE),SWDEF(LINTE)/'INTE*RPOLATION  ',.TRUE./
      DATA CHSW(LEXTR),SWDEF(LEXTR)/'EXTR*APOLATION  ',.FALSE./
      DATA CHSW(LQUST),SWDEF(LQUST)/'QUOT*E_STAT     ',.TRUE./
      DATA CHSW(LLUMP),SWDEF(LLUMP)/'LUMP*S          ',.TRUE./
      DATA (CHSW(I),SWDEF(I),I=MSW+1,MSW3)/'OFF      ',.FALSE.,
     &                                     'ON       ',.TRUE.,
     &                                     'TOGGLE   ',.TRUE./
*RVK add switches SYSREL and SYSABS and STATREL and STATABS...and MAXCORR
      DATA CHSW(LSYRE),SWDEF(LSYRE)/'SYSREL*         ',.FALSE./
      DATA CHSW(LSYAB),SWDEF(LSYAB)/'SYSABS*         ',.FALSE./
      DATA CHSW(LSTSQ),SWDEF(LSTSQ)/'STATSQRT*       ',.FALSE./
      DATA CHSW(LSTAB),SWDEF(LSTAB)/'STATABS*        ',.FALSE./
      DATA CHSW(LMXCR),SWDEF(LMXCR)/'MAXCORR*        ',.TRUE./
*
      XCOMB(0)=.TRUE.
*
      CALL GET_NEXT_VALUE(1,NVAL,VAL)
      IF(NWORDS.LE.0) RETURN
      NSW=1
      ALL=.FALSE.
      DO I=1,NWORDS
        IF(WORDS(I).EQ.'*') THEN 
          ALL=.TRUE.
          DO J=1,MSW
            ISW(J)=J
          ENDDO
          NSW=MSW+1
        ELSE
          ISW(NSW)=ICNTH(WORDS(I),CHSW,MSW3)
          IF(ISW(NSW).LE.0) THEN 
            XWORDS = WORDS(I)
            CALL COMBOS_ERROR(1,'Invalid switch name',
     &       'string '//XWORDS(:LENOCC(WORDS(I)))//' ignored')
          ELSE IF(ISW(NSW).GT.MSW) THEN
            DO J=1,NSW-1
              IF(CHSW(ISW(NSW)).EQ.'TOGGLE') THEN
                XON(ISW(J))=.NOT.XON(ISW(J))
              ELSE
                XON(ISW(J))=SWDEF(ISW(NSW))
              ENDIF
            ENDDO
            NSW=1
            ALL=.FALSE.
          ELSE IF(.NOT.ALL) THEN
            IF(NSW.GE.MSW3) CALL COMBOS_ERROR(-1,
     &        'Too many tokens on SWITCH logical line',
     &        'please increase parameter MSW3')
            NSW=NSW+1
          ENDIF
        ENDIF
      ENDDO
      NSW=NSW-1
      DO J=1,NSW
        IF(NVAL.GT.0) THEN
          XON(ISW(J))=NINT(VAL).NE.0
        ELSE
          XON(ISW(J))=SWDEF(ISW(J)) ! default used if no value specified
        ENDIF
      ENDDO
      IF (XON(LSYRE).AND.XON(LSYAB))
     &  STOP 'SYSREL and SYSABS are incompatible options'
      IF (XON(LSTSQ).AND.XON(LSTAB))
     &  STOP 'STATSQRT and STATABS are incompatible options'

      RETURN
*
      ENTRY RESET_SWIT
*     ----------------
      DO I=1,MSW
        XON(I)=SWDEF(I)
      ENDDO
      RETURN
*
      ENTRY PRINT_SWIT
*     ----------------
      DO I=1,MSW
        IF(XON(I).NEQV.SWDEF(I)) THEN
          DO J=MSW+1,MSW+2
            IF(XON(I).EQV.SWDEF(J)) THEN
              K=INDEX(CHSW(I),'*')
              WRITE(LUNOUT,1001) CHSW(I)(:K-1)//CHSW(I)(K+1:),CHSW(J)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
 1001 FORMAT('SWITCH',T17,6(1X,A))
      END
*
************************************************************************
*
      SUBROUTINE STORE_MINO(NWORDS,WORDS)
*     ===================================
*
*     Routine called for the "MINO" keyword
*     HCJS July 24 1997
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS
      CHARACTER*(*) WORDS(*)
*
*     Local variables
*
      INTEGER i
*
      INCLUDE 'combos.inc'
*
c$$$      LOGICAL FOUND
c$$$      INTEGER ICINQ
c$$$      INTEGER JMEAS,JSYS

      DO I=1,NWORDS
        IF(I.LE.NWORDS) THEN
          CHMINOS(I)=WORDS(I)
        ELSE
          CHMINOS(I)=' '
        ENDIF
      ENDDO
      XMINOS = .TRUE.
      NMINOS = NWORDS

c$$$      XMINOS = .FALSE.
c$$$*     Check that the input for minos errors makes sense.
c$$$      IF (nwords .EQ. 1 .AND. words(1) .EQ. 'ALL') THEN
c$$$        XMINOS = .TRUE.
c$$$      ELSE
c$$$*       Check that all parameters for which minos is requested are known
c$$$*       to the system.
c$$$        k = 1
c$$$        DO 100, I = 1, NWORDS
c$$$          IF (WORDS(I) .EQ. 'ALL') THEN
c$$$            K = 1
c$$$            CHMINOS(K) = 'ALL'
c$$$            XMINOS = .TRUE.
c$$$            GOTO 1000
c$$$          ELSE
c$$$            FOUND = .FALSE.
c$$$            JMEAS = ICINQ(WORDS(I),CHMEAS,1)
c$$$            IF (JMEAS .GT. 0) THEN
c$$$              FOUND = .TRUE.
c$$$              K = K + 1
c$$$              CHMINOS(K) = WORDS(I)
c$$$              XMINOS = .TRUE.
c$$$              GOTO 100
c$$$            ENDIF
c$$$            JSYS = ICINQ(WORDS(I),CHNAM,NCSYS)
c$$$            IF (JSYS .GT. 0) THEN
c$$$              FOUND = .TRUE.
c$$$              K = K + 1
c$$$              CHMINOS(K) = WORDS(I)
c$$$              XMINOS = .TRUE.
c$$$              GOTO 100
c$$$            ENDIF
c$$$            IF (.NOT. FOUND) THEN
c$$$              XMINOS = .FALSE.
c$$$              GOTO 1000
c$$$            ENDIF
c$$$          ENDIF
c$$$ 100    ENDDO
c$$$      ENDIF
c$$$
c$$$      
c$$$ 1000 CONTINUE 
c$$$      print *, 'names', chminos
c$$$      print *, 'mino flag', XMINOS
*
      END
*
************************************************************************
*
      SUBROUTINE CHECK_ANAL(NWORDS,WORDS,IERR)
*     ========================================
*
*     Routine called for the "END" keyword in case of an "input" analysis 
*     (i.e. an analysis defined with DATA keywords): 
*     various checks are preformed now that all the information for 
*     the analysis have been read in 
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS,IERR
      CHARACTER*(*) WORDS(*)
*
*     Externals
*  
      INTEGER LENOCC
*
*     Local variable
*
      INTEGER I,K,KK,L1,L2,NSUP,IDISP,JDISP,ICALL
      INTEGER DUMMY,NZERO(-1:+1),NSMALL,NLARGE,NEWPAR,IPAR,ICONT
      LOGICAL REPLACED(-1:1),WARNREPL
      DATA WARNREPL/.TRUE./
      REAL SYST(-1:1,-1:1),SL1(-1:1),SL2(-1:1),RHO,ETOT(-1:1)
      CHARACTER*16 CH16
*
      INCLUDE 'combos.inc'
*
*     Check steps
*
      CALL CHECK_STEPS(1,NANAL,IERR)
*
*     Check that a single measurement name has been specified
*
*OSOS start
      IF(NQUANT.NE.1) THEN
        IERR=IERR+1
        IF(NQUANT.LE.0) THEN
          CALL COMBOS_ERROR(1,'No MEASUREMENT keyword provided',
     &                        'this analysis will be ignored')
        ELSE
          CALL COMBOS_ERROR(1,'Several MEASUREMENT keywords provided',
     &                        'this analysis will be ignored')
        ENDIF
      ENDIF
*OSOS end
*
*     Do not allow two analyses with same experiment,
*     method, status and measurement
*
      DO I=1,NANAL-1
        IF(KSTAT(I).EQ.KSTAT(NANAL).AND.KMETH(I).EQ.KMETH(NANAL).AND.
     &     KEXP(I).EQ.KEXP(NANAL).AND.KCONT(0,I).EQ.KCONT(0,NANAL)) THEN
          CALL COMBOS_ERROR(1,'Another analysis already uses the '//
     &     'same experiment, method, status and measurement names',
     &     'this analysis will be ignored')
          IERR=IERR+1
        ENDIF
      ENDDO
*
*     Check statistical errors
*
      DO K=-1,1,2
        NZERO(K)=0
      ENDDO
      DO I=1,NSTEP(NANAL)
        DO K=-1,1
          IF(CONT(I,K*1,NANAL).EQ.0.) NZERO(K)=NZERO(K)+1
        ENDDO
      ENDDO
      DO K=-1,1,2
        IF(NZERO(K).EQ.NSTEP(NANAL)) THEN 
          WRITE(CH16,'(A,SP,I2)',IOSTAT=DUMMY) 
     &     CHNAM(KCONT(1,NANAL))(:LENOCC(CHNAM(KCONT(1,NANAL)))),K
          CALL COMBOS_ERROR(0,
     &     'No value provided for '//CH16(:LENOCC(CH16)-1),
     &     'assume no statistical uncertainty')
        ENDIF
      ENDDO
*
*     Compute total systematic errors
*
      IF(KCONT(2,NANAL).GT.0) THEN
      DO K=-1,1,2
        NZERO(K)=0
        REPLACED(K)=.FALSE.
      ENDDO
      DO I=1,NSTEP(NANAL)
        DO K=-1,1
          IF(CONT(I,K*2,NANAL).EQ.0.) NZERO(K)=NZERO(K)+1
        ENDDO
      ENDDO
      NSMALL=0
      NLARGE=0
      DO I=1,NSTEP(NANAL)
        DO KK=-1,1
        DO K=-1,1,2
          SYST(K,KK)=0.
        ENDDO
        ENDDO
        DO L1=FCONT,NCONT(NANAL)
          SL1(+1)=
     &      ABS(AMAX1(0.,AMAX1(CONT(I,-L1,NANAL),CONT(I,L1,NANAL))))
          SL1(-1)=
     &      ABS(AMIN1(0.,AMIN1(CONT(I,-L1,NANAL),CONT(I,L1,NANAL))))
        DO L2=FCONT,NCONT(NANAL)
          SL2(+1)=
     &      ABS(AMAX1(0.,AMAX1(CONT(I,-L2,NANAL),CONT(I,L2,NANAL))))
          SL2(-1)=
     &      ABS(AMIN1(0.,AMIN1(CONT(I,-L2,NANAL),CONT(I,L2,NANAL))))
          DO KK=-1,1
            IF(L1.EQ.L2) THEN 
              RHO=1.
            ELSE
              RHO=FLOAT(KK)
            ENDIF
          DO K=-1,1,2
            SYST(K,KK)=SYST(K,KK)+SL1(K)*SL2(K)*RHO
          ENDDO
          ENDDO
        ENDDO
        ENDDO
        DO K=-1,1,2
          DO KK=-1,1
            IF(KK.EQ.-1.AND.SYST(K,KK).LT.0.) THEN 
              SYST(K,KK)=0.
            ELSE
              SYST(K,KK)=FLOAT(K)*SQRT(SYST(K,KK))
            ENDIF
          ENDDO
          IF(NZERO(K).EQ.NSTEP(NANAL)) THEN
            REPLACED(K)=REPLACED(K).OR.CONT(I,K*2,NANAL).NE.SYST(K,0)
            CONT(I,K*2,NANAL)=SYST(K,0)
          ENDIF
          IF(ABS(CONT(I,K*2,NANAL)).GT.ABS(SYST(K,+1))) NLARGE=NLARGE+1
          IF(ABS(CONT(I,K*2,NANAL)).LT.ABS(SYST(K,-1))) NSMALL=NSMALL+1
        ENDDO
      ENDDO
      IF(WARNREPL) THEN
        DO K=-1,1,2
          IF(REPLACED(K)) THEN
            WRITE(CH16,'(A,SP,I2)',IOSTAT=DUMMY) 
     &       CHNAM(KCONT(2,NANAL))(:LENOCC(CHNAM(KCONT(2,NANAL)))),K
            CALL COMBOS_ERROR(1,'No value provided for '//
     &       CH16(:LENOCC(CH16)-1),'will compute total systematic '//
     &       'uncertainty from individual contributions')
          ENDIF
        ENDDO
      ENDIF
      IF(NSMALL.GT.0) CALL COMBOS_ERROR(1,'Total systematic error '//
     & 'is sometimes smaller than sum of individual contributions',' ')
      IF(NLARGE.GT.0) CALL COMBOS_ERROR(1,'Total systematic error '//
     & 'is sometimes larger than sum of individual contributions',' ')
      ENDIF
*
*     Compute total (statistical+systematic) errors
*
      IF(KCONT(3,NANAL).GT.0) THEN
      DO K=-1,1,2
        NZERO(K)=0
        REPLACED(K)=.FALSE.
      ENDDO
      DO I=1,NSTEP(NANAL)
        DO K=-1,1
          IF(CONT(I,K*3,NANAL).EQ.0.) NZERO(K)=NZERO(K)+1
        ENDDO
      ENDDO
      NSMALL=0
      NLARGE=0
      DO I=1,NSTEP(NANAL)
        DO K=-1,1,2
          ETOT(K)=FLOAT(K)*SQRT(CONT(I,K,NANAL)**2+CONT(I,2*K,NANAL)**2)
          IF(NZERO(K).EQ.NSTEP(NANAL)) THEN
            REPLACED(K)=REPLACED(K).OR.CONT(I,K*3,NANAL).NE.ETOT(K)
            CONT(I,K*3,NANAL)=ETOT(K)
          ENDIF
          IF(ABS(CONT(I,K*3,NANAL)).GT.ABS(ETOT(K))) NLARGE=NLARGE+1
          IF(ABS(CONT(I,K*3,NANAL)).LT.ABS(ETOT(K))) NSMALL=NSMALL+1
        ENDDO
      ENDDO
      IF(WARNREPL) THEN
        DO K=-1,1,2
          IF(REPLACED(K)) THEN
            WRITE(CH16,'(A,SP,I2)',IOSTAT=DUMMY) 
     &       CHNAM(KCONT(3,NANAL))(:LENOCC(CHNAM(KCONT(3,NANAL)))),K
            CALL COMBOS_ERROR(1,'No value provided for '//
     &       CH16(:LENOCC(CH16)-1),'will compute total '//
     &       'uncertainty from stat. and total syst. uncertaintites')
          ENDIF
        ENDDO
      ENDIF
      IF(NSMALL.GT.0) CALL COMBOS_ERROR(1,'Total error '//
     & 'is sometimes smaller than sum of stat. and syst. contributions',
     & ' ')
      IF(NLARGE.GT.0) CALL COMBOS_ERROR(1,'Total systematic error '//
     &  'is sometimes larger than sum of stat. and syst. contributions',
     & ' ')
      ENDIF
*
*     Create links from the parameters to the systematic uncertainties
*
      NEWPAR=0
      DO IPAR=1,NPAR(NANAL)
        PCONT(IPAR,NANAL)=0
        DO ICONT=FCONT,NCONT(NANAL)
          IF(KCONT(ICONT,NANAL).EQ.KPAR(IPAR,NANAL))
     &     PCONT(IPAR,NANAL)=ICONT
        ENDDO
        DO ICALL=1,MCALLP
          IF(CHCALP(ICALL).EQ.CHNAM(KPAR(IPAR,NANAL)))
     &     PCONT(IPAR,NANAL)=-ICALL
        ENDDO
        IF(PCONT(IPAR,NANAL).EQ.0) THEN 
          CALL COMBOS_ERROR(1,'Parameter '//
     &     CHNAM(KPAR(IPAR,NANAL))(:LENOCC(CHNAM(KPAR(IPAR,NANAL))))//
     &     ' is not associated with a systematic contribution'//
     &     ' or preparation routine','this parameter ignored')
        ELSE
          NEWPAR=NEWPAR+1
          KPAR(NEWPAR,NANAL)=KPAR(IPAR,NANAL)
          DO K=-1,1
            PAR(K,NEWPAR,NANAL)=PAR(K,IPAR,NANAL)
          ENDDO
          PCONT(NEWPAR,NANAL)=PCONT(IPAR,NANAL)
        ENDIF
      ENDDO
      NPAR(NANAL)=NEWPAR
*
*     Check that this analysis supersedes if required by DISPLAY logical line
*
      IF(NDISP(NANAL).GE.0.AND.KDISP(0,NANAL).EQ.-1) THEN
        NSUP=0
        DO I=1,NLINK(NANAL)
          IF(RHOSTA(I,NANAL).GT.1.) NSUP=NSUP+1
        ENDDO
        IF(NSUP.EQ.0) THEN
          CALL COMBOS_ERROR(1,'Analysis does not supersede anything',
     &                        'DISPLAY logical line ignored')
          NDISP(NANAL)=-1
        ENDIF
      ENDIF
*
*     If display values are missing, use default values
*
      IF(NDISP(NANAL).GE.0.AND.KDISP(0,NANAL).EQ.-1) THEN
        JDISP=0
      ELSE
        JDISP=NDISP(NANAL)+1
      ENDIF
      DO IDISP=JDISP,MDISP
        IF(KDISP(IDISP,NANAL).GE.0)
     &   KDISP(IDISP,NANAL)=KCONT(IDISP,NANAL)
        DO I=1,NSTEP(NANAL)
          DISP(I,+IDISP,NANAL)=CONT(I,+IDISP,NANAL)
          DISP(I,-IDISP,NANAL)=CONT(I,-IDISP,NANAL)
        ENDDO
      ENDDO
*
*     Check if uncertainties are symmetric
*
      XSYMC(0,NANAL)=.TRUE.
      DO ICONT=1,NCONT(NANAL)
        XSYMC(ICONT,NANAL)=.TRUE.
        DO I=1,NSTEP(NANAL)
          XSYMC(ICONT,NANAL)=XSYMC(ICONT,NANAL).AND.
     &      CONT(I,ICONT,NANAL).EQ.-CONT(I,-ICONT,NANAL)
        ENDDO
      ENDDO
      DO IDISP=0,MDISP
        XSYMD(IDISP,NANAL)=.TRUE.
        IF(IDISP.GE.1) THEN
          DO I=1,NSTEP(NANAL)
            XSYMD(IDISP,NANAL)=XSYMD(IDISP,NANAL).AND.
     &        DISP(I,IDISP,NANAL).EQ.-DISP(I,-IDISP,NANAL)
          ENDDO
        ENDIF
      ENDDO
      END
*
************************************************************************
*
      SUBROUTINE CHECK_STEPS(NCOMB,ICOMB,IERR)
*     ========================================
*
*     Routine called by CHECK_ANAL and COMB_ANAL to checks steps of an 
*     individual analysis or a combined analysis.
*
*     Input:   NCOMB = number of analyses to be used in the combined analysis,
*     -----            or 1 in case of an individual analysis
*              ICOMB = array of dimension NCOMB, with the numbers of the 
*                      analyses to be used in the combined analysis,
*                      or the number of the individual analysis 
*
*     Output:  IERR  = error flag; if not zero, the analysis will be ignored
*     ------
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER IERR,NCOMB,ICOMB(NCOMB)
*
*     Externals
*  
      INTEGER LENOCC
      REAL VMIN,VMAX
*
*     Local variable
*
      INTEGER I,J,K,FOUND,NEWSTEP
      LOGICAL OK,CHOOSE
      CHARACTER*10 CH10
      CHARACTER*20 CH20
      CHARACTER*80 CH80
      LOGICAL REQUIRE_EACH
      DATA REQUIRE_EACH/.FALSE./
*
      INCLUDE 'combos.inc'
      INTEGER INDX(MSTEP)
      REAL TEMP(MSTEP),SMAX(MANAL),SMIN(MANAL),SMAXIM,SMINIM
*
      IERR=0
*
*     For an individual analysis, check that the number of steps is
*     at least 1; return if less than 2 steps; check that step name is present.
*
      IF(.NOT.XCOMB(0)) THEN
        IF(NSTEP(NANAL).LE.0) THEN
          CALL COMBOS_ERROR(1,'No valid DATA provided',
     &                        'this analysis will be ignored')
          IERR=IERR+1
        ENDIF
        IF(NSTEP(NANAL).LE.1) RETURN
        IF(KSTEP(NANAL).LE.0) THEN 
          CALL COMBOS_ERROR(1,'No STEP keyword provided',
     &                        'this analysis will be ignored')
          IERR=IERR+1
        ENDIF
      ENDIF
*
*     For a combined analysis, choose step values if none were specified
*
      CHOOSE=NSTEP(NANAL).LE.0.AND.XCOMB(0)
      IF(CHOOSE) THEN
        NSTEP(NANAL)=0
        DO I=1,NCOMB
          DO J=1,NSTEP(ICOMB(I))
            FOUND=0
            DO K=1,NSTEP(NANAL)
              IF(STEP(K,NANAL).EQ.STEP(J,ICOMB(I))) FOUND=K
            ENDDO
            IF(FOUND.EQ.0) THEN
              IF(NSTEP(NANAL).GE.MSTEP) CALL COMBOS_ERROR(-1,
     &         'Too many steps','please increase parameter MSTEP')
              NSTEP(NANAL)=NSTEP(NANAL)+1
              STEP(NSTEP(NANAL),NANAL)=STEP(J,ICOMB(I))
            ENDIF
          ENDDO
        ENDDO
      ENDIF
*
*     Sort steps in increasing step values
*   
      OK=.TRUE.
      DO I=2,NSTEP(NANAL)
        OK=OK.AND.STEP(I,NANAL).GE.STEP(I-1,NANAL)
      ENDDO
      IF(.NOT.OK) THEN 
        IF(.NOT.CHOOSE) CALL COMBOS_ERROR(0,
     &   'Step values not in increasing order','steps reordered')
        CALL SORTZV(STEP(1,NANAL),INDX,NSTEP(NANAL),1,0,0)
        CALL UCOPY(STEP(1,NANAL),TEMP,NSTEP(NANAL))
        DO I=1,NSTEP(NANAL)
          STEP(I,NANAL)=TEMP(INDX(I))
        ENDDO
        DO K=-NCONT(NANAL),NCONT(NANAL)
          CALL UCOPY(CONT(1,K,NANAL),TEMP,NSTEP(NANAL))
          DO I=1,NSTEP(NANAL)
            CONT(I,K,NANAL)=TEMP(INDX(I))
          ENDDO
crvk
          CALL UCOPY(CONT(1,SIGN(1,K)*MCONT+K,NANAL),TEMP,NSTEP(NANAL))
          DO I=1,NSTEP(NANAL)
            CONT(I,SIGN(1,K)*MCONT+K,NANAL)=TEMP(INDX(I))
          ENDDO
        ENDDO
      ENDIF
*
*     Check that step values are not all zero (only for individual analysis)
*
      IF(.NOT.XCOMB(0).AND.STEP(1,NANAL).EQ.0.AND.
     &                     STEP(NSTEP(NANAL),NANAL).EQ.0.) THEN 
        CALL COMBOS_ERROR(1,'No step values specified',
     &                      'this analysis will be ignored')
        IERR=IERR+1
        RETURN
      ENDIF
*
*     Check that same step value was not repeated; in case of duplicates:
*       - ignore individual analysis
*       - remove duplicates in combined analysis
*
      IF(XCOMB(0)) THEN 
        CH20='this step is ignored'
      ELSE
        CH20=' '
      ENDIF
      NEWSTEP=1
      TEMP(NEWSTEP)=STEP(1,NANAL)
      DO I=2,NSTEP(NANAL)
        IF(STEP(I,NANAL).NE.STEP(I-1,NANAL)) THEN
          NEWSTEP=NEWSTEP+1
          TEMP(NEWSTEP)=STEP(I,NANAL)
        ELSE IF(.NOT.CHOOSE) THEN
          WRITE(CH10,'(F10.4)') STEP(I,NANAL)
          CALL COMBOS_ERROR(0,'Step value '//
     &     CHNAM(KSTEP(NANAL))(:LENOCC(CHNAM(KSTEP(NANAL))))//' ='//
     &     CH10//' specified more than once',CH20)
        ENDIF
      ENDDO
      IF(XCOMB(0)) THEN
        NSTEP(NANAL)=NEWSTEP
        CALL UCOPY(TEMP,STEP(1,NANAL),NSTEP(NANAL))
      ELSE IF(NEWSTEP.NE.NSTEP(NANAL)) THEN
        CALL COMBOS_ERROR(1,'Step values are not all different',
     &                      'this analysis will be ignored')
        IERR=IERR+1
      ENDIF
      IF(.NOT.XCOMB(0)) RETURN
*
*     Determine step ranges
*
      DO I=1,NCOMB
        SMAX(I)=VMAX(STEP(1,ICOMB(I)),NSTEP(ICOMB(I)))
        SMIN(I)=VMIN(STEP(1,ICOMB(I)),NSTEP(ICOMB(I)))
      ENDDO
      SMAXIM=VMAX(SMAX,NCOMB)
      SMINIM=VMIN(SMIN,NCOMB)
*
*     Remove steps that are outside of all ranges
*
      NEWSTEP=0
      DO I=1,NSTEP(NANAL)
        IF(STEP(I,NANAL).LT.SMINIM.OR.STEP(I,NANAL).GT.SMAXIM) THEN
          WRITE(CH10,'(F10.4)') STEP(I,NANAL)
          CALL COMBOS_ERROR(0,'Cannot extrapolate all measurements',
     &                         'no combination will be performed at '//
     &        CHNAM(KSTEP(NANAL))(:LENOCC(CHNAM(KSTEP(NANAL))))//' ='//
     &        CH10(:LENOCC(CH10)))
        ELSE
          NEWSTEP=NEWSTEP+1
          STEP(NEWSTEP,NANAL)=STEP(I,NANAL)
        ENDIF
      ENDDO
      NSTEP(NANAL)=NEWSTEP
*
*     Remove step values that are outside at least one of the ranges
*
      IF(REQUIRE_EACH) THEN
        NEWSTEP=0
        DO I=1,NSTEP(NANAL)
          FOUND=0
          DO J=1,NCOMB
            IF(STEP(I,NANAL).LT.SMIN(J).OR.STEP(I,NANAL).GT.SMAX(J))THEN
              FOUND=FOUND+1
              IF(FOUND.EQ.1.AND..NOT.CHOOSE) THEN
                CALL ANALYSIS_NAME(.TRUE.,ICOMB(J),CH80)
                WRITE(CH10,'(F10.4)') STEP(I,NANAL)
                CALL COMBOS_ERROR(0,
     &  'Cannot extrapolate measurements in '//CH80(:LENOCC(CH80)),
     &  'no combination will be performed at '//
     &  CHNAM(KSTEP(NANAL))(:LENOCC(CHNAM(KSTEP(NANAL))))//' ='//
     &  CH10(:LENOCC(CH10)))
              ENDIF
            ENDIF
          ENDDO
          IF(FOUND.EQ.0) THEN 
            NEWSTEP=NEWSTEP+1
            STEP(NEWSTEP,NANAL)=STEP(I,NANAL)
          ENDIF
        ENDDO
        NSTEP(NANAL)=NEWSTEP
      ENDIF
*
*     Store desired step values for combination in array DEST
*
      IF(XCOMB(0).AND..NOT.CHOOSE) THEN
        NDEST=NSTEP(NANAL)
        CALL UCOPY(STEP(1,NANAL),DEST,NDEST)
      ENDIF
      END
*
************************************************************************
*
      SUBROUTINE COMB_ANAL(NWORDS,WORDS,IERR)
*     =======================================
*
*     Routine called for the "END" keyword in case of an "ouput" analysis 
*     (i.e. an analysis defined with COMB keywords, in other words a
*     combined analysis): 
*     various checks are preformed on the set of analyses to combine 
*     and the combination is done
*
*     DR Version 3.0 modif to allow several MEASUREMENT for combination
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NWORDS,IERR
      CHARACTER*(*) WORDS(*)
*
*     Externals
*  
      INTEGER LENOCC
*
*     Local variable
*
      INTEGER IANAL,I,J,K,ICONT,JCONT,ICO,JCO,JANAL,IPAR,JPAR
*DR 
      INTEGER NQUANTMAX,IQUANT
      LOGICAL OKQUANT
*DR END
      LOGICAL EXIST
      CHARACTER*80 CH80
*
      INCLUDE 'combos.inc'
      LOGICAL SUPER(MANAL)
*
      INTEGER NCOMB ! number of analysis to combine
      INTEGER ICOMB(MANAL) ! list of analyses to combine
*             ICOMB(I) = analysis number of the Ith analysis in the combination,
*                        I=1,NCOMB
*
*     Remove from combination analyses that are superseded by other ones
*     (note: recursion is allowed, i.e. a superseded analysis can supersede 
*            other analyses)
*
      DO JANAL=1,NANAL-1
        SUPER(JANAL)=.FALSE.
      ENDDO
      DO IANAL=1,NANAL-1
        IF(XCOMB(IANAL)) THEN
          DO I=1,NLINK(IANAL)
            IF(ABS(RHOSTA(I,IANAL)).GT.1.) THEN ! IANAL supersedes LINK(I,IANAL)
              DO JANAL=1,NANAL-1
                SUPER(JANAL)=SUPER(JANAL).OR.
     &                       LINK(0,JANAL).EQ.LINK(I,IANAL)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      DO IANAL=1,NANAL-1
        XCOMB(IANAL)=XCOMB(IANAL).AND..NOT.SUPER(IANAL)
      ENDDO
*
*     Count number of analyses to combine and
*     check that all analysis have same measured quantity and step variable
*
      NCOMB=0
*DR   
    1 DO IANAL=1,NANAL-1
        IF(XCOMB(IANAL)) THEN 
          IF(KCONT(0,NANAL).EQ.0) THEN
            CALL ANALYSIS_NAME(.TRUE.,IANAL,CH80)
            CALL COMBOS_ERROR(1,'No MEASUREMENT keyword',
     &      'measurement name inherited from '//CH80(:LENOCC(CH80)))
            CALL STORE_MEAS(1,CHNAM(KSYN(KCONT(0,IANAL))))
            IF(KCONT(0,NANAL).EQ.0) THEN
              CALL COMBOS_ERROR(1,'No MEASUREMENT keyword',
     &                            'this combination not performed')
              IERR=1
              RETURN
            ENDIF
          ENDIF
          IF(KSTEP(NANAL).EQ.-1.AND.KSYN(KSTEP(IANAL)).NE.0) THEN
            KSTEP(NANAL)=KSYN(KSTEP(IANAL))
            IF(KSTEP(NANAL).NE.0) THEN
              CALL ANALYSIS_NAME(.TRUE.,IANAL,CH80)
              CALL COMBOS_ERROR(1,'No STEP keyword',
     &        'step variable inherited from '//CH80(:LENOCC(CH80)))
            ENDIF
          ENDIF
*DR properly count number of quantities to be fit
*DR          IF(KSYN(KSTEP(IANAL)).EQ.KSYN(KSTEP(NANAL)).AND.
*DR     &       KSYN(KCONT(0,IANAL)).EQ.KSYN(KCONT(0,NANAL))) THEN
*DR            NCOMB=NCOMB+1
*DR            ICOMB(NCOMB)=IANAL
*DR          ENDIF
          IF(KSYN(KSTEP(IANAL)).EQ.KSYN(KSTEP(NANAL))) THEN
            NQUANTMAX=NQUANT
            OKQUANT=.FALSE.
            IF (.NOT.XCOMB(0)) THEN
              NQUANTMAX=MIN(NQUANTMAX,1)
            ENDIF
            DO IQUANT=1,NQUANTMAX
              IF (KSYN(KCONT(0,IANAL)).EQ.
     &            KSYN(KCONT(1-IQUANT,NANAL))) OKQUANT=.TRUE.
            ENDDO
            IF (OKQUANT) THEN
              NCOMB=NCOMB+1
              ICOMB(NCOMB)=IANAL
            ENDIF
          ENDIF
*DR end
        ENDIF
      ENDDO
      IF(NCOMB.EQ.0) THEN
        IF(KSTEP(NANAL).EQ.0) THEN 
          KSTEP(NANAL)=-1
          GOTO 1
        ENDIF
        CALL COMBOS_ERROR(1,'No analysis to combine',
     &                      'this combination not performed')
        IERR=2
        RETURN
      ELSE IF(NCOMB.EQ.1) THEN
*OS        CALL ANALYSIS_NAME(.TRUE.,ICOMB(1),CH80)
*OS        CALL COMBOS_ERROR(0,'Only one single analysis to combine: '//
*OS     &                      CH80(:LENOCC(CH80)),
*OS     &                      ' ')
      ENDIF
*
*     Check steps
*
      CALL CHECK_STEPS(NCOMB,ICOMB,IERR)
      IF(IERR.NE.0) RETURN
*
*     Don't allow more than one step in case NQUANT > 1
*
      IF(NSTEP(NANAL).GT.1.AND.NQUANT.GT.1) THEN
        CALL COMBOS_ERROR(1,
     &'More than one quantity to average and more than one step',
     &                      'this combination not performed')
        IERR=3
        RETURN
      ENDIF
*
*     Inherit parameters from the first analysis
*
      DO I=1,NCOMB
        IANAL=ICOMB(I)
        DO J=1,NPAR(IANAL)
          EXIST=.FALSE.
          DO K=1,NPAR(NANAL)
            IF(KPAR(J,IANAL).EQ.KPAR(K,NANAL)) EXIST=.TRUE.
          ENDDO
          IF(.NOT.EXIST.AND.PCONT(J,IANAL).GT.0) THEN
            IF(NPAR(NANAL).EQ.MPAR) CALL COMBOS_ERROR(-1,
     &       'Too many parameters','please increase parameter MPAR')
            NPAR(NANAL)=NPAR(NANAL)+1
            KPAR(NPAR(NANAL),NANAL)=KPAR(J,IANAL)
            DO K=-1,1
              PAR(K,NPAR(NANAL),NANAL)=PAR(K,J,IANAL)
            ENDDO
            CALL ANALYSIS_NAME(.TRUE.,IANAL,CH80)
            CALL COMBOS_ERROR(0,'Missing parameter '//
     &       CHNAM(KPAR(J,IANAL))(:LENOCC(CHNAM(KPAR(J,IANAL)))),
     &       'parameter inherited from '//CH80(:LENOCC(CH80)))
          ENDIF
        ENDDO
      ENDDO
*
*     Check that two different parameters have not been "synonymed";
*     if a parameter has synonyms, then main name should be parameter name
*
      DO I=1,NPAR(NANAL)
        IPAR=KPAR(I,NANAL)
        DO J=1,NPAR(NANAL)
          JPAR=KPAR(J,NANAL)
          IF(KSYN(IPAR).EQ.KSYN(JPAR).AND.I.NE.J) THEN
            CALL COMBOS_ERROR(1,'Parameters '//
     &         CHNAM(IPAR)(:LENOCC(CHNAM(IPAR)))//' and '//
     &         CHNAM(JPAR)(:LENOCC(CHNAM(JPAR)))//
     &         ' cannot be synonymed',
     &         'this combination not performed')
            IERR=4
            RETURN
          ENDIF
        ENDDO
        K=KSYN(IPAR)
        DO J=1,NNAM
          IF(KSYN(J).EQ.K) KSYN(J)=IPAR
        ENDDO
        IF(KSYN(IPAR).NE.IPAR) THEN
          PRINT * ,'Fatal problem in COMB_ANAL:',IPAR,KSYN(IPAR)
          STOP 7766
        ENDIF
      ENDDO
*
*     Check that, in each analysis to be combined, two different systematic 
*     contributions have not been synonymed
*
      DO I=1,NCOMB
        IANAL=ICOMB(I)
        DO ICO=0,NCONT(IANAL)
          ICONT=KCONT(ICO,IANAL)
          DO JCO=0,NCONT(IANAL)
            JCONT=KCONT(JCO,IANAL)
            IF(ICO.NE.JCO.AND.KSYN(ICONT).EQ.KSYN(JCONT).AND.
     &         ICONT.NE.0.AND.JCONT.NE.0) THEN
              CALL ANALYSIS_NAME(.TRUE.,IANAL,CH80)
              CALL COMBOS_ERROR(1,'Contributions '//
     &           CHNAM(ICONT)(:LENOCC(CHNAM(ICONT)))//' and '//
     &           CHNAM(JCONT)(:LENOCC(CHNAM(JCONT)))//' from '//
     &           CH80(:LENOCC(CH80))//' cannot be synonymed',
     &           'this combination not performed')
              IERR=5
              RETURN
            ENDIF
          ENDDO
        ENDDO
      ENDDO
*
*     Choose a default routine to be called if none was specified
*
      IF(NCALLC.EQ.0) THEN 
        NCALLC=1
        KCALLC(1)=1
        CALL COMBOS_ERROR(0,
     &   'No combination routine specified after CALL keyword',
     &   'routine '//CHCALC(1)(:LENOCC(CHCALC(1)))//
     &   ' will be called by default')
      ENDIF
*
*     Let now the master routine handle the combination
*
      CALL MASTER(NCOMB,ICOMB,NSTEP(NANAL),STEP(1,NANAL),
     &            CONT(1,-(FCONT-1),NANAL),IERR)
*
*     Check if uncertainties are symmetric
*
      XSYMC(0,NANAL)=.TRUE.
      DO ICONT=1,NCONT(NANAL)
        XSYMC(ICONT,NANAL)=.TRUE.
        DO I=1,NSTEP(NANAL)
          XSYMC(ICONT,NANAL)=XSYMC(ICONT,NANAL).AND.
     &      CONT(I,ICONT,NANAL).EQ.-CONT(I,-ICONT,NANAL)
        ENDDO
      ENDDO
      END
*
************************************************************************
*
      SUBROUTINE CHOOSE_NAME(NWORDS,NAME,WHAT,LETTER,NUMBER)
*     ======================================================
*
      IMPLICIT NONE
*
      INTEGER LENOCC
*
*     Arguments (all input except NAME)
*
      INTEGER NWORDS
      CHARACTER*(*) NAME
      CHARACTER*(*) WHAT,LETTER
      CHARACTER*256 XNAME
      CHARACTER*256 XWHAT
      INTEGER NUMBER
*
*     Local variable
*
      INTEGER STATUS
*
      IF(NWORDS.LT.1.OR.NAME.EQ.' ') THEN
        IF(NUMBER.GT.0) THEN 
          WRITE(NAME,'(A,I3.3)',IOSTAT=STATUS) LETTER,NUMBER
        ELSE
          NAME=LETTER
        ENDIF
        XWHAT = WHAT
        XNAME = NAME
        CALL COMBOS_ERROR(0,'No '//XWHAT(:LENOCC(WHAT))
     &  //' name specified',
     &  XWHAT(:LENOCC(WHAT))//' name set to '//XNAME(:LENOCC(NAME)))
      ENDIF

      XWHAT = WHAT
      XNAME = NAME
      IF(NWORDS.GT.1) CALL COMBOS_ERROR(0,
     &                'More than one '//XWHAT(:LENOCC(WHAT))
     &                //' name specified',
     &                XWHAT(:LENOCC(WHAT))//' name set to '
     &                //XNAME(:LENOCC(NAME)))
      END
*
************************************************************************
*
      SUBROUTINE CHECK_NAME(NAME,IGNORE_SIGN,INDX,LSIGN)
*     ==================================================
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER INDX,LSIGN
      CHARACTER*(*) NAME
      LOGICAL IGNORE_SIGN
*
*     Local variables
*
      INTEGER NNAMOLD
*
      INCLUDE 'combos.inc'
*
      NNAMOLD=NNAM
      CALL CHECK_CHAR(NAME,MNAM,NNAM,CHNAM,' ','MNAM',
     &                IGNORE_SIGN,INDX,LSIGN)
      IF(NNAM.EQ.NNAMOLD+1) KSYN(NNAM)=NNAM
      END
*
************************************************************************
*
      SUBROUTINE CHECK_CHAR(NAME,MNAMES,NNAMES,NAMES,WHAT,CHPAR,
     &                      IGNORE_SIGN,INDX,LSIGN)
*     ===============================================================
*
      IMPLICIT NONE
*
      INTEGER MNAMES,NNAMES,INDX,LSIGN
      CHARACTER*(*) NAME,WHAT,CHPAR
      CHARACTER*(*) NAMES(MNAMES)
      CHARACTER*256 XWHAT,XNAME,XCHPAR
      LOGICAL IGNORE_SIGN
*
*     Externals
*
***      INTEGER ICNTH
      INTEGER LENOCC
*
*     Local variable
*
      INTEGER I,L
*
*     Look for sign at end of name
*
*DR look for a % at the end of name, which would indicate that 
* error is given in %
* store in standalone common
*DR awful solution not to change all calling sequences
      LOGICAL XPERCENT,XAMPER,XDOLLAR
      COMMON /PERCENT/XPERCENT,XAMPER,XDOLLAR
      LOGICAL QUITLOOP

      XPERCENT=.FALSE.
      XAMPER=.FALSE.
      XDOLLAR=.FALSE.
      QUITLOOP=.FALSE. 
      L=0
      LSIGN=0
*DR had to rewrite the loop a bit
      DO I=LEN(NAME),1,-1
        IF(NAME(I:I).NE.' ') THEN 
          IF(NAME(I:I).EQ.'+') THEN 
            L=I-1
            LSIGN=+1 
*DR
            QUITLOOP=.TRUE.
          ELSE IF(NAME(I:I).EQ.'-') THEN 
            L=I-1
            LSIGN=-1
*DR
            QUITLOOP=.TRUE.
          ELSE IF(NAME(I:I).EQ.'%') THEN 
            XPERCENT=.TRUE.
*DR end
          ELSE IF(NAME(I:I).EQ.'&') THEN 
            XAMPER=.TRUE.
          ELSE IF(NAME(I:I).EQ.'$') THEN 
            XDOLLAR=.TRUE.
          ELSE
            L=I
            LSIGN=0
            QUITLOOP=.TRUE.
          ENDIF
*DR
          IF (QUITLOOP) GOTO 1
        ENDIF
      ENDDO
    1 IF(L.LE.0) L=MAX0(LENOCC(NAME),1)
      IF(L.GE.LEN(NAMES(1))) THEN 
        L=L-1
        XWHAT = WHAT
        XNAME = NAME
        CALL COMBOS_ERROR(1,XWHAT(:LENOCC(WHAT))//' name too long',
     &                      XWHAT(:LENOCC(WHAT))//' name truncated to '
     &                      //XNAME(:L))
      ENDIF
*
*     Is that name already used ?
*
***      INDX=ICNTH(NAME(:L),NAMES,MIN0(NNAMES,MNAMES)) 
***      Note that ICNTH does not work with names containing blanks !
      INDX=0
      DO I=1,MIN0(NNAMES,MNAMES)
        IF(NAME(:L).EQ.NAMES(I)) INDX=I
      ENDDO
      IF(INDX.EQ.0) THEN
        XWHAT = WHAT
        XCHPAR = CHPAR
        IF(NNAMES.EQ.MNAMES) CALL COMBOS_ERROR(-1,
     &  'Too many '//XWHAT(:LENOCC(WHAT))//' names',
     &  'please increase parameter '//XCHPAR(:LENOCC(CHPAR)))
        NNAMES=NNAMES+1
        INDX=NNAMES
        NAMES(NNAMES)=NAME(:L)
      ENDIF
*
*     Ignore sign ?
*
      IF(IGNORE_SIGN.AND.LSIGN.NE.0) THEN
        XNAME = NAME
        CALL COMBOS_ERROR(0,
     & 'Invalid name '//XNAME(:L+1),'name is taken as '//XNAME(:L))
      ENDIF

      END






