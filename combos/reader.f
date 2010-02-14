************************************************************************
*
*     "Reader" routines for COMBOS program
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 1.00, December  6, 1996:
*     -- original release
*     Version 1.10, December 10, 1996:
*     -- added 4 keywords: SUPERSEDE, STAT_CORR_WITH, NOSYSTEMATICS, SYMMETRIZE
*     -- changed characteristics of NOCORRELATIONS keyword
*     Version 1.20, December 13, 1996:
*     -- allow special instructions EXIT and QUIT in input files
*     Version 1.30, January 27, 1997:
*     -- added new keyword: SYNONYMS
*     Version 1.40, February 6, 1997:
*     -- added new keyword: LUMP
*     -- changed logic to set the "state" of the reader:
*          an "output state" now corresponds to each possible "input state"
*     Version 1.50, February 12, 1997:
*     -- new format of STEP logical line: numerical values allowed
*     Version 2.00, February 20, 1997:
*     -- new format for MEASUREMENTS and CALL logical lines
*     -- added new keyword: SWITCH
*     -- removed 4 obsolete keywords: NOADJUSTMENTS, NOCORRELATIONS, 
*                                     NOSYSTEMATICS, SYMMETRIZE
*     Version 2.30, May 12, 1997:
*     -- add the possibility to change the comment characters 
*        with the COMMENTS command
*     Version 2.31, Jun 19, 1997:
*     -- treat tab characters as blanks (i.e. spaces)
*     Version 2.32, Jul 16, 1997:
*     -- print error message if input string truncated
*     Version 2.40, October 14, 1997:
*     -- add logic for minos error option (HCJS).
*     Version 2.50, March 9, 1998:
*     -- added new keyword: DISPLAY
*
************************************************************************
*
      SUBROUTINE GET_NEXT_KEY(KEY,MWORDS,NWORDS,WORDS)
*     ================================================
*
*     Get next keyword from the input files
*
*     Input:   MWORDS = dimension of array WORDS
*     ------   
*
*     Output:  KEY    = string containing the next keyword
*     -------           without leading blanks and converted to upper case
*              NWORDS = nmber of words following the keyword
*              WORDS  = words following the keyword
*
*     Note:  When no more input is available, KEY is returned as ' '
*     -----  and NWORDS is returned as 0.
*
      IMPLICIT NONE
*
*     Arguments
*
      CHARACTER*(*) KEY
      CHARACTER*256 XKEY
      INTEGER MWORDS,NWORDS
      CHARACTER*(*) WORDS(MWORDS)
*
*     Externals
*
      INTEGER ICNTH,LENOCC
*
*     Local variables
*
      INTEGER IPOS,IKEY,INDX
      LOGICAL XVALUE
      REAL VALUE
      CHARACTER*16 WORD
*
      INTEGER NKEYWORDS,I
      PARAMETER(NKEYWORDS=16)
      CHARACTER*16 KEYWORDS(NKEYWORDS)
      INTEGER NWEXP(NKEYWORDS)
      LOGICAL NUMER(NKEYWORDS)
      CHARACTER*8 STATES(NKEYWORDS)
      CHARACTER*8 SET_STATE(NKEYWORDS)
      DATA (KEYWORDS(I),STATES(I),SET_STATE(I),NWEXP(I),NUMER(I),
     &     I=1,NKEYWORDS)/
*
*           KEYWORDS         STATES   SET_STATE    NWEXP  NUMER
     &     'OUTP*UT',        'ABCLlD','ABCLlD',       0, .TRUE.,
     &     'BEGI*N',         'A',     'B',           -1, .TRUE.,
     &     'COMB*INE',       'BCLl',  'LLLL',        -1, .FALSE.,
     &     'LUMP',           'L',     'L',           -2, .FALSE., 
     &     'MEAS*UREMENTS',  'BCLlD', 'BCllD',       -1, .FALSE.,
     &     'DISP*LAY',       'BD',    'BD',          -1, .FALSE.,
     &     'STEP*S',         'BCLlD', 'BCllD',        1, .TRUE.,
     &     'PARA*METERS',    'BCLlD', 'BCllD',        1, .TRUE.,
     &     'SWIT*CH',        'BCLl',  'CCCC',        -1, .TRUE.,
     &     'SYNO*NYMS      ','BCLl',  'CCCC',        -2, .FALSE., 
     &     'CALL',           'BCLl',  'CCCC',        -1, .TRUE.,
     &     'DATA',           'BD',    'DD',          -1, .TRUE., 
     &     'SUPE*RSEDE',     'BD',    'DD',          -1, .FALSE., 
     &     'STAT*_CORR_WITH','BD',    'DD',          -1, .TRUE., 
     &     'MINO*S',         'BCLl',  'CCCC',        -1, .FALSE.,
     &     'END',            'BCLlD', 'AAAAA',        0, .FALSE./
*
*     Note: a negative value of NWEXP(I) indicates that the number of
*           expected words is at least ABS(NWEXP(I))
*           
*
      INTEGER JKEY
      SAVE    JKEY
      DATA    JKEY/0/
      LOGICAL IGNORE
      SAVE    IGNORE
      DATA    IGNORE/.FALSE./
      CHARACTER*1 STATE
      SAVE STATE
      DATA STATE/'A'/
*
*     Default output arguments
*
      KEY=' '
      NWORDS=0
*
*     Get next word from input until next keyword
*
    1 CALL GET_NEXT_WORD(IPOS,KEY,XVALUE,VALUE)
      IF(KEY.EQ.' ') THEN 
        IF(STATE.NE.'A') THEN 
          CALL COMBOS_ERROR(1,
     &    'Input ends inside a BEGIN-END block',
     &    'this incompletely defined analysis ignored')
        ENDIF
        NWORDS=0
        RETURN
      ENDIF
      IKEY=0
      IF(IPOS.EQ.1) IKEY=ICNTH(KEY,KEYWORDS,NKEYWORDS)
      IF(IKEY.GT.0) IGNORE=.FALSE.
      IF(IGNORE) GOTO 1
      IF(IKEY.EQ.0) THEN 
        IF(JKEY.EQ.0) THEN 
          XKEY = KEY
          CALL COMBOS_ERROR(1,
     &      'Input does not start with recognized keyword: '//
     &      XKEY(:LENOCC(KEY)),'this word ignored')
          GOTO 1
        ENDIF
        CALL SET_SAME_WORD
        IKEY=JKEY
      ENDIF
      KEY=KEYWORDS(IKEY)(:4)
*
*     Check that keyword can be used now
*
      INDX=INDEX(STATES(IKEY),STATE)
      IF(INDX.LE.0) THEN
        XKEY = KEY
        CALL COMBOS_ERROR(1,'Keyword '//XKEY(:4)
     &                      //' cannot be used here',
     &                      'this logical line ignored')
        IGNORE=.TRUE.
        GOTO 1
      ENDIF
      IF(SET_STATE(IKEY)(INDX:INDX).NE.' ')
     & STATE=SET_STATE(IKEY)(INDX:INDX)
      JKEY=IKEY
*
*     Read and count the non-numerical words following the keyword
*
    2 CALL GET_NEXT_WORD(IPOS,WORD,XVALUE,VALUE)
      IF(WORD.EQ.' ') GOTO 3
      IF(XVALUE) THEN ! numerical value
        IF(NUMER(IKEY)) THEN
          CALL SET_SAME_WORD
        ELSE
          XKEY = KEY
          CALL COMBOS_ERROR(0,'Unexpected numerical data after '//
     &                        XKEY(:LENOCC(KEY))//' keyword',
     &                        'these numerical data ignored')
   21     CALL GET_NEXT_VALUE(1,I,VALUE)
          IF(I.GT.0) GOTO 21
        ENDIF
        GOTO 3
      ENDIF
      I=0
      IF(IPOS.EQ.1) I=ICNTH(WORD,KEYWORDS,NKEYWORDS)
      IF(I.NE.0) THEN ! following keyword
        CALL SET_SAME_WORD
        GOTO 3
      ENDIF
      IF(NWORDS.EQ.MWORDS) CALL COMBOS_ERROR(-1,
     &   'Too many consecutive names on logical line',
     &   'please increase parameter MWORDS')
      IF(NWORDS.EQ.NWEXP(IKEY)) THEN 
        XKEY = KEY
        CALL COMBOS_ERROR(0,'Unexpected string '//WORD(:LENOCC(WORD))//
     &                      ' after '//XKEY(:LENOCC(KEY))//' keyword',
     &                      'this string ignored')
      ELSE
        NWORDS=NWORDS+1
        WORDS(NWORDS)=WORD
      ENDIF
      GOTO 2
*
    3 IF(NWORDS.LT.IABS(NWEXP(IKEY))) THEN
        XKEY = KEY
        CALL COMBOS_ERROR(0,
     &  'Too few strings found after '//XKEY(:LENOCC(KEY))
     &  //' keyword',' ')
      ENDIF

      END
*
************************************************************************
*
      SUBROUTINE GET_NEXT_VALUE(MVAL,NVAL,VAL)
*     ========================================
*
*     Get array of next values from the input files
*
*     Input:   MVAL   = size of array VAL = number of desired values
*     ------
*
*     Output:  NVAL   = number of values returned
*     -------  VAL    = next values 
*
*     Note:  When no more input is available, NVAL and VAL are
*     -----  returned as 0.
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER MVAL,NVAL
      REAL VAL(MVAL)
*
*     Local variables
*
      INTEGER IPOS,I
      CHARACTER*16 WORD
      LOGICAL XVALUE
*
*     Default output arguments
*
      NVAL=0
      CALL UZERO(VAL,1,MVAL)
      DO I=1,MVAL
        CALL GET_NEXT_WORD(IPOS,WORD,XVALUE,VAL(I))
        IF(WORD.EQ.' ') GOTO 1
        IF(.NOT.XVALUE) THEN 
          CALL SET_SAME_WORD
          GOTO 1
        ENDIF
        NVAL=NVAL+1
      ENDDO
    1 IF(NVAL.GT.0.AND.NVAL.NE.MVAL) CALL COMBOS_ERROR(0,
     & 'Too few values specified','missing values set to 0')
      END
*
************************************************************************
*
      SUBROUTINE GET_NEXT_WORD(IPOS,NEXTWORD,XVALUE,VALUE)
*     ====================================================
*
*     Get next word from the input files
*
*     Output:  IPOS   = position of next word in physical input line
*     -------  WORD   = string containing the next word
*                       without leading blanks and converted to upper case
*              LOGX   = logical to indicate if word contains a value
*              VALUE  = associated value (or 0 if XVALUE=.FALSE.)
*
*     Note:  When no more input is available, WORD is returned as ' '
*     -----  and IPOS is returned as 0
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER IPOS
      CHARACTER*(*) NEXTWORD
      LOGICAL XVALUE
      REAL VALUE
*
*     Arguments of entries GET_WHERE and SET_INPUT_NAME
*
      INTEGER WDEPTH,WUNIT,WLINE,WWORD 
      CHARACTER*(*) WFILE
      INTEGER XLEN
*
*     Externals
*
      INTEGER LENOCC 
*
*     Local variables
*
      INTEGER LUNDAT,MDEPTH
      PARAMETER(LUNDAT=90,MDEPTH=9)
      INTEGER ODEPTH,OUNIT,OLINE,OWORD 
      SAVE    ODEPTH,OUNIT,OLINE,OWORD 
      CHARACTER*80 OFILE
      SAVE         OFILE
      INTEGER IDEPTH
      SAVE    IDEPTH
      DATA    IDEPTH/0/
      INTEGER LUNIT(0:MDEPTH)
      SAVE    LUNIT
      DATA    LUNIT/5,MDEPTH*0/
      CHARACTER*80 CFILE(0:MDEPTH) ! current file
      SAVE         CFILE
      DATA         CFILE/'standard input',MDEPTH*' '/
      INTEGER JWORD ! current word
      SAVE    JWORD
      DATA    JWORD/0/
      INTEGER JLINE(0:MDEPTH) ! current line
      SAVE    JLINE
      DATA    JLINE/0,MDEPTH*0/
      CHARACTER*256 LINE
      SAVE          LINE
      DATA          LINE/' '/
      CHARACTER*16  WORD
      SAVE          WORD
      DATA          WORD/' '/
      LOGICAL NEXT
      SAVE    NEXT
      DATA    NEXT/.TRUE./
      INTEGER STATUS,IND,I
      INTEGER NCOMCHARS
      PARAMETER(NCOMCHARS=4)
      CHARACTER*1 COMCHARS(NCOMCHARS)
      SAVE        COMCHARS
      DATA COMCHARS/'*','#','|','!'/ ! characters used for comments
      CHARACTER*(NCOMCHARS) CCHARS
      CHARACTER*8 CHINCL
      DATA CHINCL/'INCLUDE'/
      CHARACTER*8 CHCOMM
      DATA CHCOMM/'COMMENTS'/
      LOGICAL TRUNC
*
      IF(.NOT.NEXT) THEN 
        NEXT=.TRUE.
        GOTO 9
      ENDIF
      ODEPTH=IDEPTH
      OUNIT=LUNIT(IDEPTH)
      OFILE=CFILE(IDEPTH)
      OLINE=JLINE(IDEPTH)
      OWORD=JWORD
    1 IF(LINE.EQ.' ') THEN
*
*       Read next line
*
        JWORD=0
        JLINE(IDEPTH)=JLINE(IDEPTH)+1
        READ(LUNIT(IDEPTH),'(A)',ERR=80,END=90) LINE
        CALL GET_FIRST_WORD(LINE,WORD,TRUNC)
        CALL CLTOU(WORD)
        IF(WORD.EQ.
     &     CHINCL(:MIN0(MAX0(LENOCC(WORD),4),LEN(CHINCL)))) THEN
          IF(IDEPTH.EQ.MDEPTH) 
     &      CALL READER_ERROR(CFILE(IDEPTH),JLINE(IDEPTH),-1,
     &      'Too many INCLUDE nested','please rearrange input files')
          IF(LINE.EQ.' ') THEN
            CALL READER_ERROR(CFILE(IDEPTH),JLINE(IDEPTH),0,
     &      'No include file specified','this line ignored')
            GOTO 1
          ENDIF
          IDEPTH=IDEPTH+1
          JLINE(IDEPTH)=0
          LUNIT(IDEPTH)=LUNDAT+IDEPTH
          CALL GET_FIRST_WORD(LINE,CFILE(IDEPTH),TRUNC)
          IF(TRUNC) CALL COMBOS_ERROR(1,'word too long',
     &      'string truncated to '//WORD//'; unpredictable results !')
          LINE=' '
          OPEN(UNIT=LUNIT(IDEPTH),FILE=CFILE(IDEPTH),STATUS='OLD',
     &         FORM='FORMATTED',ACCESS='SEQUENTIAL',IOSTAT=STATUS)
          IF(STATUS.NE.0) THEN 
            CALL READER_ERROR(CFILE(IDEPTH-1),JLINE(IDEPTH-1),1,
     &     'Cannot open file '//CFILE(IDEPTH)(:LENOCC(CFILE(IDEPTH))),
     &     'this line ignored')
            IDEPTH=IDEPTH-1
          ENDIF
          GOTO 1
        ELSE IF(WORD.EQ.'EXIT') THEN
          GOTO 90
        ELSE IF(WORD.EQ.'QUIT') THEN
          IF(IDEPTH.GT.0) THEN
            CALL READER_ERROR(CFILE(IDEPTH),JLINE(IDEPTH),0,
     &      'QUIT instruction found','all input files closed')
            DO I=1,IDEPTH
              CLOSE(LUNIT(I))
            ENDDO
            IDEPTH=0
          ENDIF
          GOTO 90
        ELSE IF(WORD.EQ.
     &     CHCOMM(:MIN0(MAX0(LENOCC(WORD),4),LEN(CHCOMM)))) THEN
          CALL GET_FIRST_WORD(LINE,CCHARS,TRUNC)
          IF(TRUNC) CALL COMBOS_ERROR(1,'word too long',
     &      'string truncated to '//WORD//'; unpredictable results !')
          LINE=' '
          DO I=1,NCOMCHARS
            IF(CCHARS(I:I).NE.' ') THEN
              COMCHARS(I)=CCHARS(I:I)
            ELSE IF(I.GT.1) THEN
              COMCHARS(I)=COMCHARS(I-1)
            ENDIF
          ENDDO
          GOTO 1
        ENDIF
        LINE=WORD//LINE
*
*       Remove comments
*
        DO I=1,NCOMCHARS
          IND=INDEX(LINE,COMCHARS(I))
          IF(IND.EQ.1) LINE=' '
          IF(IND.GT.1.AND.I.NE.1) LINE=LINE(:IND-1)
        ENDDO
*
*       Convert line to upper case
*
        CALL CLTOU(LINE)
        GOTO 1
      ENDIF
*
*     Get next word in line
*
      JWORD=JWORD+1
      CALL GET_FIRST_WORD(LINE,WORD,TRUNC)
      IF(TRUNC) CALL COMBOS_ERROR(1,'word too long',
     &  'string truncated to '//WORD//'; unpredictable results !')
*
*     Return output arguments
*
    9 IPOS=JWORD
      NEXTWORD=WORD
      XVALUE=INDEX('0123456789+-.',WORD(:1)).GT.0
      IF(XVALUE) THEN 
        READ(WORD,'(F16.0)',IOSTAT=STATUS) VALUE
        XVALUE=STATUS.EQ.0
      ENDIF
      IF(.NOT.XVALUE) VALUE=0.
      RETURN
*
*     Handle read errors
* 
   80 CALL READER_ERROR(CFILE(IDEPTH),JLINE(IDEPTH),1,
     &                  'FORTRAN read error','this line ignored')
      LINE=' '
      GOTO 1
*
*     Handle end-of-file conditions
*
   90 CLOSE(LUNIT(IDEPTH))
      IF(IDEPTH.GT.0) THEN
        IDEPTH=IDEPTH-1
        LINE=' '
        GOTO 1
      ENDIF
      WORD=' '
      NEXT=.FALSE.
      GOTO 9
*
      ENTRY SET_SAME_WORD
*     ===================
*
      NEXT=.FALSE.
      RETURN
*
      ENTRY GET_WHERE(WDEPTH,WUNIT,WFILE,WLINE,WWORD)
*     ===============================================
*
      IF(NEXT) THEN
        WDEPTH=IDEPTH
        WUNIT=LUNIT(IDEPTH)
        WFILE=CFILE(IDEPTH)
        WLINE=JLINE(IDEPTH)
        WWORD=JWORD
      ELSE
        WDEPTH=ODEPTH
        WUNIT=OUNIT
        WFILE=OFILE
        WLINE=OLINE
        WWORD=OWORD
      ENDIF
      RETURN
*
      ENTRY SET_INPUT_NAME(XLEN,WFILE)
*     ===========================
*
      CFILE(0)=WFILE(1:XLEN)
      END
*
************************************************************************
*
      SUBROUTINE GET_FIRST_WORD(LINE,WORD,TRUNC)
*     ==========================================
*
*     Input:  LINE = line of text
*     ------  
*
*     Output: LINE = line without first word
*     ------- WORD = first word in line (without leading blanks)
*             TRUNC= .TRUE. the first word has been truncated
*
*     Note: tabs are treated as blanks
*
      IMPLICIT NONE
*
*     Arguments
*
      CHARACTER*(*) LINE,WORD
      LOGICAL TRUNC
*
*     Local variables
*
      INTEGER I
      INTEGER*8 IITAB
      CHARACTER*1 TAB
      EQUIVALENCE(IITAB,TAB)
      DATA IITAB/'9'X/ ! tab character, hexadecimal ASCII code = 9
*
      IF(LINE.EQ.' ') THEN 
        WORD=' '
        RETURN
      ENDIF
    1 IF(LINE(1:1).EQ.' '.OR.LINE(1:1).EQ.TAB) THEN
        LINE=LINE(2:)
        GOTO 1
      ENDIF
      DO I=1,LEN(LINE)
        IF(LINE(I:I).EQ.' '.OR.LINE(I:I).EQ.TAB) THEN
          WORD=LINE(:I)
          TRUNC=WORD.NE.LINE(:I)
          LINE=LINE(I:)
          RETURN
        ENDIF
      ENDDO
      WORD=LINE
      TRUNC=WORD.NE.LINE
      LINE=' '
      END
