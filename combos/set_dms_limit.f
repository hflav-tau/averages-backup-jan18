************************************************************************
*
*     "SET_DMS_LIMIT" termination routine for COMBOS program
*     "KUMAC_DMS_LIMIT" termination routine for COMBOS program
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 2.00, February 20, 1997:
*     -- original release
*     Version 2.10, February 27, 1997:
*     -- determine "sensitive regions" in addition to excluded regions
*     -- determine "sensitivity" in addition to limit
*     -- ability to write results in a KUMAC file
*     -- new entry KUMAC_DMS_LIMIT; with this new entry, routine SET_DMS_LIMIT
*        is now treated as two different termination routines, called 
*        SET_DMS_LIMIT and KUMAC_DMS_LIMIT (see routine MASTER for the special
*        way to call termination routine KUMAC_DMS_LIMIT)
*     Version 2.50, March 9, 1998:
*     -- print also asymmetric errors (in addition to symmetrized error)
*        in KUMAC file
*     Version 2.97, July 21, 2000:
*     -- parameter MU increased from 100 to MQUAN**2 (=12**2=144)
*     Version 2.98, September 7, 2000:
*     -- new entry SET_DMS_LIMIT_USER to receive user values 
*        for the sensitivities
*     Version 3.01, February 8, 2002
*     -- increase parameter MEXCL from 9 to 15 
*     Version 3.32, January 17, 2010
*     -- put variables UVAL and XVAL in a common block,
*        and replace entry SET_DMS_LIMIT_USER with a subroutine
*        (requested by A. Lusiani to compile with gfortran on SL5)
*
************************************************************************
*
      BLOCK DATA UVAL_CM_INIT
      REAL UVAL(3,2)
      LOGICAL XVAL
      COMMON/UVAL_CM/UVAL,XVAL
      DATA UVAL/6*-1./
      DATA XVAL/.FALSE./
      END
*
************************************************************************
*
      SUBROUTINE SET_DMS_LIMIT(NSTEP,NC,MSTEP,
     &                         DMS,AMP,ERR2P,ERR2N,ERR2,CL,IERR)
*     ==========================================================
*
*     Termination routine to compute dms limit
*
*     Input:  NSTEP  = number of steps
*     -----   NC     = number of cases
*                      NC=1: means only total error available
*                      NC=2: means also statistical error available
*             MSTEP  = first dimension of arrays DMS, AMP, ... 
*             DMS    = array of DMS values (i.e. step values)
*             AMP    = array of amplitude values
*             ERR2P  = array of positive errors**2 on AMP
*             ERR2N  = array of negative errors**2 on AMP
*             ERR2   = array of symmetric errors**2 on AMP
*             CL     = array of confidence levels of AMP
*             IERR   = array of error codes (0 means OK)
*       
*     Output:  none
*     -------
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER MC,NC,NSTEP,MSTEP
      PARAMETER(MC=2)
      DOUBLE PRECISION DMS(MSTEP),AMP(MSTEP,MC),CL(MSTEP,MC),
     &                 ERR2P(MSTEP,MC),ERR2N(MSTEP,MC),ERR2(MSTEP,MC)
      INTEGER IERR(0:MSTEP,MC)
*
*     Argument of entry KUMAC_DMS_LIMIT
*
      LOGICAL KUM
*
*     Argument of entry SET_DMS_LIMIT_USER
*
      REAL USER(3)
      LOGICAL XUSER
*  
      INCLUDE 'master.inc' ! input data of last combination performed
*
*     Externals
*
      INTEGER LENOCC
      DOUBLE PRECISION DGAUSN
*
*     Local variables
*
      INTEGER IC,LTOT,LSTA
      PARAMETER(LTOT=1,LSTA=2) ! cases
      CHARACTER*16 TEXT(MC)
      DATA TEXT(LTOT)/'(syst. included)'/
      DATA TEXT(LSTA)/'(stat. only)'/
      CHARACTER*4 SUFFIX(MC)
      DATA SUFFIX(LTOT)/'    '/
      DATA SUFFIX(LSTA)/'stat'/
      CHARACTER*12 CHLIM(2)
      DATA CHLIM/'lower limit','sensitivity'/
      INTEGER MNAMES
      CHARACTER*7 NAME
      PARAMETER(MNAMES=9)
      CHARACTER*3 NAMES(MNAMES)
      DATA NAMES/'amp','erp','ern','err','prb','irc','exc','lim','sen'/
      DOUBLE PRECISION CL_DEF,CONLEV,NSIGM
      PARAMETER(CL_DEF=0.95D0) ! 95% confidence level by default
      INTEGER IREF,IDIR,I,J,IFAIL
      INTEGER MU
      PARAMETER(MU=MQUAN**2)
      DOUBLE PRECISION U(MU),DAMPL
      INTEGER MEXCL,MZERO,NZERO,NZ
      PARAMETER(MEXCL=15,MZERO=2*MEXCL)
      INTEGER DERI(MZERO),SORT(MZERO),ITEMP(MZERO),NEXCL
      REAL LIM(2),ZERO(MZERO),RTEMP(MZERO),EXCL(2,MEXCL)
      REAL ULIM(2)
*
*     Saved local variable
*
      LOGICAL KUMAC
      SAVE    KUMAC
      DATA    KUMAC/.FALSE./
      REAL UVAL(3,2)
*      SAVE UVAL ! commented on Jan 17, 2010
*      DATA UVAL/6*-1./ ! commented on Jan 17, 2010
      LOGICAL XVAL
*      SAVE    XVAL ! commented on Jan 17, 2010
*      DATA    XVAL/.FALSE./ ! commented on Jan 17, 2010
      COMMON/UVAL_CM/UVAL,XVAL ! added on Jan 17, 2010
*
      IF(.NOT.KUMAC.AND.LUNIT.GT.0)
     & WRITE(LUNIT,2000) CHROUT(:LENOCC(CHROUT))
 2000 FORMAT(/,' SET_DMS_LIMIT called for results obtained with ',A,':')
*
*     Perform checks
*
      IFAIL=0
      IF(NC.LE.0.OR.NC.GT.MC) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'SET_DMS_LIMIT: invalid number of cases = ',NC
        IFAIL=IFAIL+1
      ENDIF
      IF(NSTEP.LE.1) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'SET_DMS_LIMIT: less than 2 steps provided'
        IFAIL=IFAIL+1
      ENDIF
      IF(NSTEP.GT.MU) THEN 
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'SET_DMS_LIMIT: parameter MU =',MU,' less than NSTEP = ',NSTEP
        IFAIL=IFAIL+1
      ENDIF
      IF(CHSTEP.NE.'DMS') THEN 
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'SET_DMS_LIMIT: step name is not DMS'
        IFAIL=IFAIL+1
      ENDIF
      IF(CHMEAS.NE.'A'.AND.CHMEAS.NE.'AMP'.AND.CHMEAS.NE.'AMPLITUDE')
     &THEN 
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'SET_DMS_LIMIT: amplitude name is not recognized'
        IFAIL=IFAIL+1
      ENDIF
      IF(IFAIL.NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &                 'SET_DMS_LIMIT: cannot do the job !'
        RETURN
      ENDIF
*
*     Get desired confidence level
*
      CONLEV=CL_DEF ! default confidence level
      DO I=1,NPARA
        IF(CHPARA(I).EQ.'CONFIDENCE') CONLEV=PARA(I)
      ENDDO
      NSIGM=DGAUSN(CONLEV)! number of sigmas corresponding to a given 1-sided CL
*
*     Write KUMAC file
*
      IF(KUMAC.AND.LUNIT.GT.0) THEN
        WRITE(LUNIT,6001) CHROUT(:LENOCC(CHROUT)),
     &                    CHCOMB(:LENOCC(CHCOMB))
        WRITE(LUNIT,6002) 'global/delete name_*'
        WRITE(LUNIT,6003) 'name_rout',CHROUT(:LENOCC(CHROUT))
        WRITE(LUNIT,6003) 'name_comb',CHCOMB(:LENOCC(CHCOMB))
        DO I=1,NMEAS
          WRITE(NAME,'(I7)') I
          CALL CLEFT(NAME,1,LEN(NAME))
          NAME='name_'//NAME
          WRITE(LUNIT,6003) NAME(:LENOCC(NAME)),
     &                      CHANAL(I)(:LENOCC(CHANAL(I)))
        ENDDO
        WRITE(LUNIT,6007) 'ncomb',FLOAT(NMEAS)
        WRITE(LUNIT,6004) 'cl',1,CONLEV
        WRITE(LUNIT,6004) 'dms',NSTEP,(DMS(I),I=1,NSTEP)
      ENDIF
*
*     Loop on cases (total error, then stat. only)
*    
      DO IC=1,NC
        IF(LUNIT.GT.0) THEN
          IF(KUMAC) THEN
            NAME=NAMES(1)//SUFFIX(IC)
            WRITE(LUNIT,6004) NAME(:LENOCC(NAME)),
     &                        NSTEP,(AMP(I,IC),I=1,NSTEP)
            NAME=NAMES(2)//SUFFIX(IC)
            WRITE(LUNIT,6004) NAME(:LENOCC(NAME)),
     &                        NSTEP,(+DSQRT(ERR2P(I,IC)),I=1,NSTEP)
            NAME=NAMES(3)//SUFFIX(IC)
            WRITE(LUNIT,6004) NAME(:LENOCC(NAME)),
     &                        NSTEP,(-DSQRT(ERR2N(I,IC)),I=1,NSTEP)
            NAME=NAMES(4)//SUFFIX(IC)
            WRITE(LUNIT,6004) NAME(:LENOCC(NAME)),
     &                        NSTEP,(DSQRT(ERR2(I,IC)),I=1,NSTEP)
            NAME=NAMES(5)//SUFFIX(IC)
            WRITE(LUNIT,6004) NAME(:LENOCC(NAME)),
     &                        NSTEP,(CL(I,IC),I=1,NSTEP)
            NAME=NAMES(6)//SUFFIX(IC)
            WRITE(LUNIT,6004) NAME(:LENOCC(NAME)),
     &                        NSTEP,(FLOAT(IERR(I,IC)),I=1,NSTEP)
          ELSE
            WRITE(LUNIT,1000)
     &       CHSTEP(:LENOCC(CHSTEP)),TEXT(IC)(:LENOCC(TEXT(IC)))
 1000       FORMAT(/,' Determine ',A,' limit ',A,':')
          ENDIF
        ENDIF
        CALL UZERO(LIM,1,2)
        DO IREF=0,2 ! IREF=2 for the sensitivity
          NZERO=0
          DO IDIR=-1,1,2
            DO I=1,NSTEP
              IF(IREF.EQ.2) THEN ! sensitivity
                DAMPL=-1.D0
              ELSE
                DAMPL=AMP(I,IC)-DBLE(IREF)
              ENDIF
              IF(IDIR.EQ.-1) THEN
                U(I)=DAMPL-NSIGM*DSQRT(ERR2N(I,IC))
              ELSE 
                U(I)=DAMPL+NSIGM*DSQRT(ERR2P(I,IC))
              ENDIF
            ENDDO
            CALL FIND_ZERO(NSTEP,DMS,U,MZERO-NZERO,NZ,ZERO(NZERO+1),
     &                                                DERI(NZERO+1))
            IF(NZ+NZERO+2.GT.MZERO) THEN 
              IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &         'SET_DMS_LIMIT: parameter MZERO too small ',
     &         MZERO,NZ+NZERO+2
              STOP 46
            ENDIF
            DO J=1,NZ
              DERI(NZERO+J)=IDIR*DERI(NZERO+J)
            ENDDO
            NZERO=NZERO+NZ
          ENDDO
          DO I=1,NSTEP,NSTEP-1
            NZERO=NZERO+1
            ZERO(NZERO)=SNGL(DMS(I))
            IF(IREF.EQ.2) THEN ! sensitivity
              DAMPL=-1.D0
            ELSE
              DAMPL=AMP(I,IC)-DBLE(IREF)
            ENDIF
            IF((DAMPL+NSIGM*DSQRT(ERR2P(I,IC)).LE.0.D0).OR.
     &         (DAMPL-NSIGM*DSQRT(ERR2N(I,IC)).GE.0.D0)) THEN
              DERI(NZERO)=-1
            ELSE
              DERI(NZERO)=1
            ENDIF
          ENDDO
          DERI(NZERO)=-DERI(NZERO)
          CALL SORTZV(ZERO,SORT,NZERO,1,0,0)
          CALL UCOPY(ZERO,RTEMP,NZERO)
          CALL UCOPY(DERI,ITEMP,NZERO)
          DO I=1,NZERO
            ZERO(I)=RTEMP(SORT(I))
            DERI(I)=ITEMP(SORT(I))
          ENDDO
          NEXCL=0
          DO I=1,NZERO-1
            IF(DERI(I).EQ.-1.AND.DERI(I+1).EQ.1) THEN 
              IF(LUNIT.GT.0.AND.IREF.LT.2) THEN
                IF(.NOT.KUMAC) THEN
                  WRITE(LUNIT,1001)
     &             CHMEAS(:LENOCC(CHMEAS)),IREF,CONLEV,
     &             CHSTEP(:LENOCC(CHSTEP)),ZERO(I),ZERO(I+1)
 1001             FORMAT(1X,A,' = ',I1,' excluded at ',2P,F3.0,0P,TL1,
     &            '% CL for ',A,' values from',F10.4,' to',F10.4)
                ELSE IF(IREF.EQ.1) THEN
                  NEXCL=NEXCL+1
                  EXCL(1,NEXCL)=ZERO(I)
                  EXCL(2,NEXCL)=ZERO(I+1)
                ENDIF
              ENDIF
              IF(LUNIT.GT.0.AND.IREF.EQ.2.AND..NOT.KUMAC)
     &         WRITE(LUNIT,1002) CHMEAS(:LENOCC(CHMEAS)),CONLEV,
     &                         CHSTEP(:LENOCC(CHSTEP)),ZERO(I),ZERO(I+1)
 1002         FORMAT(1X,A,' "sensitive"  at ',2P,F3.0,0P,TL1,
     &        '% CL for ',A,' values from',F10.4,' to',F10.4)
              IF(IREF.NE.0.AND.ZERO(I).EQ.SNGL(DMS(1)))
     &         LIM(IREF)=ZERO(I+1)
            ENDIF
          ENDDO
          IF(IREF.EQ.1.AND.KUMAC.AND.LUNIT.GT.0) THEN
            NAME=NAMES(7)//SUFFIX(IC)
            WRITE(LUNIT,6005) NAME(:LENOCC(NAME)),
     &                        NEXCL,((EXCL(J,I),J=1,2),I=1,NEXCL)
          ENDIF
        ENDDO
        IF(LUNIT.GT.0) THEN 
          IF(KUMAC) THEN
            NAME=NAMES(8)//SUFFIX(IC)
            IF(XVAL.AND.UVAL(IC,1).GT.0.) THEN
              WRITE(LUNIT,6014) NAME(:LENOCC(NAME)),1,LIM(1)
              WRITE(LUNIT,6024) NAME(:LENOCC(NAME)),1,UVAL(IC,1)
            ELSE
              WRITE(LUNIT,6004) NAME(:LENOCC(NAME)),1,LIM(1)
            ENDIF
            NAME=NAMES(9)//SUFFIX(IC)
            IF(XVAL.AND.UVAL(IC,2).GT.0.) THEN
              WRITE(LUNIT,6014) NAME(:LENOCC(NAME)),1,LIM(2)
              WRITE(LUNIT,6024) NAME(:LENOCC(NAME)),1,UVAL(IC,2)
            ELSE
              WRITE(LUNIT,6004) NAME(:LENOCC(NAME)),1,LIM(2)
            ENDIF
          ELSE
            WRITE(LUNIT,1003)
     &                 (CONLEV,CHLIM(I),CHSTEP(:LENOCC(CHSTEP)),
     &                  LIM(I),TEXT(IC)(:LENOCC(TEXT(IC))),I=1,2)
 1003       FORMAT(1X,2P,F3.0,0P,TL1,'% CL ',A,'on ',A,' =',F5.1,1X,A)
            IF(XVAL.AND.(UVAL(IC,1).GT.0..OR.UVAL(IC,2).GT.0.)) THEN
              WRITE(LUNIT,2003) 
 2003         FORMAT(' !!! Above results overwritten by user values:')
              DO I=1,2
                IF(UVAL(IC,I).GT.0.) THEN
                  ULIM(I)=UVAL(IC,I)
                ELSE
                  ULIM(I)=LIM(I)
                ENDIF
              ENDDO
              WRITE(LUNIT,1003)
     &                 (CONLEV,CHLIM(I),CHSTEP(:LENOCC(CHSTEP)),
     &                  ULIM(I),TEXT(IC)(:LENOCC(TEXT(IC))),I=1,2)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
*
      IF(KUMAC.AND.LUNIT.GT.0) THEN 
        IF(NC.LT.LSTA) THEN 
          DO I=1,MNAMES
            NAME=NAMES(I)//SUFFIX(LSTA)
            WRITE(LUNIT,6006) (NAME(:LENOCC(NAME)),J=1,2)
          ENDDO
        ENDIF
        WRITE(LUNIT,6002) 'exitm; endif'
      ENDIF
      RETURN
*
      ENTRY KUMAC_DMS_LIMIT(KUM)
      KUMAC=KUM
 6001 FORMAT('if (([1].eq.'' '').OR.([1].eq.''',A,''')).AND. _',/,
     &       '   (([2].eq.'' '').OR.([2].eq.''',A,''')) then')
 6002 FORMAT(2X,A,1X,A)
 6003 FORMAT(2X,'global/create ',A,' ''',A,'''')
 6004 FORMAT(2X,'vector/create ',A,'(',I3.3,') r ',
     &       (T32,5(F7.3,1X),:,'_'))
 6014 FORMAT('* vector/create ',A,'(',I3.3,') r ',
     &       (T32,F7.3,1X,'| value computed by SET_DMS_LIMIT'))
 6024 FORMAT(2X,'vector/create ',A,'(',I3.3,') r ',
     &       (T32,F7.3,1X,'| user-provided value'))
 6005 FORMAT(2X,'vector/create ',A,'(2,',I1,') r ',
     &       (T32,2(F7.3,1X),:,'_'))
 6006 FORMAT(2X,'if $VEXIST(',A,') then; vector/delete ',A,'; endif')
 6007 FORMAT(2X,'sigma ',A,' = ',F7.3)
      RETURN
*
*** Replace entry with subroutine (Jan 17, 2010)
*
***   ENTRY SET_DMS_LIMIT_USER(XUSER,USER)
*
      END
*
************************************************************************
*
      SUBROUTINE SET_DMS_LIMIT_USER(XUSER,USER)
*
      REAL USER(3)
      LOGICAL XUSER
      REAL UVAL(3,2)
      LOGICAL XVAL
      COMMON/UVAL_CM/UVAL,XVAL
*
*** End replace (Jan 17, 2010)
*
      DO I=1,3
        UVAL(I,1)=-1.
        UVAL(I,2)=-1.
      ENDDO
      XVAL=XUSER
      IF(XVAL) THEN
        UVAL(1,1)=-1.  ! lim syst incl
        UVAL(2,1)=-1.  ! lim stat only
        UVAL(3,1)=-1.  ! lim
        UVAL(1,2)=USER(1) ! sens syst incl
        UVAL(2,2)=USER(2) ! sens stat only
        UVAL(3,2)=USER(3) ! sens
      ENDIF
      END
*
************************************************************************
*
      SUBROUTINE FIND_ZERO(NSTEP,X,F,MZERO,NZERO,ZERO,DERIV)
*     ========================================================
*
*     Input:  NSTEP  = number of steps
*     ------  X      = array of X values (i.e. step values)
*             F      = array of values at the NSTEP values of X
*             MZERO  = dimension of arrays ZERO and DERIV
*
*     Output: NZERO  = number of zeroes of function F of variable X
*     ------- ZERO   = array of zeroes
*             DERIV  = array of the signs of the derivatives of F at the zeroes
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NSTEP,MZERO,NZERO,DERIV(MZERO)
      DOUBLE PRECISION X(NSTEP),F(NSTEP)
      REAL ZERO(MZERO)
*
*     Externals
*
      
*
*     Local variables
*
      INTEGER I,IDER
*
      NZERO=0
      IF(NSTEP.LE.0) RETURN
*
*zzz  PRINT * ,'FIND_ZERO:'
*zzz  DO I=1,NSTEP
*zzz    PRINT * ,'DMS,F = ',X(I),F(I)
*zzz  ENDDO
*
      IDER=NINT(DSIGN(1.D0,-F(1)))
      DO I=1,NSTEP
        IF(IDER.EQ.NINT(DSIGN(1.D0,F(I)))) THEN
          NZERO=NZERO+1
          IF(NZERO.LE.MZERO) THEN 
            ZERO(NZERO)=(F(I-1)*X(I)-F(I)*X(I-1))/
     &                  (F(I-1)     -F(I)       )
            DERIV(NZERO)=IDER
          ENDIF
          IDER=-IDER
        ENDIF
      ENDDO
*zzz      WRITE(*,*) 'FIND_ZERO: I,ZERO(I),DERIV(I)'
*zzz      DO I=1,NZERO
*zzz        WRITE(*,*) I,ZERO(I),DERIV(I)
*zzz      ENDDO
      END
